#pragma once

#include "nmopt/contract/complementarity.hpp"
#include "nmopt/contract/quadratic_kkt.hpp"
#include "nmopt/contract/quadratic_kkt_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::contract
{
  template <typename Backend>
  inline bool
  same_block_shape(const LayoutPtr &left, const LayoutPtr &right)
  {
    if (left->n_blocks() != right->n_blocks())
      return false;
    for (std::size_t block = 0; block < left->n_blocks(); ++block)
      if (left->dimension(block) != right->dimension(block))
        return false;
    return true;
  }

  template <typename Backend>
  inline std::vector<typename Backend::Vector>
  copy_blocks(const BlockValuesT<Backend> &values)
  {
    std::vector<typename Backend::Vector> blocks;
    blocks.reserve(values.n_blocks());
    for (std::size_t block = 0; block < values.n_blocks(); ++block)
      blocks.push_back(values.block(block));
    return blocks;
  }

  template <typename Backend>
  inline std::string
  unique_appended_space_id(const LayoutPtr &layout, std::string stem)
  {
    std::string candidate = std::move(stem);
    std::size_t suffix = 0;
    while (true)
      {
        bool collision = false;
        for (std::size_t block = 0; block < layout->n_blocks(); ++block)
          if (layout->space(block).value == candidate)
            {
              collision = true;
              break;
            }
        if (!collision)
          return candidate;
        ++suffix;
        candidate = "pdas_active_box_" + std::to_string(suffix);
      }
  }

  template <typename Backend>
  inline LayoutPtr
  append_layout_block(const LayoutPtr &base,
                      const std::string &label,
                      const std::string &space_stem,
                      const std::size_t dimension)
  {
    require(dimension > 0, "PDAS active layout needs a positive dimension");
    std::vector<SpaceId> spaces;
    std::vector<std::size_t> dimensions;
    spaces.reserve(base->n_blocks() + 1);
    dimensions.reserve(base->n_blocks() + 1);
    for (std::size_t block = 0; block < base->n_blocks(); ++block)
      {
        spaces.push_back(base->space(block));
        dimensions.push_back(base->dimension(block));
      }
    spaces.push_back(
      SpaceId{unique_appended_space_id<Backend>(base, space_stem)});
    dimensions.push_back(dimension);
    return std::make_shared<const BlockLayout>(
      label, std::move(spaces), std::move(dimensions));
  }

  template <typename Backend>
  class ActiveSetKKTSubproblemT final
  {
  public:
    using Product = EqualityConstrainedQuadraticKKTProductT<Backend>;
    using Complementarity = BoxComplementarityT<Backend>;
    using Selection = ActiveSetSelectionT<Backend>;
    using Primal = typename Product::Primal;
    using Covector = typename Product::Covector;
    using Point = typename Product::Point;
    using Vector = typename Backend::Vector;

    ActiveSetKKTSubproblemT(
      const Product &       base,
      const Complementarity &complementarity,
      const Selection &     selection,
      const std::size_t     control_block,
      QuadraticKKTAssumptions active_set_assumptions)
      : base_(base)
      , selection_(selection)
      , control_block_(control_block)
    {
      validate_base(complementarity);
      active_values_ = make_active_values(complementarity);
      if (selection_.active_size() > 0)
        active_product_.emplace(
          make_active_product(std::move(active_set_assumptions)));
    }

    const Product &
    product() const
    {
      return active_product_ ? *active_product_ : base_;
    }

    bool
    has_active_constraints() const
    {
      return selection_.active_size() > 0;
    }

    const Selection &
    selection() const
    {
      return selection_;
    }

    const Vector &
    active_values() const
    {
      return active_values_;
    }

    Point
    to_base_point(const Point &point) const
    {
      require(point.primal.layout()->compatible_with(
                *base_.layout().primal),
              "PDAS KKT solution has an incompatible primal layout");
      if (!active_product_)
        {
          require(point.multiplier.layout()->compatible_with(
                    *base_.layout().multiplier),
                  "PDAS KKT solution has an incompatible multiplier layout");
          return point;
        }

      require(point.multiplier.layout()->compatible_with(
                *active_product_->layout().multiplier),
              "PDAS active KKT solution has an incompatible multiplier layout");
      std::vector<Vector> multiplier_blocks;
      multiplier_blocks.reserve(base_.layout().multiplier->n_blocks());
      for (std::size_t block = 0;
           block < base_.layout().multiplier->n_blocks();
           ++block)
        multiplier_blocks.push_back(point.multiplier.block(block));
      return Point{
        point.primal,
        Primal(base_.layout().multiplier, std::move(multiplier_blocks))};
    }

  private:
    void
    validate_base(const Complementarity &complementarity) const
    {
      require(control_block_ < base_.layout().primal->n_blocks(),
              "PDAS control block exceeds the KKT primal layout");
      require(control_block_ < base_.layout().stationarity->n_blocks(),
              "PDAS control block exceeds the KKT stationarity layout");
      require(same_block_shape<Backend>(base_.layout().primal,
                                        base_.layout().stationarity),
              "PDAS KKT primal and stationarity block shapes differ");
      require(base_.layout().primal->dimension(control_block_) ==
                complementarity.layout()->dimension(0),
              "PDAS control block has an incompatible complementarity layout");
      require(selection_.layout()->compatible_with(
                *complementarity.layout()),
              "PDAS selection has an incompatible complementarity layout");
    }

    Vector
    make_active_values(const Complementarity &complementarity) const
    {
      Vector values = Backend::zeros(selection_.active_size());
      for (std::size_t selected = 0;
           selected < selection_.active_indices().size();
           ++selected)
        {
          const std::size_t index = selection_.active_indices()[selected];
          const double value =
            selection_.activities()[index] == BoxActivity::upper
              ? Backend::value(complementarity.bounds().upper().block(0), index)
              : Backend::value(complementarity.bounds().lower().block(0), index);
          Backend::set_value(values, selected, value);
        }
      return values;
    }

    static std::vector<Vector>
    negate_blocks(const BlockValuesT<Backend> &values)
    {
      std::vector<Vector> blocks = copy_blocks(values);
      for (auto &block : blocks)
        Backend::scale(block, -1.0);
      return blocks;
    }

    Product
    make_active_product(QuadraticKKTAssumptions assumptions) const
    {
      const LayoutPtr equality_layout = append_layout_block<Backend>(
        base_.layout().equality,
        "pdas_equality",
        "pdas_active_box",
        selection_.active_size());
      const LayoutPtr multiplier_layout = append_layout_block<Backend>(
        base_.layout().multiplier,
        "pdas_multiplier",
        "pdas_box_multiplier",
        selection_.active_size());
      const typename Product::Layout layout(base_.layout().primal,
                                            multiplier_layout,
                                            base_.layout().adjoint,
                                            base_.layout().stationarity,
                                            equality_layout);

      const auto base_zero = Point{
        Primal::zeros(base_.layout().primal),
        Primal::zeros(base_.layout().multiplier)};
      const auto base_zero_residual = base_.residual(base_zero);
      const Vector active_values = active_values_;
      const Selection selection = selection_;
      const std::size_t control_block = control_block_;

      const auto quadratic_action = [this](const Primal &primal) {
        return base_.apply_q(primal);
      };

      const auto equality_action = [this,
                                    equality_layout,
                                    selection,
                                    control_block](const Primal &primal) {
        const Covector base_value = base_.apply_d(primal);
        std::vector<Vector> blocks = copy_blocks(base_value);
        blocks.push_back(selection.restrict_active(
          primal.block(control_block)));
        return Covector(equality_layout, std::move(blocks));
      };

      const auto multiplier_action = [this,
                                     multiplier_layout,
                                     selection,
                                     control_block](const Primal &multiplier) {
        const std::size_t base_blocks =
          base_.layout().multiplier->n_blocks();
        std::vector<Vector> base_multiplier_blocks;
        base_multiplier_blocks.reserve(base_blocks);
        for (std::size_t block = 0; block < base_blocks; ++block)
          base_multiplier_blocks.push_back(multiplier.block(block));
        const Primal base_multiplier(base_.layout().multiplier,
                                     std::move(base_multiplier_blocks));
        const Covector base_value = base_.apply_d_transpose(base_multiplier);
        std::vector<Vector> blocks = copy_blocks(base_value);
        Vector &control = blocks.at(control_block);
        const Vector &active = multiplier.block(base_blocks);
        for (std::size_t selected = 0;
             selected < selection.active_indices().size();
             ++selected)
          Backend::set_value(control,
                             selection.active_indices()[selected],
                             Backend::value(active, selected));
        return Covector(base_.layout().stationarity, std::move(blocks));
      };

      const auto transpose_action = [this,
                                     equality_layout,
                                     multiplier_layout,
                                     selection,
                                     control_block](const typename Product::Seed &seed) {
        const std::size_t base_equality_blocks =
          base_.layout().equality->n_blocks();
        std::vector<Vector> base_equality_values;
        base_equality_values.reserve(base_equality_blocks);
        for (std::size_t block = 0; block < base_equality_blocks; ++block)
          base_equality_values.push_back(seed.equality.block(block));
        const Primal base_stationarity(
          base_.layout().stationarity,
          copy_blocks(seed.stationarity));
        const Primal base_equality(base_.layout().equality,
                                   std::move(base_equality_values));
        const auto base_result = base_.apply_kkt_transpose(
          typename Product::Seed{base_stationarity, base_equality});

        std::vector<Vector> primal_blocks = copy_blocks(base_result.primal);
        Vector &control = primal_blocks.at(control_block);
        const Vector &active_seed = seed.equality.block(base_equality_blocks);
        for (std::size_t selected = 0;
             selected < selection.active_indices().size();
             ++selected)
          {
            const std::size_t index = selection.active_indices()[selected];
            const double value = Backend::value(control, index) +
                                 Backend::value(active_seed, selected);
            Backend::set_value(control, index, value);
          }

        std::vector<Vector> multiplier_blocks =
          copy_blocks(base_result.multiplier);
        multiplier_blocks.push_back(selection.restrict_active(
          seed.stationarity.block(control_block)));
        return typename Product::TransposeResult{
          Covector(base_.layout().primal, std::move(primal_blocks)),
          Covector(multiplier_layout, std::move(multiplier_blocks))};
      };

      std::vector<Vector> stationarity_rhs =
        negate_blocks(base_zero_residual.stationarity);
      std::vector<Vector> equality_rhs =
        negate_blocks(base_zero_residual.equality);
      equality_rhs.push_back(active_values);

      const auto multiplier_conversion =
        typename Product::MultiplierConversion{
          "base KKT multiplier conversion with separate PDAS box multiplier",
          [this](const Primal &multiplier) {
            const std::size_t base_blocks =
              base_.layout().multiplier->n_blocks();
            std::vector<Vector> blocks;
            blocks.reserve(base_blocks);
            for (std::size_t block = 0; block < base_blocks; ++block)
              blocks.push_back(multiplier.block(block));
            return base_.multiplier_to_adjoint(
              Primal(base_.layout().multiplier, std::move(blocks)));
          },
          [this, multiplier_layout](const Primal &adjoint) {
            const Primal base_multiplier =
              base_.adjoint_to_multiplier(adjoint);
            std::vector<Vector> blocks = copy_blocks(base_multiplier);
            blocks.push_back(Backend::zeros(selection_.active_size()));
            return Primal(multiplier_layout, std::move(blocks));
          }};

      return Product(
        layout,
        quadratic_action,
        equality_action,
        multiplier_action,
        transpose_action,
        Covector(base_.layout().stationarity, std::move(stationarity_rhs)),
        Covector(equality_layout, std::move(equality_rhs)),
        multiplier_conversion,
        std::move(assumptions),
        base_.symmetry());
    }

    const Product &              base_;
    Selection                    selection_;
    std::size_t                  control_block_;
    Vector                       active_values_;
    std::optional<Product>       active_product_;
  };

  enum class PDASStoppingReason
  {
    converged,
    maximum_iterations,
    kkt_solve_failed
  };

  struct PDASPolicy
  {
    std::size_t maximum_iterations = 50;
    double      classification_parameter = 1.0;
    double      primal_feasibility_tolerance = 1e-10;
    double      dual_feasibility_tolerance = 1e-10;
    double      complementarity_tolerance = 1e-10;
    double      stationarity_tolerance = 1e-10;
    double      equality_tolerance = 1e-10;
    QuadraticKKTAssumptions active_set_assumptions;
  };

  inline bool
  valid(const PDASPolicy &policy)
  {
    return policy.maximum_iterations > 0 &&
           std::isfinite(policy.classification_parameter) &&
           policy.classification_parameter > 0.0 &&
           std::isfinite(policy.primal_feasibility_tolerance) &&
           policy.primal_feasibility_tolerance > 0.0 &&
           std::isfinite(policy.dual_feasibility_tolerance) &&
           policy.dual_feasibility_tolerance > 0.0 &&
           std::isfinite(policy.complementarity_tolerance) &&
           policy.complementarity_tolerance > 0.0 &&
           std::isfinite(policy.stationarity_tolerance) &&
           policy.stationarity_tolerance > 0.0 &&
           std::isfinite(policy.equality_tolerance) &&
           policy.equality_tolerance > 0.0;
  }

  template <typename Backend>
  inline double
  pdas_block_norm(const BlockValuesT<Backend> &values)
  {
    double squared_norm = 0.0;
    for (std::size_t block = 0; block < values.n_blocks(); ++block)
      squared_norm += Backend::dot(values.block(block), values.block(block));
    return std::sqrt(squared_norm);
  }

  template <typename Backend>
  struct PDASIterationReportT
  {
    using Selection = ActiveSetSelectionT<Backend>;

    std::size_t             iteration = 0;
    Selection               selection;
    std::size_t             active_set_changes = 0;
    bool                    active_set_stable = false;
    double                  primal_violation = 0.0;
    double                  dual_violation = 0.0;
    double                  complementarity_residual = 0.0;
    double                  stationarity_residual = 0.0;
    double                  equality_residual = 0.0;
    bool                    primal_feasible = false;
    bool                    dual_feasible = false;
    bool                    complementarity_converged = false;
    bool                    kkt_residuals_converged = false;
    QuadraticKKTSolveReport kkt_solve;
  };

  template <typename Backend>
  struct PDASSolveResultT
  {
    using Product = EqualityConstrainedQuadraticKKTProductT<Backend>;
    using Covector = CovectorBlockT<Backend>;
    using IterationReport = PDASIterationReportT<Backend>;

    typename Product::Point solution;
    Covector                box_multiplier;
    std::vector<IterationReport> iterations;
    PDASStoppingReason      stopping_reason;

    bool
    converged() const
    {
      return stopping_reason == PDASStoppingReason::converged;
    }
  };

  template <typename Backend>
  class PDASSolverT final
  {
  public:
    using Product = EqualityConstrainedQuadraticKKTProductT<Backend>;
    using Complementarity = BoxComplementarityT<Backend>;
    using Primal = PrimalBlockT<Backend>;
    using Covector = CovectorBlockT<Backend>;
    using Point = typename Product::Point;
    using SolveResult = QuadraticKKTSolveResultT<Backend>;
    using IterationReport = PDASIterationReportT<Backend>;
    using Result = PDASSolveResultT<Backend>;
    using SolveAction = std::function<SolveResult(const Product &)>;

    PDASSolverT(const Product &       product,
                const Complementarity &complementarity,
                const std::size_t     control_block,
                SolveAction           solve_action)
      : product_(product)
      , complementarity_(complementarity)
      , control_block_(control_block)
      , solve_action_(std::move(solve_action))
    {
      require(static_cast<bool>(solve_action_),
              "PDAS solver needs a KKT solve action");
      require(control_block_ < product_.layout().primal->n_blocks(),
              "PDAS control block exceeds the KKT primal layout");
      require(control_block_ < product_.layout().stationarity->n_blocks(),
              "PDAS control block exceeds the KKT stationarity layout");
      require(same_block_shape<Backend>(product_.layout().primal,
                                        product_.layout().stationarity),
              "PDAS KKT primal and stationarity block shapes differ");
      require(product_.layout().primal->dimension(control_block_) ==
                complementarity_.layout()->dimension(0),
              "PDAS control block has an incompatible complementarity layout");
    }

    Result
    solve(const Point &initial_point,
          const Covector &initial_box_multiplier,
          const PDASPolicy &policy) const
    {
      require(valid(policy), "PDAS policy is invalid");
      require(initial_point.primal.layout()->compatible_with(
                *product_.layout().primal),
              "PDAS initial point has an incompatible primal layout");
      require(initial_point.multiplier.layout()->compatible_with(
                *product_.layout().multiplier),
              "PDAS initial point has an incompatible multiplier layout");
      require_complementarity_layout(
        initial_box_multiplier,
        complementarity_.layout(),
        "PDAS initial box multiplier");
      const Primal initial_control(
        complementarity_.layout(),
        {initial_point.primal.block(control_block_)});
      require(complementarity_.bounds().is_feasible(initial_control),
              "PDAS initial control must be feasible");

      Point current = initial_point;
      Covector current_box_multiplier = initial_box_multiplier;
      auto selection = complementarity_.classify(
        Primal(complementarity_.layout(),
               {current.primal.block(control_block_)}),
        current_box_multiplier,
        policy.classification_parameter);
      std::vector<IterationReport> reports;

      for (std::size_t iteration = 0;
           iteration < policy.maximum_iterations;
           ++iteration)
        {
          ActiveSetKKTSubproblemT<Backend> subproblem(
            product_,
            complementarity_,
            selection,
            control_block_,
            policy.active_set_assumptions);
          const SolveResult kkt_result = solve_action_(subproblem.product());
          const Point base_solution =
            subproblem.to_base_point(kkt_result.solution);
          const auto residual = product_.residual(base_solution);
          const Covector box_multiplier = make_box_multiplier(
            residual.stationarity,
            selection);
          const auto next_selection = complementarity_.classify(
            Primal(complementarity_.layout(),
                   {base_solution.primal.block(control_block_)}),
            box_multiplier,
            policy.classification_parameter);
          IterationReport report = make_report(
            iteration,
            selection,
            next_selection,
            base_solution.primal.block(control_block_),
            box_multiplier,
            residual,
            kkt_result.report,
            policy);
          reports.push_back(report);

          const bool converged =
            kkt_result.report.converged() &&
            report.active_set_stable &&
            report.kkt_residuals_converged;
          if (converged)
            return Result{base_solution,
                          box_multiplier,
                          std::move(reports),
                          PDASStoppingReason::converged};
          if (!kkt_result.report.converged())
            return Result{base_solution,
                          box_multiplier,
                          std::move(reports),
                          PDASStoppingReason::kkt_solve_failed};

          current = base_solution;
          current_box_multiplier = box_multiplier;
          selection = next_selection;
        }

      return Result{current,
                    current_box_multiplier,
                    std::move(reports),
                    PDASStoppingReason::maximum_iterations};
    }

  private:
    using Vector = typename Backend::Vector;
    using Residual = typename Product::Residual;

    Covector
    make_box_multiplier(const Covector &stationarity,
                        const ActiveSetSelectionT<Backend> &selection) const
    {
      require(stationarity.n_blocks() > control_block_,
              "PDAS stationarity has no control block");
      Vector values = Backend::zeros(complementarity_.layout()->dimension(0));
      const Vector &raw = stationarity.block(control_block_);
      for (std::size_t selected = 0;
           selected < selection.active_indices().size();
           ++selected)
        {
          const std::size_t index = selection.active_indices()[selected];
          Backend::set_value(values,
                             index,
                             -Backend::value(raw, index));
        }
      return Covector(complementarity_.layout(), {std::move(values)});
    }

    IterationReport
    make_report(const std::size_t                   iteration,
                const ActiveSetSelectionT<Backend> &selection,
                const ActiveSetSelectionT<Backend> &next_selection,
                const Vector &                      point_control,
                const Covector &                    box_multiplier,
                const Residual &                    residual,
                const QuadraticKKTSolveReport &     kkt_solve,
                const PDASPolicy &                  policy) const
    {
      require(residual.stationarity.n_blocks() > control_block_,
              "PDAS residual stationarity has no control block");
      const Primal representative =
        complementarity_.multiplier_to_primal(box_multiplier);
      double primal_violation = 0.0;
      double dual_violation = 0.0;
      double complementarity_residual = 0.0;
      for (std::size_t index = 0;
           index < complementarity_.layout()->dimension(0);
           ++index)
        {
          const double value = Backend::value(point_control, index);
          const double multiplier =
            Backend::value(representative.block(0), index);
          const double lower = Backend::value(
            complementarity_.bounds().lower().block(0), index);
          const double upper = Backend::value(
            complementarity_.bounds().upper().block(0), index);
          primal_violation = std::max(
            primal_violation,
            std::max(lower - value, value - upper));
          if (value <= lower)
            dual_violation = std::max(dual_violation, multiplier);
          else if (value >= upper)
            dual_violation = std::max(dual_violation, -multiplier);
          else
            dual_violation = std::max(dual_violation, std::abs(multiplier));
          complementarity_residual = std::max(
            complementarity_residual,
            std::max(std::abs(std::max(multiplier, 0.0) * (upper - value)),
                     std::abs(std::min(multiplier, 0.0) * (value - lower))));
        }

      std::vector<Vector> stationarity_blocks =
        copy_blocks(residual.stationarity);
      Vector &control_stationarity = stationarity_blocks.at(control_block_);
      Backend::add_scaled(control_stationarity,
                          1.0,
                          box_multiplier.block(0));
      const Covector constrained_stationarity(
        residual.stationarity.layout(), std::move(stationarity_blocks));
      const double stationarity_residual =
        pdas_block_norm(constrained_stationarity);
      const double equality_residual = pdas_block_norm(residual.equality);
      const bool primal_feasible =
        primal_violation <= policy.primal_feasibility_tolerance;
      const bool dual_feasible =
        dual_violation <= policy.dual_feasibility_tolerance;
      const bool complementarity_converged =
        complementarity_residual <= policy.complementarity_tolerance;
      const bool kkt_residuals_converged =
        primal_feasible && dual_feasible && complementarity_converged &&
        stationarity_residual <= policy.stationarity_tolerance &&
        equality_residual <= policy.equality_tolerance;

      std::size_t active_set_changes = 0;
      for (std::size_t index = 0; index < selection.activities().size(); ++index)
        if (selection.activities()[index] != next_selection.activities()[index])
          ++active_set_changes;

      return IterationReport{iteration,
                             selection,
                             active_set_changes,
                             selection == next_selection,
                             primal_violation,
                             dual_violation,
                             complementarity_residual,
                             stationarity_residual,
                             equality_residual,
                             primal_feasible,
                             dual_feasible,
                             complementarity_converged,
                             kkt_residuals_converged,
                             kkt_solve};
    }

    const Product &        product_;
    const Complementarity &complementarity_;
    std::size_t            control_block_;
    SolveAction            solve_action_;
  };

  using ActiveSetKKTSubproblem = ActiveSetKKTSubproblemT<DenseBackend>;
  using PDASSolver = PDASSolverT<DenseBackend>;
} // namespace nmopt::contract
