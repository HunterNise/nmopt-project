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
      : base_(std::make_shared<const Product>(base))
      , selection_(selection)
      , control_block_(control_block)
    {
      validate_base(complementarity);
      active_values_ = make_active_values(complementarity);
      if (selection_.active_size() > 0)
        {
          initialize_restriction_layouts();
          active_product_.emplace(
            make_active_product(std::move(active_set_assumptions)));
        }
    }

    const Product &
    product() const
    {
      return active_product_ ? *active_product_ : *base_;
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
      if (!active_product_)
        {
          require(point.primal.layout()->compatible_with(
                    *base_->layout().primal),
                  "PDAS KKT solution has an incompatible primal layout");
          require(point.multiplier.layout()->compatible_with(
                    *base_->layout().multiplier),
                  "PDAS KKT solution has an incompatible multiplier layout");
          return point;
        }

      require(point.primal.layout()->compatible_with(
                *restricted_primal_layout_),
              "PDAS active KKT solution has an incompatible primal layout");
      require(point.multiplier.layout()->compatible_with(
                *base_->layout().multiplier),
              "PDAS active KKT solution has an incompatible multiplier layout");
      return Point{
        expand_primal(point.primal, active_values_),
        point.multiplier};
    }

  private:
    void
    validate_base(const Complementarity &complementarity) const
    {
      require(control_block_ < base_->layout().primal->n_blocks(),
              "PDAS control block exceeds the KKT primal layout");
      require(control_block_ < base_->layout().stationarity->n_blocks(),
              "PDAS control block exceeds the KKT stationarity layout");
      require(same_block_shape<Backend>(base_->layout().primal,
                                        base_->layout().stationarity),
              "PDAS KKT primal and stationarity block shapes differ");
      require(base_->layout().primal->dimension(control_block_) ==
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

    void
    initialize_restriction_layouts()
    {
      reduced_block_for_base_.clear();
      reduced_block_for_base_.reserve(base_->layout().primal->n_blocks());
      std::size_t next_block = 0;
      for (std::size_t block = 0;
           block < base_->layout().primal->n_blocks();
           ++block)
        {
          if (block == control_block_ && selection_.free_size() == 0)
            {
              reduced_block_for_base_.push_back(std::nullopt);
              continue;
            }
          reduced_block_for_base_.push_back(next_block++);
        }

      restricted_primal_layout_ = make_restricted_layout(
        base_->layout().primal, "pdas_restricted_primal");
      restricted_stationarity_layout_ = make_restricted_layout(
        base_->layout().stationarity, "pdas_restricted_stationarity");
    }

    LayoutPtr
    make_restricted_layout(const LayoutPtr &base,
                           const std::string &label) const
    {
      std::vector<SpaceId> spaces;
      std::vector<std::size_t> dimensions;
      spaces.reserve(base->n_blocks());
      dimensions.reserve(base->n_blocks());

      std::string free_control_space = "pdas_free_control";
      std::size_t suffix = 0;
      while (true)
        {
          bool collision = false;
          for (std::size_t block = 0; block < base->n_blocks(); ++block)
            if (block != control_block_ &&
                base->space(block).value == free_control_space)
              {
                collision = true;
                break;
              }
          if (!collision)
            break;
          free_control_space = "pdas_free_control_" +
                               std::to_string(++suffix);
        }

      for (std::size_t block = 0; block < base->n_blocks(); ++block)
        {
          if (block == control_block_ && selection_.free_size() == 0)
            continue;
          if (block == control_block_)
            {
              spaces.push_back(SpaceId{free_control_space});
              dimensions.push_back(selection_.free_size());
            }
          else
            {
              spaces.push_back(base->space(block));
              dimensions.push_back(base->dimension(block));
            }
        }

      require(!spaces.empty(),
              "PDAS restricted KKT product needs a free or state block");
      return std::make_shared<const BlockLayout>(
        label, std::move(spaces), std::move(dimensions));
    }

    QuadraticKKTBlockPairing
    make_restricted_pairing(
      const QuadraticKKTBlockPairing &base_pairing) const
    {
      QuadraticKKTBlockPairing pairing = base_pairing;
      pairing.id += "_restricted";
      pairing.domain_blocks.clear();
      pairing.range_blocks.clear();
      pairing.pairing_ids.clear();

      for (std::size_t pair = 0;
           pair < base_pairing.domain_blocks.size();
           ++pair)
        {
          const std::size_t domain_block = base_pairing.domain_blocks[pair];
          const std::size_t range_block = base_pairing.range_blocks[pair];
          if (selection_.free_size() == 0 &&
              (domain_block == control_block_ ||
               range_block == control_block_))
            continue;
          require(domain_block < reduced_block_for_base_.size() &&
                    range_block < reduced_block_for_base_.size() &&
                    reduced_block_for_base_[domain_block].has_value() &&
                    reduced_block_for_base_[range_block].has_value(),
                  "PDAS restricted KKT pairing lost a free block");
          pairing.domain_blocks.push_back(
            *reduced_block_for_base_[domain_block]);
          pairing.range_blocks.push_back(
            *reduced_block_for_base_[range_block]);
          pairing.pairing_ids.push_back(base_pairing.pairing_ids[pair]);
        }
      return pairing;
    }

    Primal
    expand_primal(const Primal &restricted,
                  const Vector &active_values) const
    {
      require(restricted.layout()->compatible_with(
                *restricted_primal_layout_),
              "PDAS restricted primal has an incompatible layout");
      require(Backend::size(active_values) == selection_.active_size(),
              "PDAS active values have the wrong size");

      std::vector<Vector> blocks;
      blocks.reserve(base_->layout().primal->n_blocks());
      for (std::size_t block = 0;
           block < base_->layout().primal->n_blocks();
           ++block)
        {
          if (block == control_block_)
            {
              Vector free_values = Backend::zeros(selection_.free_size());
              if (selection_.free_size() > 0)
                free_values = restricted.block(
                  *reduced_block_for_base_[block]);
              blocks.push_back(
                selection_.prolong(free_values, active_values));
            }
          else
            blocks.push_back(
              restricted.block(*reduced_block_for_base_[block]));
        }
      return Primal(base_->layout().primal, std::move(blocks));
    }

    Covector
    restrict_stationarity(const Covector &value) const
    {
      require(value.layout()->compatible_with(
                *base_->layout().stationarity),
              "PDAS base stationarity has an incompatible layout");
      std::vector<Vector> blocks;
      blocks.reserve(restricted_stationarity_layout_->n_blocks());
      for (std::size_t block = 0;
           block < base_->layout().stationarity->n_blocks();
           ++block)
        {
          if (block == control_block_)
            {
              if (selection_.free_size() > 0)
                blocks.push_back(selection_.restrict_free(
                  value.block(block)));
            }
          else
            blocks.push_back(value.block(block));
        }
      return Covector(restricted_stationarity_layout_, std::move(blocks));
    }

    Covector
    restrict_primal(const Covector &value) const
    {
      require(value.layout()->compatible_with(*base_->layout().primal),
              "PDAS base primal covector has an incompatible layout");
      std::vector<Vector> blocks;
      blocks.reserve(restricted_primal_layout_->n_blocks());
      for (std::size_t block = 0;
           block < base_->layout().primal->n_blocks();
           ++block)
        {
          if (block == control_block_)
            {
              if (selection_.free_size() > 0)
                blocks.push_back(selection_.restrict_free(
                  value.block(block)));
            }
          else
            blocks.push_back(value.block(block));
        }
      return Covector(restricted_primal_layout_, std::move(blocks));
    }

    struct ActiveProductState
    {
      std::shared_ptr<const Product>              base;
      Selection                                  selection;
      std::size_t                                control_block;
      Vector                                     active_values;
      std::vector<std::optional<std::size_t>>   reduced_block_for_base;
      LayoutPtr                                  restricted_primal_layout;
      LayoutPtr                                  restricted_stationarity_layout;

      Primal
      expand_primal(const Primal &restricted,
                    const Vector &active) const
      {
        require(restricted.layout()->compatible_with(
                  *restricted_primal_layout),
                "PDAS restricted primal has an incompatible layout");
        require(Backend::size(active) == selection.active_size(),
                "PDAS active values have the wrong size");

        std::vector<Vector> blocks;
        blocks.reserve(base->layout().primal->n_blocks());
        for (std::size_t block = 0;
             block < base->layout().primal->n_blocks();
             ++block)
          {
            if (block == control_block)
              {
                Vector free_values = Backend::zeros(selection.free_size());
                if (selection.free_size() > 0)
                  free_values = restricted.block(
                    *reduced_block_for_base[block]);
                blocks.push_back(selection.prolong(free_values, active));
              }
            else
              blocks.push_back(
                restricted.block(*reduced_block_for_base[block]));
          }
        return Primal(base->layout().primal, std::move(blocks));
      }

      Covector
      restrict_stationarity(const Covector &value) const
      {
        require(value.layout()->compatible_with(
                  *base->layout().stationarity),
                "PDAS base stationarity has an incompatible layout");
        std::vector<Vector> blocks;
        blocks.reserve(restricted_stationarity_layout->n_blocks());
        for (std::size_t block = 0;
             block < base->layout().stationarity->n_blocks();
             ++block)
          {
            if (block == control_block)
              {
                if (selection.free_size() > 0)
                  blocks.push_back(selection.restrict_free(value.block(block)));
              }
            else
              blocks.push_back(value.block(block));
          }
        return Covector(restricted_stationarity_layout, std::move(blocks));
      }

      Covector
      restrict_primal(const Covector &value) const
      {
        require(value.layout()->compatible_with(*base->layout().primal),
                "PDAS base primal covector has an incompatible layout");
        std::vector<Vector> blocks;
        blocks.reserve(restricted_primal_layout->n_blocks());
        for (std::size_t block = 0;
             block < base->layout().primal->n_blocks();
             ++block)
          {
            if (block == control_block)
              {
                if (selection.free_size() > 0)
                  blocks.push_back(selection.restrict_free(value.block(block)));
              }
            else
              blocks.push_back(value.block(block));
          }
        return Covector(restricted_primal_layout, std::move(blocks));
      }
    };

    Product
    make_active_product(QuadraticKKTAssumptions assumptions) const
    {
      const typename Product::Layout layout(restricted_primal_layout_,
                                            base_->layout().multiplier,
                                            base_->layout().adjoint,
                                            restricted_stationarity_layout_,
                                            base_->layout().equality,
                                            make_restricted_pairing(
                                              base_->layout().primal_stationarity_pairing),
                                            base_->layout().multiplier_equality_pairing);

      const auto base_active_point = Point{
        expand_primal(Primal::zeros(restricted_primal_layout_), active_values_),
        Primal::zeros(base_->layout().multiplier)};
      const auto base_zero_residual = base_->residual(base_active_point);

      const auto state = std::make_shared<const ActiveProductState>(
        ActiveProductState{base_,
                           selection_,
                           control_block_,
                           active_values_,
                           reduced_block_for_base_,
                           restricted_primal_layout_,
                           restricted_stationarity_layout_});

      const auto quadratic_action = [state](const Primal &primal) {
        const Vector zero_active = Backend::zeros(state->selection.active_size());
        return state->restrict_stationarity(
          state->base->apply_q(state->expand_primal(primal, zero_active)));
      };

      const auto equality_action = [state](const Primal &primal) {
        const Vector zero_active = Backend::zeros(state->selection.active_size());
        return state->base->apply_d(
          state->expand_primal(primal, zero_active));
      };

      const auto multiplier_action = [state](const Primal &multiplier) {
        return state->restrict_stationarity(
          state->base->apply_d_transpose(multiplier));
      };

      const auto transpose_action = [state](const typename Product::Seed &seed) {
        const Vector zero_active = Backend::zeros(state->selection.active_size());
        const auto base_result = state->base->apply_kkt_transpose(
          typename Product::Seed{
            state->expand_primal(seed.stationarity, zero_active),
            seed.equality});
        return typename Product::TransposeResult{
          state->restrict_primal(base_result.primal),
          base_result.multiplier};
      };

      const Covector restricted_stationarity_residual =
        restrict_stationarity(base_zero_residual.stationarity);
      const std::vector<Vector> stationarity_rhs =
        negate_blocks(restricted_stationarity_residual);
      const std::vector<Vector> equality_rhs =
        negate_blocks(base_zero_residual.equality);

      const auto multiplier_conversion =
        typename Product::MultiplierConversion{
          "base KKT multiplier conversion on restricted free coordinates",
          [state](const Primal &multiplier) {
            return state->base->multiplier_to_adjoint(multiplier);
          },
          [state](const Primal &adjoint) {
            return state->base->adjoint_to_multiplier(adjoint);
          }};

      return Product(
        layout,
        quadratic_action,
        equality_action,
        multiplier_action,
        transpose_action,
        Covector(restricted_stationarity_layout_,
                 std::move(stationarity_rhs)),
        Covector(base_->layout().equality, std::move(equality_rhs)),
        multiplier_conversion,
        std::move(assumptions),
        base_->symmetry(),
        state);
    }

    std::shared_ptr<const Product> base_;
    Selection                    selection_;
    std::size_t                  control_block_;
    Vector                       active_values_;
    std::vector<std::optional<std::size_t>> reduced_block_for_base_;
    LayoutPtr                    restricted_primal_layout_;
    LayoutPtr                    restricted_stationarity_layout_;
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
      : product_(std::make_shared<const Product>(product))
      , complementarity_(std::make_shared<const Complementarity>(complementarity))
      , control_block_(control_block)
      , solve_action_(std::move(solve_action))
    {
      require(static_cast<bool>(solve_action_),
              "PDAS solver needs a KKT solve action");
      require(control_block_ < product_->layout().primal->n_blocks(),
              "PDAS control block exceeds the KKT primal layout");
      require(control_block_ < product_->layout().stationarity->n_blocks(),
              "PDAS control block exceeds the KKT stationarity layout");
      require(same_block_shape<Backend>(product_->layout().primal,
                                        product_->layout().stationarity),
              "PDAS KKT primal and stationarity block shapes differ");
      require(product_->layout().primal->dimension(control_block_) ==
                complementarity_->layout()->dimension(0),
              "PDAS control block has an incompatible complementarity layout");
    }

    Result
    solve(const Point &initial_point,
          const Covector &initial_box_multiplier,
          const PDASPolicy &policy) const
    {
      require(valid(policy), "PDAS policy is invalid");
      require(initial_point.primal.layout()->compatible_with(
                *product_->layout().primal),
              "PDAS initial point has an incompatible primal layout");
      require(initial_point.multiplier.layout()->compatible_with(
                *product_->layout().multiplier),
              "PDAS initial point has an incompatible multiplier layout");
      require_finite_block_values(initial_point.primal,
                                  "PDAS initial primal point");
      require_finite_block_values(initial_point.multiplier,
                                  "PDAS initial equality multiplier");
      require_complementarity_layout(
        initial_box_multiplier,
        complementarity_->layout(),
        "PDAS initial box multiplier");
      require_finite_block_values(initial_box_multiplier,
                                  "PDAS initial box multiplier");
      const Primal initial_control(
        complementarity_->layout(),
        {initial_point.primal.block(control_block_)});
      require(complementarity_->bounds().is_feasible(initial_control),
              "PDAS initial control must be feasible");

      Point current = initial_point;
      Covector current_box_multiplier = initial_box_multiplier;
      auto selection = complementarity_->classify(
        Primal(complementarity_->layout(),
               {current.primal.block(control_block_)}),
        current_box_multiplier,
        policy.classification_parameter);
      std::vector<IterationReport> reports;

      for (std::size_t iteration = 0;
           iteration < policy.maximum_iterations;
           ++iteration)
        {
          ActiveSetKKTSubproblemT<Backend> subproblem(
            *product_,
            *complementarity_,
            selection,
            control_block_,
            policy.active_set_assumptions);
          const SolveResult kkt_result = solve_action_(subproblem.product());
          require_finite_block_values(kkt_result.solution.primal,
                                      "PDAS KKT solution primal");
          require_finite_block_values(kkt_result.solution.multiplier,
                                      "PDAS KKT solution equality multiplier");
          require_finite_solve_report(kkt_result.report);
          const Point base_solution =
            subproblem.to_base_point(kkt_result.solution);
          const auto residual = product_->residual(base_solution);
          require_finite_block_values(residual.stationarity,
                                      "PDAS stationarity residual");
          require_finite_block_values(residual.equality,
                                      "PDAS equality residual");
          const Covector box_multiplier = make_box_multiplier(
            residual.stationarity,
            selection);
          const auto next_selection = complementarity_->classify(
            Primal(complementarity_->layout(),
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

    static void
    require_finite_solve_report(const QuadraticKKTSolveReport &report)
    {
      require(std::isfinite(report.stationarity_residual),
              "PDAS KKT stationarity residual is non-finite");
      require(std::isfinite(report.equality_residual),
              "PDAS KKT equality residual is non-finite");
      require(std::isfinite(report.linear_solve.relative_tolerance),
              "PDAS linear solve relative tolerance is non-finite");
      require(std::isfinite(report.linear_solve.absolute_tolerance),
              "PDAS linear solve absolute tolerance is non-finite");
      require(std::isfinite(report.linear_solve.requested_tolerance),
              "PDAS linear solve requested tolerance is non-finite");
      require(std::isfinite(report.linear_solve.achieved_residual),
              "PDAS linear solve achieved residual is non-finite");
    }

    Covector
    make_box_multiplier(const Covector &stationarity,
                        const ActiveSetSelectionT<Backend> &selection) const
    {
      require(stationarity.n_blocks() > control_block_,
              "PDAS stationarity has no control block");
      Vector values = Backend::zeros(complementarity_->layout()->dimension(0));
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
      Covector result(complementarity_->layout(), {std::move(values)});
      require_finite_block_values(result, "PDAS box multiplier");
      return result;
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
        complementarity_->multiplier_to_primal(box_multiplier);
      require_finite_block_values(
        box_multiplier, "PDAS reported box multiplier");
      require_finite_block_values(
        representative, "PDAS represented box multiplier");
      double primal_violation = 0.0;
      double dual_violation = 0.0;
      double complementarity_residual = 0.0;
      for (std::size_t index = 0;
           index < complementarity_->layout()->dimension(0);
           ++index)
        {
          const double value = Backend::value(point_control, index);
          const double multiplier =
            Backend::value(representative.block(0), index);
          const double lower = Backend::value(
            complementarity_->bounds().lower().block(0), index);
          const double upper = Backend::value(
            complementarity_->bounds().upper().block(0), index);
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
      require(std::isfinite(primal_violation),
              "PDAS primal feasibility residual is non-finite");
      require(std::isfinite(dual_violation),
              "PDAS dual feasibility residual is non-finite");
      require(std::isfinite(complementarity_residual),
              "PDAS complementarity residual is non-finite");
      require(std::isfinite(stationarity_residual),
              "PDAS stationarity residual is non-finite");
      require(std::isfinite(equality_residual),
              "PDAS equality residual is non-finite");
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

    std::shared_ptr<const Product>         product_;
    std::shared_ptr<const Complementarity> complementarity_;
    std::size_t            control_block_;
    SolveAction            solve_action_;
  };

  using ActiveSetKKTSubproblem = ActiveSetKKTSubproblemT<DenseBackend>;
  using PDASSolver = PDASSolverT<DenseBackend>;
} // namespace nmopt::contract
