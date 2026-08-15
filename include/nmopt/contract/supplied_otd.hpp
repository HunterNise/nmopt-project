#pragma once

#include "nmopt/contract/linear_solve.hpp"

#include <cstddef>
#include <functional>
#include <string>
#include <utility>

namespace nmopt::contract
{
  // The selected P6.2 product has one state, one supplied adjoint, and one
  // decision block.  The residual layout contains the corresponding state,
  // adjoint, and stationarity equation blocks.  Space identifiers remain
  // caller-owned; this descriptor only records their block roles.
  struct SuppliedOTDBlockSelection
  {
    std::size_t state_variable = 0;
    std::size_t adjoint_variable = 1;
    std::size_t control_variable = 2;
    std::size_t state_equation = 0;
    std::size_t adjoint_equation = 1;
    std::size_t control_stationarity = 2;
  };

  class SuppliedOTDLayout
  {
  public:
    SuppliedOTDLayout(LayoutPtr                    variable_layout,
                      LayoutPtr                    residual_layout,
                      SuppliedOTDBlockSelection selection = {})
      : variable_layout_(std::move(variable_layout))
      , residual_layout_(std::move(residual_layout))
      , selection_(selection)
    {
      require(static_cast<bool>(variable_layout_),
              "Supplied OTD layout needs a variable layout");
      require(static_cast<bool>(residual_layout_),
              "Supplied OTD layout needs a residual layout");
      require(variable_layout_->n_blocks() == 3,
              "Supplied OTD variable layout needs state, adjoint, and control blocks");
      require(residual_layout_->n_blocks() == 3,
              "Supplied OTD residual layout needs state, adjoint, and stationarity blocks");
      validate_selection(selection_.state_variable,
                         selection_.adjoint_variable,
                         selection_.control_variable,
                         variable_layout_->n_blocks(),
                         "variable");
      validate_selection(selection_.state_equation,
                         selection_.adjoint_equation,
                         selection_.control_stationarity,
                         residual_layout_->n_blocks(),
                         "residual");
    }

    const LayoutPtr &
    variable_layout() const
    {
      return variable_layout_;
    }

    const LayoutPtr &
    residual_layout() const
    {
      return residual_layout_;
    }

    const SuppliedOTDBlockSelection &
    selection() const
    {
      return selection_;
    }

  private:
    static void
    validate_selection(const std::size_t first,
                       const std::size_t second,
                       const std::size_t third,
                       const std::size_t block_count,
                       const char *     layout_name)
    {
      require(first < block_count && second < block_count &&
                third < block_count,
              std::string("Supplied OTD ") + layout_name +
                " block selection is outside its layout");
      require(first != second && first != third && second != third,
              std::string("Supplied OTD ") + layout_name +
                " block selection must be distinct");
    }

    LayoutPtr                  variable_layout_;
    LayoutPtr                  residual_layout_;
    SuppliedOTDBlockSelection  selection_;
  };

  template <typename Backend>
  class SuppliedOTDSystemT final
  {
  public:
    using Primal = PrimalBlockT<Backend>;
    using Covector = CovectorBlockT<Backend>;
    using SolveResult = FormulationSolveResultT<Backend>;

    using ResidualAction = std::function<Covector(const Primal &)>;
    using LinearizedAction =
      std::function<Covector(const Primal &, const Primal &)>;
    using TransposeAction =
      std::function<Covector(const Primal &, const Primal &)>;
    using SolveAction = std::function<SolveResult(const Primal &)>;

    SuppliedOTDSystemT(SuppliedOTDLayout layout,
                       ResidualAction residual,
                       LinearizedAction residual_jvp,
                       TransposeAction residual_vjp,
                       SolveAction solve)
      : layout_(std::move(layout))
      , residual_(std::move(residual))
      , residual_jvp_(std::move(residual_jvp))
      , residual_vjp_(std::move(residual_vjp))
      , solve_(std::move(solve))
    {
      require(static_cast<bool>(residual_),
              "Supplied OTD system needs a residual action");
      require(static_cast<bool>(residual_jvp_),
              "Supplied OTD system needs a residual JVP action");
      require(static_cast<bool>(residual_vjp_),
              "Supplied OTD system needs a residual VJP action");
      require(static_cast<bool>(solve_),
              "Supplied OTD system needs a solve action");
    }

    const LayoutPtr &
    variable_layout() const
    {
      return layout_.variable_layout();
    }

    const LayoutPtr &
    residual_layout() const
    {
      return layout_.residual_layout();
    }

    const SuppliedOTDBlockSelection &
    block_selection() const
    {
      return layout_.selection();
    }

    Covector
    residual(const Primal &point) const
    {
      require_variable(point, "Supplied OTD residual");
      return checked_residual(residual_(point), "Supplied OTD residual");
    }

    Covector
    residual_jvp(const Primal &point, const Primal &tangent) const
    {
      require_variable(point, "Supplied OTD residual JVP point");
      require_variable(tangent, "Supplied OTD residual JVP tangent");
      return checked_residual(residual_jvp_(point, tangent),
                              "Supplied OTD residual JVP");
    }

    Covector
    residual_vjp(const Primal &point, const Primal &residual_seed) const
    {
      require_variable(point, "Supplied OTD residual VJP point");
      require(residual_seed.layout()->compatible_with(*residual_layout()),
              "Supplied OTD residual VJP seed has an incompatible layout");
      return checked_variable_covector(residual_vjp_(point, residual_seed),
                                       "Supplied OTD residual VJP");
    }

    SolveResult
    solve(const Primal &initial_point) const
    {
      require_variable(initial_point, "Supplied OTD solve initial point");
      SolveResult result = solve_(initial_point);
      require_variable(result.solution, "Supplied OTD solve result");
      return result;
    }

    Covector
    state_residual(const Primal &point) const
    {
      return extract_covector_block(residual(point),
                                    block_selection().state_equation,
                                    "state_equation");
    }

    Covector
    adjoint_residual(const Primal &point) const
    {
      return extract_covector_block(residual(point),
                                    block_selection().adjoint_equation,
                                    "adjoint_equation");
    }

    Covector
    control_stationarity(const Primal &point) const
    {
      return extract_covector_block(residual(point),
                                    block_selection().control_stationarity,
                                    "control_stationarity");
    }

  private:
    void
    require_variable(const Primal &value, const char *operation) const
    {
      require(value.layout()->compatible_with(*variable_layout()),
              std::string(operation) + " has an incompatible variable layout");
    }

    Covector
    checked_residual(Covector value, const char *operation) const
    {
      require(value.layout()->compatible_with(*residual_layout()),
              std::string(operation) + " returned an incompatible layout");
      return value;
    }

    Covector
    checked_variable_covector(Covector value, const char *operation) const
    {
      require(value.layout()->compatible_with(*variable_layout()),
              std::string(operation) + " returned an incompatible layout");
      return value;
    }

    SuppliedOTDLayout layout_;
    ResidualAction    residual_;
    LinearizedAction  residual_jvp_;
    TransposeAction   residual_vjp_;
    SolveAction       solve_;
  };

  using SuppliedOTDSystem = SuppliedOTDSystemT<DenseBackend>;
} // namespace nmopt::contract
