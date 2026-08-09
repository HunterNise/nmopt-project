#pragma once

#include "nmopt/dealii/mass_metric.hpp"
#include "nmopt/dealii/serial_spd_solver.hpp"
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <variant>

namespace nmopt::compiler::v1
{
  struct DealiiBindingProvenance
  {
    std::string forcing;
    std::string desired_state;
    std::string fixed_dirichlet_data;
  };

  // Concrete values are supplied after semantic validation. They are not
  // stored in ProblemSpec, which remains independent of deal.II.
  template <int dim>
  struct DealiiDataBindings
  {
    DealiiDataBindings(
      const dealii::Function<dim> &forcing_function,
      const dealii::Function<dim> &desired_state_function,
      std::optional<double>        diffusion_value,
      const double                 reaction_value,
      const double                 regularisation_value,
      DealiiBindingProvenance      binding_provenance,
      std::optional<std::reference_wrapper<const dealii::Function<dim>>>
        fixed_data = std::nullopt)
      : forcing(forcing_function)
      , desired_state(desired_state_function)
      , diffusion(diffusion_value)
      , reaction(reaction_value)
      , regularisation_weight(regularisation_value)
      , fixed_dirichlet_data(fixed_data)
      , provenance(std::move(binding_provenance))
    {}

    const dealii::Function<dim> &forcing;
    const dealii::Function<dim> &desired_state;
    // Constant diffusion is required by the linear source-control targets.
    // The coefficient-identification target supplies diffusion through its
    // decision variable and therefore leaves this binding disengaged.
    std::optional<double>        diffusion;
    double                       reaction;
    double                       regularisation_weight;
    // Present only when the semantic graph declares a fixed-Dirichlet
    // reconstruction. Its boundary interpolation is the selected discrete
    // lifting rule, rather than an implicit use of the state data.
    std::optional<std::reference_wrapper<const dealii::Function<dim>>>
      fixed_dirichlet_data = std::nullopt;
    // Function objects cannot describe their own data source. These labels
    // are required compiler provenance, not executable configuration.
    DealiiBindingProvenance provenance;
  };

  using CellwiseBoundValue = std::variant<double, dealii::Vector<double>>;

  struct CellwiseBoxDataBindings
  {
    CellwiseBoundValue lower;
    CellwiseBoundValue upper;
  };

  // This is intentionally separate from CellwiseBoxDataBindings: a boundary
  // control has one coefficient per selected face rather than one per cell.
  using FacewiseBoundValue = std::variant<double, dealii::Vector<double>>;

  struct FacewiseBoxDataBindings
  {
    FacewiseBoundValue lower;
    FacewiseBoundValue upper;
  };

  struct DealiiDiscretisationPolicy
  {
    enum class Execution
    {
      assembled,
      matrix_free
    };

    unsigned int                                state_degree = 1;
    Execution                                   execution = Execution::assembled;
    dealii_backend::MassMetricSolveParameters control_metric_solve = {};
    dealii_backend::SPDLinearSolvePolicy       state_solve = {};
    dealii_backend::SPDLinearSolvePolicy       adjoint_solve = {};
  };
} // namespace nmopt::compiler::v1
