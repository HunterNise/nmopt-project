#pragma once

#include "nmopt/dealii/mass_metric.hpp"

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

#include <functional>
#include <optional>
#include <variant>

namespace nmopt::compiler::v1
{
  // Concrete values are supplied after semantic validation. They are not
  // stored in ProblemSpec, which remains independent of deal.II.
  template <int dim>
  struct DealiiDataBindings
  {
    const dealii::Function<dim> &forcing;
    const dealii::Function<dim> &desired_state;
    double                       diffusion;
    double                       reaction;
    double                       regularisation_weight;
    // Present only when the semantic graph declares a fixed-Dirichlet
    // reconstruction. Its boundary interpolation is the selected discrete
    // lifting rule, rather than an implicit use of the state data.
    std::optional<std::reference_wrapper<const dealii::Function<dim>>>
      fixed_dirichlet_data = std::nullopt;
  };

  using CellwiseBoundValue = std::variant<double, dealii::Vector<double>>;

  struct CellwiseBoxDataBindings
  {
    CellwiseBoundValue lower;
    CellwiseBoundValue upper;
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
  };
} // namespace nmopt::compiler::v1
