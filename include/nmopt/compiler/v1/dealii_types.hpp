#pragma once

#include "nmopt/dealii/mass_metric.hpp"

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

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
