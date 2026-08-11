#pragma once

#include "nmopt/dealii/mass_metric.hpp"
#include "nmopt/dealii/serial_spd_solver.hpp"
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/vector.h>

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <variant>

namespace nmopt::compiler::v1
{
  class DealiiCompiler;

  struct DealiiBindingProvenance
  {
    std::string forcing;
    std::string desired_state;
    std::string fixed_dirichlet_data;
  };

  struct DealiiGeneralScalarBindingProvenance
  {
    std::string diffusion_tensor;
    std::string conservative_transport;
    std::string advective_transport;
    std::string reaction;
    std::string robin_coefficient;
    std::string robin_source;
  };

  struct DealiiConservativeTransportBindingProvenance
  {
    std::string conservative_transport;
  };

  // C5.6 consumes only the conservative transport coefficient in addition
  // to the registered constant diffusion-reaction data.  Keeping this
  // narrow binding separate from the full P5.1 coefficient bundle prevents
  // unused Robin and advective ports from becoming accidental requirements.
  template <int dim>
  struct DealiiConservativeTransportDataBindings
  {
    DealiiConservativeTransportDataBindings(
      const dealii::TensorFunction<1, dim> &transport_function,
      DealiiConservativeTransportBindingProvenance binding_provenance)
      : conservative_transport(transport_function)
      , provenance(std::move(binding_provenance))
    {}

    const dealii::TensorFunction<1, dim> &conservative_transport;
    DealiiConservativeTransportBindingProvenance provenance;
  };

  template <int dim>
  struct DealiiGeneralScalarDataBindings
  {
    DealiiGeneralScalarDataBindings(
      const dealii::TensorFunction<2, dim> &diffusion_tensor_function,
      const dealii::TensorFunction<1, dim> &conservative_transport_function,
      const dealii::TensorFunction<1, dim> &advective_transport_function,
      const dealii::Function<dim> &          reaction_function,
      const dealii::Function<dim> &          robin_coefficient_function,
      const dealii::Function<dim> &          robin_source_function,
      DealiiGeneralScalarBindingProvenance   binding_provenance)
      : diffusion_tensor(diffusion_tensor_function)
      , conservative_transport(conservative_transport_function)
      , advective_transport(advective_transport_function)
      , reaction(reaction_function)
      , robin_coefficient(robin_coefficient_function)
      , robin_source(robin_source_function)
      , provenance(std::move(binding_provenance))
    {}

    const dealii::TensorFunction<2, dim> &diffusion_tensor;
    const dealii::TensorFunction<1, dim> &conservative_transport;
    const dealii::TensorFunction<1, dim> &advective_transport;
    const dealii::Function<dim> &          reaction;
    const dealii::Function<dim> &          robin_coefficient;
    const dealii::Function<dim> &          robin_source;
    DealiiGeneralScalarBindingProvenance   provenance;
  };

  template <int dim>
  struct DealiiWeightedTraceDataBindings
  {
    DealiiWeightedTraceDataBindings(
      const dealii::Function<dim> &weight_function,
      std::string                  binding_provenance)
      : weight(weight_function)
      , provenance(std::move(binding_provenance))
    {}

    const dealii::Function<dim> &weight;
    std::string                  provenance;
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
        fixed_data = std::nullopt,
      std::optional<DealiiGeneralScalarDataBindings<dim>>
        general_scalar_data = std::nullopt,
      std::optional<DealiiWeightedTraceDataBindings<dim>>
        weighted_trace_data = std::nullopt,
      std::optional<DealiiConservativeTransportDataBindings<dim>>
        conservative_transport_data = std::nullopt)
      : forcing(forcing_function)
      , desired_state(desired_state_function)
      , diffusion(diffusion_value)
      , reaction(reaction_value)
      , regularisation_weight(regularisation_value)
      , fixed_dirichlet_data(fixed_data)
      , provenance(std::move(binding_provenance))
      , general_scalar(std::move(general_scalar_data))
      , weighted_trace(std::move(weighted_trace_data))
      , conservative_transport(std::move(conservative_transport_data))
    {}

    const dealii::Function<dim> &forcing;
    const dealii::Function<dim> &desired_state;
    // Constant diffusion is required by the coefficient-bound linear targets.
    // Coefficient identification supplies it through the decision variable;
    // the normalized Chapter 5.11.2 Laplacian fixes it to one. Both leave
    // this binding disengaged.
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
    // Present only for P5.1's bounded general scalar component target. The
    // rank-specific TensorFunction types make coefficient shape explicit at
    // the compiler boundary rather than interpreting scalar components.
    std::optional<DealiiGeneralScalarDataBindings<dim>> general_scalar;
    // Present only when an observation explicitly consumes a fixed boundary
    // weight. It remains separate from target data and loss configuration.
    std::optional<DealiiWeightedTraceDataBindings<dim>> weighted_trace;
    // Present only for the C5.6 Neumann-control composition. Its narrow
    // surface avoids requiring the unused P5.1 coefficient ports.
    std::optional<DealiiConservativeTransportDataBindings<dim>>
      conservative_transport;
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

  // Owns one static triangulation exclusively for the lifetime of compiled
  // products. Moving a unique_ptr into the session prevents caller mutation;
  // the compiler is the only layer allowed mutable access while lowering.
  template <int dim>
  class DealiiCompilationSession final
  {
  public:
    DealiiCompilationSession(
      std::unique_ptr<dealii::Triangulation<dim>> triangulation,
      std::string                                  mesh_provenance)
      : triangulation_(std::move(triangulation))
      , mesh_provenance_(std::move(mesh_provenance))
    {
      contract::require(static_cast<bool>(triangulation_),
                        "A deal.II compilation session needs a triangulation");
      contract::require(!mesh_provenance_.empty(),
                        "A deal.II compilation session needs mesh provenance");
    }

    const dealii::Triangulation<dim> &
    triangulation() const
    {
      return *triangulation_;
    }

    const std::string &
    mesh_provenance() const
    {
      return mesh_provenance_;
    }

  private:
    dealii::Triangulation<dim> &
    mutable_triangulation() const
    {
      return *triangulation_;
    }

    std::unique_ptr<dealii::Triangulation<dim>> triangulation_;
    std::string                                  mesh_provenance_;

    friend class DealiiCompiler;
  };
} // namespace nmopt::compiler::v1
