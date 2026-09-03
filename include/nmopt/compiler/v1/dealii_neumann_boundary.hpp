#pragma once

#include "nmopt/compiler/v1/dealii_neumann_control_realisation.hpp"
#include "nmopt/compiler/v1/dealii_reference_cell.hpp"
#include "nmopt/compiler/v1/dealii_volume_observation.hpp"
#include "nmopt/contract/executable_model.hpp"
#include "nmopt/dealii/serial_backend.hpp"
#include "nmopt/dealii/serial_spd_solver.hpp"

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <locale>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nmopt::compiler::v1::detail
{
  // V1-only Neumann target. It owns the selected boundary-control realization
  // and lowers
  // exactly the declared natural-boundary contribution
  //
  //   r(y,u) = A y - f_h - C_Gamma u,
  //
  // where each column of C_Gamma is the FEFaceValues realization of one
  // facewise coefficient or continuous trace basis function. It intentionally
  // does not reuse the volume-control model or present this term as a
  // Dirichlet/control lifting.
  template <int dim>
  class NeumannBoundaryControlModel final
    : public contract::ExecutableModelT<dealii_backend::SerialBackend>
  {
  public:
    using Backend = dealii_backend::SerialBackend;
    using Vector = dealii::Vector<double>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;
    using SolveResult = contract::FormulationSolveResultT<Backend>;

    enum class StateGauge
    {
      fixed_dirichlet,
      mean_zero_multiplier
    };

    enum class StateObservation
    {
      boundary_trace,
      volume_restriction
    };

    enum class WeightedTraceRealisation
    {
      fe_q_face_quadrature
    };

    NeumannBoundaryControlModel(
      dealii::Triangulation<dim> &                triangulation,
      const dealii::Function<dim> &               forcing,
      const dealii::Function<dim> &               desired_state,
      const double                                  diffusion,
      const double                                  reaction,
      const double                                  regularisation_weight,
      const unsigned int                            state_degree,
      std::set<dealii::types::boundary_id>         dirichlet_boundary_ids,
      std::set<dealii::types::boundary_id>         control_boundary_ids,
      std::set<dealii::types::boundary_id>         observation_boundary_ids,
      const unsigned int                            volume_observation_quadrature_order,
      const VolumeObservationTargetRealisation      volume_observation_target_realisation,
      const dealii::Function<dim> *                 observation_weight = nullptr,
      const StateGauge                              state_gauge =
        StateGauge::fixed_dirichlet,
      const StateObservation                       state_observation =
        StateObservation::boundary_trace,
      std::set<dealii::types::material_id>         observation_material_ids = {},
      const dealii::TensorFunction<1, dim> *       conservative_transport = nullptr,
      std::set<dealii::types::boundary_id>         transport_boundary_ids = {},
      const TransportBoundaryRealisation           transport_boundary_realisation =
        TransportBoundaryRealisation::total_conormal,
      const std::optional<WeightedTraceRealisation> weighted_trace_realisation =
        std::nullopt,
      std::optional<std::reference_wrapper<const dealii::Function<dim>>>
        fixed_dirichlet_data = std::nullopt,
      const NeumannControlRealisationKind control_realisation_kind =
        NeumannControlRealisationKind::facewise_constant,
      const dealii::Function<dim> *            natural_boundary_source = nullptr,
      std::set<dealii::types::boundary_id>      natural_boundary_source_ids = {})
      : state_fe_(make_scalar_lagrange_element(
          triangulation, state_degree, "The Neumann v1 target"))
      , state_dof_handler_(triangulation)
      , diffusion_(diffusion)
      , reaction_(reaction)
      , regularisation_weight_(regularisation_weight)
      , state_gauge_(state_gauge)
      , dirichlet_boundary_ids_(std::move(dirichlet_boundary_ids))
      , observation_boundary_ids_(std::move(observation_boundary_ids))
      , state_observation_(state_observation)
      , observation_material_ids_(std::move(observation_material_ids))
      , uses_conservative_transport_(conservative_transport != nullptr)
      , transport_boundary_ids_(std::move(transport_boundary_ids))
      , transport_boundary_realisation_(transport_boundary_realisation)
      , weighted_trace_realisation_(weighted_trace_realisation)
      , has_fixed_dirichlet_data_(fixed_dirichlet_data.has_value())
      , natural_boundary_source_(natural_boundary_source)
      , natural_boundary_source_ids_(std::move(natural_boundary_source_ids))
    {
      contract::require(diffusion_ > 0.0,
                        "Diffusion coefficient must be strictly positive");
      contract::require(reaction_ >= 0.0,
                        "Reaction coefficient must be non-negative");
      contract::require(regularisation_weight_ > 0.0,
                        "Control regularisation weight must be strictly positive");
      contract::require(state_degree > 0,
                        "State FE degree must be at least one");
      contract::require(
        state_gauge_ != StateGauge::fixed_dirichlet ||
          !dirichlet_boundary_ids_.empty(),
        "The Neumann v1 target needs a fixed Dirichlet boundary for this gauge");
      contract::require(
        state_gauge_ != StateGauge::mean_zero_multiplier ||
          dirichlet_boundary_ids_.empty(),
        "The pure-Neumann mean-zero gauge cannot also fix boundary DoFs");
      contract::require(
        state_gauge_ != StateGauge::mean_zero_multiplier ||
          std::abs(reaction_) <= 1e-14,
        "The pure-Neumann mean-zero gauge requires zero reaction");
      contract::require(!control_boundary_ids.empty(),
                        "The Neumann v1 target needs a marked control boundary");
      contract::require(
        observation_weight == nullptr ||
          (weighted_trace_realisation_ &&
           *weighted_trace_realisation_ ==
             WeightedTraceRealisation::fe_q_face_quadrature),
        "The weighted Neumann trace needs its selected face-quadrature realization");
      contract::require(
        state_observation_ == StateObservation::boundary_trace
          ? !observation_boundary_ids_.empty()
          : !observation_material_ids_.empty(),
        "The Neumann v1 target needs the declared observation region");
      contract::require(
        state_observation_ == StateObservation::boundary_trace ||
          observation_weight == nullptr,
        "The C5.6 volume observation does not consume boundary-weight data");
      contract::require(
        transport_boundary_realisation_ ==
            TransportBoundaryRealisation::total_conormal ||
          uses_conservative_transport_,
        "The ordinary transport boundary realization needs conservative transport data");
      contract::require(
        natural_boundary_source_ == nullptr ||
          !natural_boundary_source_ids_.empty(),
        "The natural-boundary source needs a marked source boundary");
      contract::require(
        natural_boundary_source_ != nullptr ||
          natural_boundary_source_ids_.empty(),
        "Natural-boundary source ids need a source Function");
      if (natural_boundary_source_ != nullptr)
        contract::require(natural_boundary_source_->n_components == 1,
                          "The natural-boundary source must be scalar");

      state_dof_handler_.distribute_dofs(*state_fe_);
      build_constraints(fixed_dirichlet_data);
      const double coupling_scale =
        uses_ordinary_transport_boundary_realisation() ? diffusion_ : 1.0;
      switch (control_realisation_kind)
        {
          case NeumannControlRealisationKind::facewise_constant:
            control_realisation_ =
              std::make_unique<FacewiseNeumannControlRealisation<dim>>(
                state_dof_handler_,
                *state_fe_,
                constrained_state_dofs_,
                std::move(control_boundary_ids),
                state_fe_->degree + 2,
                coupling_scale);
            break;
          case NeumannControlRealisationKind::continuous_p1_trace:
            control_realisation_ =
              std::make_unique<ContinuousNeumannControlRealisation<dim>>(
                state_dof_handler_,
                *state_fe_,
                constrained_state_dofs_,
                std::move(control_boundary_ids),
                state_fe_->degree + 2,
                coupling_scale);
            break;
        }
      initialise_storage();
      if (state_observation_ == StateObservation::volume_restriction)
        volume_observation_ =
          std::make_unique<VolumeObservationAssembly<dim>>(
            state_dof_handler_,
            *state_fe_,
            constrained_state_dofs_,
            fixed_state_values_,
            observation_material_ids_,
            desired_state,
            volume_observation_quadrature_order,
            volume_observation_target_realisation);
      assemble(forcing,
               desired_state,
               conservative_transport,
               observation_weight);
      if (uses_mean_zero_gauge())
        build_mean_zero_system();
      if (uses_conservative_transport_)
        {
          nonsymmetric_solver_ = std::make_unique<dealii::SparseDirectUMFPACK>();
          nonsymmetric_solver_->initialize(system_matrix_);
        }
    }

    const contract::LayoutPtr &
    variable_layout() const override
    {
      return variable_layout_;
    }

    const contract::LayoutPtr &
    test_layout() const override
    {
      return test_layout_;
    }

    std::size_t
    physical_state_dimension() const
    {
      return state_dof_handler_.n_dofs();
    }

    std::size_t
    independent_state_dimension() const
    {
      return static_cast<std::size_t>(std::count(constrained_state_dofs_.begin(),
                                                 constrained_state_dofs_.end(),
                                                 false));
    }

    std::size_t
    physical_control_dimension() const
    {
      return control_realisation_->dimension();
    }

    std::size_t
    independent_control_dimension() const
    {
      return control_realisation_->dimension();
    }

    const std::vector<dealii::Point<dim>> &
    control_coordinates() const
    {
      return control_realisation_->coordinates();
    }

    std::size_t
    realized_observation_dimension() const
    {
      return state_observation_ == StateObservation::boundary_trace
               ? observation_evaluations_.size()
               : physical_state_dimension();
    }

    // The boundary-trace map owns the same ordered face-quadrature samples
    // used to assemble the tracking objective.  For the weighted target the
    // stored evaluations include h, so these ports describe y -> h gamma y.
    std::vector<double>
    boundary_trace_values(const Primal &variables) const
    {
      require_variables(variables, "Boundary-trace value");
      require_boundary_trace_observation("Boundary-trace values");
      std::vector<double> values;
      values.reserve(observation_evaluations_.size());
      for (const auto &evaluation : observation_evaluations_)
        values.push_back(evaluation * variables.block(0));
      return values;
    }

    Covector
    boundary_trace_jvp(const Primal &variable_tangent) const
    {
      require_variables(variable_tangent, "Boundary-trace JVP tangent");
      require_boundary_trace_observation("Boundary-trace JVP");
      Vector values(observation_evaluations_.size());
      for (std::size_t index = 0; index < observation_evaluations_.size(); ++index)
        values[index] = observation_evaluations_[index] *
                        variable_tangent.block(0);
      const std::size_t observation_dimension = values.size();
      return Covector(std::make_shared<const contract::BlockLayout>(
                        "boundary_trace_observation",
                        std::vector<contract::SpaceId>{{"boundary_trace"}},
                        std::vector<std::size_t>{observation_dimension}),
                      {std::move(values)});
    }

    Covector
    boundary_trace_vjp(const std::vector<double> &seed) const
    {
      require_boundary_trace_observation("Boundary-trace VJP");
      contract::require(seed.size() == observation_evaluations_.size(),
                        "Boundary-trace VJP seed has the wrong dimension");
      Vector physical_covector(state_dof_handler_.n_dofs());
      for (std::size_t index = 0; index < seed.size(); ++index)
        physical_covector.add(seed[index] * observation_quadrature_weights_[index],
                              observation_evaluations_[index]);
      return Covector(state_layout_, {std::move(physical_covector)});
    }

    const std::vector<double> &
    boundary_trace_quadrature_weights() const
    {
      require_boundary_trace_observation("Boundary-trace quadrature weights");
      return observation_quadrature_weights_;
    }

    dealii_backend::MassMetric
    control_l2_metric(
      dealii_backend::MassMetricSolveParameters solve_parameters = {}) const
    {
      return control_realisation_->l2_metric(solve_parameters);
    }

    dealii_backend::FacewiseBoxConstraint
    control_l2_box_constraint(
      Vector                              lower,
      Vector                              upper,
      const dealii_backend::MassMetric & projection_metric) const
    {
      const auto *facewise = dynamic_cast<const
        FacewiseNeumannControlRealisation<dim> *>(control_realisation_.get());
      contract::require(facewise != nullptr,
                        "Continuous Neumann control has no coefficientwise box realization");
      return facewise->l2_box_constraint(std::move(lower),
                                         std::move(upper),
                                         projection_metric);
    }

    dealii_backend::FacewiseBoxConstraint
    control_l2_box_constraint(
      const double                        lower,
      const double                        upper,
      const dealii_backend::MassMetric & projection_metric) const
    {
      const auto *facewise = dynamic_cast<const
        FacewiseNeumannControlRealisation<dim> *>(control_realisation_.get());
      contract::require(facewise != nullptr,
                        "Continuous Neumann control has no coefficientwise box realization");
      return facewise->l2_box_constraint(lower, upper, projection_metric);
    }

    bool
    uses_mean_zero_gauge() const
    {
      return state_gauge_ == StateGauge::mean_zero_multiplier;
    }

    bool
    forcing_is_compatible() const
    {
      return is_compatible(forcing_load_);
    }

    double
    state_mean(const Primal &state) const
    {
      contract::require(
        state.layout()->compatible_with(*state_layout_) ||
          state.layout()->compatible_with(*test_layout_),
        "State mean received an incompatible state or adjoint layout");
      return mean_constraint_ * state.block(0);
    }

    void
    write_native_output(const std::filesystem::path &directory,
                        const Primal &                 state,
                        const Primal &                 control,
                        const Primal &                 adjoint,
                        const Primal *                 uncontrolled_state = nullptr,
                        const dealii::Function<dim> *  forcing = nullptr,
                        const dealii::Function<dim> *  desired_state = nullptr) const
    {
      static_assert(dim == 2,
                    "Neumann field output currently supports two dimensions");
      contract::require(state.layout()->compatible_with(*state_layout_),
                        "Native output state has an incompatible layout");
      contract::require(control.layout()->compatible_with(*control_layout_),
                        "Native output control has an incompatible layout");
      contract::require(adjoint.layout()->compatible_with(*test_layout_),
                        "Native output adjoint has an incompatible layout");
      if (uncontrolled_state != nullptr)
        contract::require(
          uncontrolled_state->layout()->compatible_with(*state_layout_),
          "Native output uncontrolled state has an incompatible layout");
      contract::require(control.block(0).size() ==
                          control_realisation_->dimension(),
                        "Native output control has an incompatible size");
      contract::require((forcing == nullptr) == (desired_state == nullptr),
                        "Native output needs both forcing and target functions");
      std::filesystem::create_directories(directory);

      dealii::DataOut<dim> mesh_out;
      mesh_out.attach_triangulation(state_dof_handler_.get_triangulation());
      mesh_out.build_patches();

      std::ofstream mesh_output(directory / "mesh-volume.vtu");
      if (!mesh_output)
        throw std::runtime_error("could not open Neumann volume mesh output");
      mesh_output.imbue(std::locale::classic());
      mesh_out.write_vtu(mesh_output);
      if (!mesh_output)
        throw std::runtime_error("could not write Neumann volume mesh output");

      dealii::GridOut grid_out;
      std::ofstream  svg_output(directory / "mesh-volume.svg");
      if (!svg_output)
        throw std::runtime_error(
          "could not open Neumann volume mesh SVG output");
      svg_output.imbue(std::locale::classic());
      grid_out.write_svg(state_dof_handler_.get_triangulation(), svg_output);
      if (!svg_output)
        throw std::runtime_error(
          "could not write Neumann volume mesh SVG output");

      dealii::DataOut<dim> data_out;
      data_out.attach_dof_handler(state_dof_handler_);
      data_out.add_data_vector(state.block(0), "state");
      if (uncontrolled_state != nullptr)
        data_out.add_data_vector(uncontrolled_state->block(0),
                                 "state_uncontrolled");
      data_out.add_data_vector(adjoint.block(0), "adjoint");
      Vector negative_adjoint = adjoint.block(0);
      negative_adjoint *= -1.0;
      data_out.add_data_vector(negative_adjoint, "negative_adjoint");
      if (forcing != nullptr)
        {
          Vector forcing_values(state_dof_handler_.n_dofs());
          Vector desired_state_values(state_dof_handler_.n_dofs());
          dealii::VectorTools::interpolate(state_dof_handler_,
                                            *forcing,
                                            forcing_values);
          dealii::VectorTools::interpolate(state_dof_handler_,
                                            *desired_state,
                                            desired_state_values);
          data_out.add_data_vector(forcing_values, "forcing");
          data_out.add_data_vector(desired_state_values, "target");
        }
      if (state_observation_ == StateObservation::volume_restriction)
        {
          dealii::Vector<double> observation_region(
            state_dof_handler_.get_triangulation().n_active_cells());
          observation_region = 0.0;
          for (auto cell = state_dof_handler_.begin_active();
               cell != state_dof_handler_.end();
               ++cell)
            observation_region[cell->active_cell_index()] =
              observation_material_ids_.count(cell->material_id()) != 0 ? 1.0 :
                                                                            0.0;
          data_out.add_data_vector(observation_region,
                                   "observation_region",
                                   dealii::DataOut<dim>::type_cell_data);
        }
      data_out.build_patches();

      std::ofstream fields_output(directory / "fields-volume.vtu");
      if (!fields_output)
        throw std::runtime_error("could not open Neumann volume field output");
      fields_output.imbue(std::locale::classic());
      data_out.write_vtu(fields_output);
      if (!fields_output)
        throw std::runtime_error("could not write Neumann volume field output");

      control_realisation_->write_native_output(
        directory / "control-boundary.vtu",
        control.block(0));
    }

    Covector
    residual(const Primal &variables) const override
    {
      require_variables(variables, "Residual");
      Vector value(state_dof_handler_.n_dofs());
      system_matrix_.vmult(value, variables.block(0));
      value.add(1.0, fixed_state_load_);
      value.add(-1.0, forcing_load_);

      Vector control_contribution =
        control_realisation_->coupling_action(variables.block(1));
      value.add(-1.0, control_contribution);
      zero_constrained_entries(value);
      return Covector(test_layout_, {std::move(value)});
    }

    Covector
    residual_jvp(const Primal &variables,
                 const Primal &variable_tangent) const override
    {
      require_variables(variables, "Residual JVP");
      require_variables(variable_tangent, "Residual JVP tangent");
      Vector value(state_dof_handler_.n_dofs());
      system_matrix_.vmult(value, variable_tangent.block(0));
      Vector control_contribution =
        control_realisation_->coupling_action(variable_tangent.block(1));
      value.add(-1.0, control_contribution);
      zero_constrained_entries(value);
      return Covector(test_layout_, {std::move(value)});
    }

    Covector
    residual_vjp(const Primal &variables,
                 const Primal &test_seed) const override
    {
      require_variables(variables, "Residual VJP");
      contract::require(test_seed.layout()->compatible_with(*test_layout_),
                        "Residual VJP seed has an incompatible test layout");

      Vector state(state_dof_handler_.n_dofs());
      system_matrix_.Tvmult(state, test_seed.block(0));
      zero_constrained_entries(state);
      Vector control =
        control_realisation_->coupling_transpose_action(test_seed.block(0));
      control *= -1.0;
      return Covector(variable_layout_, {std::move(state), std::move(control)});
    }

    double
    objective(const Primal &variables) const override
    {
      const auto components = objective_components(variables);
      return components.state_tracking + components.control_regularisation;
    }

    struct ObjectiveComponents
    {
      double state_tracking;
      double control_regularisation;
    };

    ObjectiveComponents
    objective_components(const Primal &variables) const
    {
      require_variables(variables, "Objective");
      Vector tracked_state(state_dof_handler_.n_dofs());
      state_tracking_matrix().vmult(tracked_state, variables.block(0));
      const double state_value =
        0.5 * (variables.block(0) * tracked_state) -
        (desired_state_load() * variables.block(0)) +
        0.5 * desired_state_norm();

      const double control_value =
        control_realisation_->regularisation_objective(
          variables.block(1),
          regularisation_weight_);
      return {state_value, control_value};
    }

    Covector
    objective_derivative(const Primal &variables) const override
    {
      require_variables(variables, "Objective derivative");
      Vector state(state_dof_handler_.n_dofs());
      state_tracking_matrix().vmult(state, variables.block(0));
      state.add(-1.0, desired_state_load());

      Vector control = control_realisation_->regularisation_derivative(
        variables.block(1),
        regularisation_weight_);
      return Covector(variable_layout_, {std::move(state), std::move(control)});
    }

    Primal
    solve_state(const Primal &control) const
    {
      auto result = solve_state_with_report(control, {});
      contract::require(result.report.converged(),
                        "State solve did not converge under its declared policy");
      return std::move(result.solution);
    }

    SolveResult
    solve_state_with_report(
      const Primal &                              control,
      const dealii_backend::SPDLinearSolvePolicy &policy) const
    {
      contract::require(control.layout()->compatible_with(*control_layout_),
                        "State solve control has an incompatible layout");
      Vector right_hand_side = forcing_load_;
      right_hand_side.add(-1.0, fixed_state_load_);
      Vector control_contribution =
        control_realisation_->coupling_action(control.block(0));
      right_hand_side.add(1.0, control_contribution);

      Vector state(state_dof_handler_.n_dofs());
      contract::LinearSolveReport report;
      if (uses_mean_zero_gauge())
        {
          contract::require(
            is_compatible(right_hand_side),
            "Pure-Neumann state load violates the discrete constant-mode compatibility condition");
          solve_mean_zero_system(state, right_hand_side);
          report = dealii_backend::direct_solve_report(
            "serial_sparse_direct_umfpack_mean_zero_saddle");
        }
      else if (uses_conservative_transport_)
        {
          (void)policy;
          nonsymmetric_solver_->vmult(state, right_hand_side);
          report = dealii_backend::direct_solve_report("serial_sparse_direct_umfpack");
        }
      else
        report = solve_symmetric_system(state, right_hand_side, policy);
      state_constraints_.distribute(state);
      return {Primal(state_layout_, {std::move(state)}), std::move(report)};
    }

    Primal
    solve_adjoint(const Primal &full_point,
                  const Covector &state_objective_derivative) const
    {
      auto result = solve_adjoint_with_report(full_point,
                                              state_objective_derivative,
                                              {});
      contract::require(result.report.converged(),
                        "Adjoint solve did not converge under its declared policy");
      return std::move(result.solution);
    }

    SolveResult
    solve_adjoint_with_report(
      const Primal &                              full_point,
      const Covector &                            state_objective_derivative,
      const dealii_backend::SPDLinearSolvePolicy &policy) const
    {
      require_variables(full_point, "Adjoint solve point");
      contract::require(
        state_objective_derivative.layout()->compatible_with(*state_layout_),
        "Adjoint solve right-hand side has an incompatible state layout");

      Vector adjoint(test_layout_->dimension(0));
      contract::LinearSolveReport report;
      if (uses_mean_zero_gauge())
        {
          solve_mean_zero_system(adjoint, state_objective_derivative.block(0));
          report = dealii_backend::direct_solve_report(
            "serial_sparse_direct_umfpack_mean_zero_saddle");
        }
      else if (uses_conservative_transport_)
        {
          (void)policy;
          nonsymmetric_solver_->Tvmult(adjoint, state_objective_derivative.block(0));
          report = dealii_backend::direct_solve_report(
            "serial_sparse_direct_umfpack_transpose");
        }
      else
        report = solve_symmetric_system(adjoint,
                                        state_objective_derivative.block(0),
                                        policy);
      if (has_fixed_dirichlet_data_)
        zero_constrained_entries(adjoint);
      else
        state_constraints_.distribute(adjoint);
      return {Primal(test_layout_, {std::move(adjoint)}), std::move(report)};
    }

  private:
    void
    zero_constrained_entries(Vector &values) const
    {
      if (!has_fixed_dirichlet_data_)
        return;
      for (dealii::types::global_dof_index index = 0;
           index < values.size();
           ++index)
        if (constrained_state_dofs_.at(index))
          values[index] = 0.0;
    }

    void
    require_boundary_trace_observation(const char *operation) const
    {
      contract::require(
        state_observation_ == StateObservation::boundary_trace,
        std::string(operation) + " needs the boundary-trace observation target");
    }

    void
    require_variables(const Primal &variables, const char *operation) const
    {
      contract::require(
        variables.layout()->compatible_with(*variable_layout_),
        std::string(operation) + " received an incompatible variable layout");
    }

    bool
    is_control_face(const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
                    const unsigned int face) const
    {
      return control_realisation_->is_control_face(cell, face);
    }

    bool
    is_observation_face(
      const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      const unsigned int face) const
    {
      return state_observation_ == StateObservation::boundary_trace &&
             cell->face(face)->at_boundary() &&
             observation_boundary_ids_.count(cell->face(face)->boundary_id()) != 0;
    }

    bool
    uses_ordinary_transport_boundary_realisation() const
    {
      return transport_boundary_realisation_ ==
             TransportBoundaryRealisation::ordinary_normal_transport;
    }

    bool
    is_transport_boundary_face(
      const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      const unsigned int face) const
    {
      return cell->face(face)->at_boundary() &&
             (is_control_face(cell, face) ||
              transport_boundary_ids_.count(cell->face(face)->boundary_id()) != 0);
    }

    const dealii::SparseMatrix<double> &
    state_tracking_matrix() const
    {
      return volume_observation_ == nullptr ?
               state_tracking_matrix_ :
               volume_observation_->state_tracking_matrix();
    }

    const Vector &
    desired_state_load() const
    {
      return volume_observation_ == nullptr ?
               desired_state_load_ :
               volume_observation_->desired_state_load();
    }

    double
    desired_state_norm() const
    {
      return volume_observation_ == nullptr ?
               desired_state_norm_ :
               volume_observation_->desired_state_norm();
    }

    void
    build_constraints(
      const std::optional<std::reference_wrapper<const dealii::Function<dim>>>
        fixed_dirichlet_data)
    {
      state_constraints_.clear();
      dealii::DoFTools::make_hanging_node_constraints(state_dof_handler_,
                                                       state_constraints_);
      dealii::Functions::ZeroFunction<dim> zero;
      const dealii::Function<dim> &physical_dirichlet_data =
        fixed_dirichlet_data ? fixed_dirichlet_data->get() : zero;
      for (const auto boundary_id : dirichlet_boundary_ids_)
        dealii::VectorTools::interpolate_boundary_values(state_dof_handler_,
                                                          boundary_id,
                                                          physical_dirichlet_data,
                                                          state_constraints_);
      state_constraints_.close();

      std::map<dealii::types::global_dof_index, double> dirichlet_values;
      for (const auto boundary_id : dirichlet_boundary_ids_)
        dealii::VectorTools::interpolate_boundary_values(state_dof_handler_,
                                                          boundary_id,
                                                          physical_dirichlet_data,
                                                          dirichlet_values);
      fixed_state_values_.reinit(state_dof_handler_.n_dofs());
      constrained_state_dofs_.assign(state_dof_handler_.n_dofs(), false);
      for (const auto &entry : dirichlet_values)
        {
          constrained_state_dofs_.at(entry.first) = true;
          fixed_state_values_[entry.first] = entry.second;
        }
      for (dealii::types::global_dof_index index = 0;
           index < state_dof_handler_.n_dofs();
           ++index)
        if (state_constraints_.is_constrained(index) &&
            !constrained_state_dofs_.at(index))
          throw contract::ContractError(
            "The Neumann v1 target does not support hanging, periodic, or "
            "other non-Dirichlet affine constraints");
    }

    void
    initialise_storage()
    {
      const auto state_size = state_dof_handler_.n_dofs();
      variable_layout_ = std::make_shared<const contract::BlockLayout>(
        "neumann_boundary_variables",
        std::vector<contract::SpaceId>{{"state"}, {"control"}},
        std::vector<std::size_t>{state_size,
                                 control_realisation_->dimension()});
      test_layout_ = std::make_shared<const contract::BlockLayout>(
        "neumann_boundary_state_test",
        std::vector<contract::SpaceId>{{"state_test"}},
        std::vector<std::size_t>{state_size});
      state_layout_ = variable_layout_->single_block(0, "state");
      control_layout_ = control_realisation_->layout();

      dealii::DynamicSparsityPattern state_dsp(state_size, state_size);
      dealii::DoFTools::make_sparsity_pattern(state_dof_handler_, state_dsp);
      state_sparsity_.copy_from(state_dsp);
      system_matrix_.reinit(state_sparsity_);
      state_tracking_matrix_.reinit(state_sparsity_);

      forcing_load_.reinit(state_size);
      fixed_state_load_.reinit(state_size);
      desired_state_load_.reinit(state_size);
      mean_constraint_.reinit(state_size);
      constant_mode_.reinit(state_size);
      for (dealii::types::global_dof_index index = 0; index < state_size; ++index)
        constant_mode_[index] = 1.0;
    }

    void
    assemble(const dealii::Function<dim> &forcing,
             const dealii::Function<dim> &desired_state,
             const dealii::TensorFunction<1, dim> *conservative_transport,
             const dealii::Function<dim> *observation_weight)
    {
      const unsigned int quadrature_order = state_fe_->degree + 2;
      const auto volume_quadrature = make_gauss_volume_quadrature(
        state_dof_handler_.get_triangulation(),
        quadrature_order,
        "The Neumann v1 target");
      dealii::FEValues<dim> state_values(
        *state_fe_,
        *volume_quadrature,
        dealii::update_values | dealii::update_gradients |
          dealii::update_quadrature_points | dealii::update_JxW_values);
      dealii::FullMatrix<double> local_system(state_fe_->dofs_per_cell,
                                              state_fe_->dofs_per_cell);
      dealii::Vector<double> local_forcing(state_fe_->dofs_per_cell);
      dealii::Vector<double> local_mean_constraint(state_fe_->dofs_per_cell);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_->dofs_per_cell);

      for (auto cell = state_dof_handler_.begin_active();
           cell != state_dof_handler_.end();
           ++cell)
        {
          state_values.reinit(cell);
          local_system = 0.0;
          local_forcing = 0.0;
          local_mean_constraint = 0.0;
          for (unsigned int q = 0; q < volume_quadrature->size(); ++q)
            {
              const double weight = state_values.JxW(q);
              const auto &point = state_values.quadrature_point(q);
              const double forcing_value =
                forcing.value(point);
              for (unsigned int i = 0; i < state_fe_->dofs_per_cell; ++i)
                {
                  const double phi_i = state_values.shape_value(i, q);
                  local_forcing(i) += forcing_value * phi_i * weight;
                  local_mean_constraint(i) += phi_i * weight;
                  for (unsigned int j = 0; j < state_fe_->dofs_per_cell; ++j)
                    {
                      const double phi_j = state_values.shape_value(j, q);
                      local_system(i, j) +=
                        (diffusion_ * (state_values.shape_grad(i, q) *
                                       state_values.shape_grad(j, q)) +
                         reaction_ * phi_i * phi_j -
                         (conservative_transport == nullptr
                            ? 0.0
                            : phi_j *
                                (conservative_transport->value(point) *
                                 state_values.shape_grad(i, q)))) *
                        weight;
                    }
                }
            }
          cell->get_dof_indices(state_indices);
          for (unsigned int i = 0; i < state_fe_->dofs_per_cell; ++i)
            {
              const auto global_i = state_indices[i];
              if (constrained_state_dofs_.at(global_i))
                continue;
              forcing_load_[global_i] += local_forcing(i);
              mean_constraint_[global_i] += local_mean_constraint(i);
              for (unsigned int j = 0; j < state_fe_->dofs_per_cell; ++j)
                {
                  const auto global_j = state_indices[j];
                  if (constrained_state_dofs_.at(global_j))
                    fixed_state_load_[global_i] +=
                      local_system(i, j) * fixed_state_values_[global_j];
                  else
                    system_matrix_.add(global_i, global_j, local_system(i, j));
                }
            }
        }

      if (!uses_mean_zero_gauge())
        for (dealii::types::global_dof_index index = 0;
             index < state_dof_handler_.n_dofs();
             ++index)
          if (constrained_state_dofs_.at(index))
            system_matrix_.set(index, index, 1.0);

      assemble_boundary_operators(desired_state,
                                  observation_weight,
                                  quadrature_order,
                                  conservative_transport);
    }

    void
    assemble_boundary_operators(
      const dealii::Function<dim> &desired_state,
      const dealii::Function<dim> *observation_weight,
      const unsigned int           quadrature_order,
      const dealii::TensorFunction<1, dim> *conservative_transport)
    {
      const auto face_quadrature = make_gauss_face_quadrature(
        state_dof_handler_.get_triangulation(),
        quadrature_order,
        "The Neumann v1 target");
      dealii::UpdateFlags face_update_flags =
        dealii::update_values | dealii::update_quadrature_points |
        dealii::update_JxW_values;
      if (uses_ordinary_transport_boundary_realisation())
        face_update_flags |= dealii::update_normal_vectors;
      dealii::FEFaceValues<dim> face_values(
        *state_fe_,
        *face_quadrature,
        face_update_flags);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe_->dofs_per_cell);
      for (auto cell = state_dof_handler_.begin_active();
           cell != state_dof_handler_.end();
           ++cell)
        {
          cell->get_dof_indices(state_indices);
          for (unsigned int face = 0; face < cell->n_faces(); ++face)
            {
              const bool observation_face = is_observation_face(cell, face);
              const bool transport_face =
                uses_ordinary_transport_boundary_realisation() &&
                is_transport_boundary_face(cell, face);
              const bool natural_source_face =
                cell->face(face)->at_boundary() &&
                natural_boundary_source_ids_.count(
                  cell->face(face)->boundary_id()) != 0;
              if (!observation_face && !transport_face &&
                  !natural_source_face)
                continue;
              face_values.reinit(cell, face);
              for (unsigned int q = 0; q < face_quadrature->size(); ++q)
                {
                  const double weight = face_values.JxW(q);
                  const double boundary_matrix_coefficient =
                    transport_face
                      ? -(diffusion_ - 1.0) *
                        (conservative_transport->value(
                           face_values.quadrature_point(q)) *
                         face_values.normal_vector(q)) *
                        weight
                      : 0.0;
                  const double natural_source_value =
                    natural_source_face
                      ? boundary_source_scale() *
                          natural_boundary_source_->value(
                            face_values.quadrature_point(q))
                      : 0.0;
                  const double desired_value = observation_face
                    ? desired_state.value(face_values.quadrature_point(q))
                    : 0.0;
                  const double observation_weight_value =
                    observation_face && observation_weight != nullptr
                      ? observation_weight->value(
                          face_values.quadrature_point(q))
                      : 1.0;
                  std::optional<Vector> observation_evaluation;
                  if (observation_face)
                    {
                      desired_state_norm_ += desired_value * desired_value * weight;
                      observation_evaluation.emplace(state_dof_handler_.n_dofs());
                    }
                  for (unsigned int i = 0; i < state_fe_->dofs_per_cell; ++i)
                    {
                      const auto global_i = state_indices[i];
                      if (constrained_state_dofs_.at(global_i))
                        continue;
                      const double phi_i = face_values.shape_value(i, q);
                      if (transport_face)
                        for (unsigned int j = 0;
                             j < state_fe_->dofs_per_cell;
                             ++j)
                          {
                            const auto global_j = state_indices[j];
                            const double contribution =
                              boundary_matrix_coefficient * phi_i *
                              face_values.shape_value(j, q);
                            if (constrained_state_dofs_.at(global_j))
                              fixed_state_load_[global_i] +=
                                contribution * fixed_state_values_[global_j];
                            else
                              system_matrix_.add(global_i,
                                                 global_j,
                                                 contribution);
                          }
                      if (observation_face)
                        {
                          (*observation_evaluation)[global_i] =
                            observation_weight_value * phi_i;
                          desired_state_load_[global_i] +=
                            observation_weight_value * desired_value * phi_i *
                            weight;
                          for (unsigned int j = 0;
                               j < state_fe_->dofs_per_cell;
                               ++j)
                            {
                              const auto global_j = state_indices[j];
                              if (!constrained_state_dofs_.at(global_j))
                                state_tracking_matrix_.add(
                                  global_i,
                                  global_j,
                                  observation_weight_value *
                                    observation_weight_value * phi_i *
                                    face_values.shape_value(j, q) * weight);
                            }
                        }
                      if (natural_source_face)
                        forcing_load_[global_i] +=
                          natural_source_value * phi_i * weight;
                    }
                  if (observation_face)
                    {
                      observation_evaluations_.push_back(
                        std::move(*observation_evaluation));
                      observation_quadrature_weights_.push_back(weight);
                    }
                }
            }
        }
    }

    contract::LinearSolveReport
    solve_symmetric_system(
      Vector &                                      solution,
      const Vector &                                right_hand_side,
      const dealii_backend::SPDLinearSolvePolicy &policy) const
    {
      return dealii_backend::solve_serial_spd(system_matrix_,
                                              solution,
                                              right_hand_side,
                                              policy);
    }

    bool
    is_compatible(const Vector &load) const
    {
      const double tolerance = std::max(1e-12, 1e-11 * load.l2_norm());
      return std::abs(constant_mode_ * load) <= tolerance;
    }

    double
    boundary_source_scale() const
    {
      return uses_ordinary_transport_boundary_realisation() ? diffusion_ : 1.0;
    }

    void
    build_mean_zero_system()
    {
      const auto state_size = state_dof_handler_.n_dofs();
      dealii::DynamicSparsityPattern augmented_dsp(state_size + 1,
                                                    state_size + 1);
      for (dealii::types::global_dof_index row = 0; row < state_size; ++row)
        {
          for (auto entry = system_matrix_.begin(row);
               entry != system_matrix_.end(row);
               ++entry)
            augmented_dsp.add(row, entry->column());
          if (mean_constraint_[row] != 0.0)
            {
              augmented_dsp.add(row, state_size);
              augmented_dsp.add(state_size, row);
            }
        }
      augmented_sparsity_.copy_from(augmented_dsp);
      augmented_system_.reinit(augmented_sparsity_);
      for (dealii::types::global_dof_index row = 0; row < state_size; ++row)
        {
          for (auto entry = system_matrix_.begin(row);
               entry != system_matrix_.end(row);
               ++entry)
            augmented_system_.set(row, entry->column(), entry->value());
          augmented_system_.set(row, state_size, mean_constraint_[row]);
          augmented_system_.set(state_size, row, mean_constraint_[row]);
        }
      augmented_solver_.initialize(augmented_system_);
    }

    void
    solve_mean_zero_system(Vector &solution, const Vector &right_hand_side) const
    {
      const auto state_size = state_dof_handler_.n_dofs();
      dealii::Vector<double> augmented_right_hand_side(state_size + 1);
      for (dealii::types::global_dof_index index = 0; index < state_size; ++index)
        augmented_right_hand_side[index] = right_hand_side[index];
      dealii::Vector<double> augmented_solution(state_size + 1);
      augmented_solver_.vmult(augmented_solution, augmented_right_hand_side);
      for (dealii::types::global_dof_index index = 0; index < state_size; ++index)
        solution[index] = augmented_solution[index];
    }

    std::unique_ptr<dealii::FiniteElement<dim>> state_fe_;
    dealii::DoFHandler<dim> state_dof_handler_;
    dealii::AffineConstraints<double> state_constraints_;
    std::vector<bool> constrained_state_dofs_;

    const double diffusion_;
    const double reaction_;
    const double regularisation_weight_;
    const StateGauge state_gauge_;
    const std::set<dealii::types::boundary_id> dirichlet_boundary_ids_;
    const std::set<dealii::types::boundary_id> observation_boundary_ids_;
    const StateObservation state_observation_;
    const std::set<dealii::types::material_id> observation_material_ids_;
    const bool uses_conservative_transport_;
    const std::set<dealii::types::boundary_id> transport_boundary_ids_;
    const TransportBoundaryRealisation transport_boundary_realisation_;
    const std::optional<WeightedTraceRealisation> weighted_trace_realisation_;
    const bool has_fixed_dirichlet_data_;
    const dealii::Function<dim> *natural_boundary_source_;
    const std::set<dealii::types::boundary_id> natural_boundary_source_ids_;
    std::unique_ptr<NeumannControlRealisation<dim>> control_realisation_;
    std::unique_ptr<VolumeObservationAssembly<dim>> volume_observation_;
    std::vector<Vector> observation_evaluations_;
    std::vector<double> observation_quadrature_weights_;

    dealii::SparsityPattern state_sparsity_;
    dealii::SparsityPattern augmented_sparsity_;
    dealii::SparseMatrix<double> system_matrix_;
    dealii::SparseMatrix<double> state_tracking_matrix_;
    dealii::SparseMatrix<double> augmented_system_;
    dealii::SparseDirectUMFPACK augmented_solver_;
    std::unique_ptr<dealii::SparseDirectUMFPACK> nonsymmetric_solver_;
    Vector forcing_load_;
    Vector fixed_state_values_;
    Vector fixed_state_load_;
    Vector desired_state_load_;
    Vector mean_constraint_;
    Vector constant_mode_;
    double desired_state_norm_ = 0.0;

    contract::LayoutPtr variable_layout_;
    contract::LayoutPtr test_layout_;
    contract::LayoutPtr state_layout_;
    contract::LayoutPtr control_layout_;
  };
} // namespace nmopt::compiler::v1::detail
