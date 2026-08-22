#pragma once

#include "nmopt/application/chapter6.hpp"
#include "nmopt/application/runner.hpp"
#include "nmopt/compiler/v1/dealii_compiler.hpp"
#include "nmopt/experiment/reduced_envelope.hpp"

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::application::chapter6::dealii
{
  using Backend = nmopt::dealii_backend::SerialBackend;
  using Envelope =
    nmopt::experiment::ReducedSearchExperimentEnvelopeT<Backend>;
  using Evidence =
    nmopt::application::benchmark::BenchmarkExecutionEvidenceT<Envelope>;

  template <int dim>
  class B2ForcingFunction final : public ::dealii::Function<dim>
  {
  public:
    B2ForcingFunction(const B2ProblemParameters::ForcingSelection selection,
                      const double value)
      : value_(selection == B2ProblemParameters::ForcingSelection::zero ? 0.0 :
                                                               value)
    {}

    double
    value(const ::dealii::Point<dim> &,
          const unsigned int = 0) const override
    {
      return value_;
    }

  private:
    double value_;
  };

  template <int dim>
  class B2DesiredStateFunction final : public ::dealii::Function<dim>
  {
  public:
    static_assert(dim >= 2,
                  "The Chapter 6 B2 target requires at least two coordinates");

    explicit B2DesiredStateFunction(const GraetzCase graetz_case)
      : graetz_case_(graetz_case)
    {}

    double
    value(const ::dealii::Point<dim> &point,
          const unsigned int          component = 0) const override
    {
      (void)component;
      switch (graetz_case_)
        {
          case GraetzCase::observation_wings_constant_target:
          case GraetzCase::observation_full_constant_target:
            return 2.0;
          case GraetzCase::observation_wings_parabolic_target:
          case GraetzCase::observation_full_parabolic_target:
            return 4.0 * point[1] * (1.0 - point[1]);
        }
      throw std::invalid_argument("B2 has an unknown Graetz target case");
    }

  private:
    GraetzCase graetz_case_;
  };

  template <int dim>
  class B2ConservativeTransportFunction final
    : public ::dealii::TensorFunction<1, dim>
  {
  public:
    static_assert(dim >= 2,
                  "The Chapter 6 B2 transport requires at least two coordinates");

    ::dealii::Tensor<1, dim>
    value(const ::dealii::Point<dim> &point) const override
    {
      ::dealii::Tensor<1, dim> transport;
      transport[0] = 1.5 * point[1] * (1.0 - point[1]);
      return transport;
    }
  };

  template <int dim>
  struct B2RuntimeDataT
  {
    B2RuntimeDataT(
      const ::dealii::Function<dim> &             forcing_function,
      const ::dealii::Function<dim> &             desired_state_function,
      const ::dealii::Function<dim> &             fixed_temperature_function,
      const ::dealii::TensorFunction<1, dim> &transport_function)
      : forcing(forcing_function)
      , desired_state(desired_state_function)
      , fixed_temperature(fixed_temperature_function)
      , conservative_transport(transport_function)
    {}

    const ::dealii::Function<dim> &             forcing;
    const ::dealii::Function<dim> &             desired_state;
    const ::dealii::Function<dim> &             fixed_temperature;
    const ::dealii::TensorFunction<1, dim> &conservative_transport;
  };

  template <int dim>
  class B2ManufacturedDataT final
  {
  public:
    explicit B2ManufacturedDataT(
      const GraetzCase graetz_case =
        GraetzCase::observation_wings_constant_target,
      const double fixed_temperature = 1.0,
      const B2ProblemParameters::ForcingSelection forcing_selection =
        B2ProblemParameters::ForcingSelection::zero,
      const double forcing_value = 0.0)
      : forcing_(forcing_selection, forcing_value)
      , desired_state_(graetz_case)
      , fixed_temperature_(fixed_temperature)
      , fixed_temperature_value_(fixed_temperature)
      , graetz_case_(graetz_case)
      , forcing_selection_(forcing_selection)
      , forcing_value_(forcing_value)
    {}

    B2RuntimeDataT<dim>
    runtime_data() const
    {
      return {forcing_,
              desired_state_,
              fixed_temperature_,
              conservative_transport_};
    }

    double
    fixed_temperature_value() const
    {
      return fixed_temperature_value_;
    }

    GraetzCase
    graetz_case() const
    {
      return graetz_case_;
    }

    B2ProblemParameters::ForcingSelection
    forcing_selection() const
    {
      return forcing_selection_;
    }

    double
    forcing_value() const
    {
      return forcing_value_;
    }

  private:
    B2ForcingFunction<dim>                 forcing_;
    B2DesiredStateFunction<dim>             desired_state_;
    ::dealii::Functions::ConstantFunction<dim> fixed_temperature_;
    B2ConservativeTransportFunction<dim>    conservative_transport_;
    double                                  fixed_temperature_value_;
    GraetzCase                              graetz_case_;
    B2ProblemParameters::ForcingSelection  forcing_selection_;
    double                                  forcing_value_;
  };

  template <int dim>
  B2RuntimeDataT<dim>
  make_b2_manufactured_runtime_data(
    const B2Scenario &scenario,
    const B2ManufacturedDataT<dim> &data)
  {
    validate_b2(scenario);
    if (scenario.compile.mesh.dimension != dim)
      throw std::invalid_argument(
        "B2 runtime data dimension does not match the scenario mesh");
    if (data.graetz_case() != scenario.problem.graetz_case)
      throw std::invalid_argument(
        "B2 manufactured target case does not match the scenario");
    if (std::abs(data.fixed_temperature_value() -
                 scenario.problem.fixed_temperature) > 1.0e-14)
      throw std::invalid_argument(
        "B2 manufactured fixed temperature does not match the scenario");
    if (data.forcing_selection() != scenario.problem.forcing_selection ||
        std::abs(data.forcing_value() - scenario.problem.forcing_value) >
          1.0e-14)
      throw std::invalid_argument(
        "B2 manufactured forcing does not match the scenario");
    return data.runtime_data();
  }

  template <int dim>
  compiler::v1::DealiiDataBindings<dim>
  make_b2_data_bindings(const B2ProblemParameters & parameters,
                        const B2RuntimeDataT<dim> &runtime)
  {
    validate_runtime_data(parameters.data);
    if (parameters.data.fixed_dirichlet_data_provenance.empty() ||
        parameters.data.conservative_transport_provenance.empty())
      throw std::invalid_argument(
        "B2 runtime bindings need fixed-data and transport provenance");

    return {runtime.forcing,
            runtime.desired_state,
            std::optional<double>{parameters.data.diffusion},
            parameters.data.reaction,
            parameters.data.regularisation_weight,
            {parameters.data.forcing_provenance,
             parameters.data.desired_state_provenance,
             parameters.data.fixed_dirichlet_data_provenance},
            std::cref(runtime.fixed_temperature),
            std::nullopt,
            std::nullopt,
            compiler::v1::DealiiConservativeTransportDataBindings<dim>{
              runtime.conservative_transport,
              {parameters.data.conservative_transport_provenance}}};
  }

  template <int dim>
  std::shared_ptr<compiler::v1::DealiiCompilationSession<dim>>
  make_b2_compilation_session(const B2Scenario &scenario)
  {
    validate_b2(scenario);
    if (dim != 2 || scenario.compile.mesh.dimension != dim)
      throw std::invalid_argument(
        "B2 compilation sessions require a two-dimensional matching mesh");
    if (scenario.problem.recipe.observed_material_id == 0)
      throw std::invalid_argument(
        "B2 needs a nonzero observed material id");

    auto mesh = std::make_unique<::dealii::Triangulation<dim>>();
    ::dealii::Point<dim> lower;
    ::dealii::Point<dim> upper;
    upper[0] = 4.0;
    upper[1] = 1.0;
    switch (scenario.compile.mesh.generation)
      {
        case MeshGeneration::framework_native:
          ::dealii::GridGenerator::hyper_rectangle(*mesh, lower, upper);
          mesh->refine_global(scenario.compile.mesh.refinement);
          break;
        case MeshGeneration::structured_simplex:
          if constexpr (dim == 2)
            ::dealii::GridGenerator::
              subdivided_hyper_rectangle_with_simplices(
                *mesh,
                scenario.compile.mesh.axis_subdivisions,
                lower,
                upper,
                false);
          else
            throw std::invalid_argument(
              "B2 simplex meshes are implemented only in two dimensions");
          break;
        case MeshGeneration::centroid_split_simplex:
          throw std::invalid_argument(
            "B2 centroid-split simplex meshes are not implemented");
      }

    const auto is_close = [](const double first, const double second) {
      return std::abs(first - second) < 1.0e-12;
    };
    for (auto cell = mesh->begin_active(); cell != mesh->end(); ++cell)
      {
        const auto center = cell->center();
        const bool downstream = center[0] > 1.0;
        const bool wings = center[1] < 0.3 || center[1] > 0.7;
        const bool observed = downstream &&
                              (scenario.problem.graetz_case ==
                                 GraetzCase::observation_full_constant_target ||
                               scenario.problem.graetz_case ==
                                 GraetzCase::observation_full_parabolic_target ||
                               wings);
        cell->set_material_id(
          observed ? scenario.problem.recipe.observed_material_id : 0);
        for (unsigned int face = 0; face < cell->n_faces(); ++face)
          if (cell->face(face)->at_boundary())
            {
              const auto face_center = cell->face(face)->center();
              const bool left_boundary = is_close(face_center[0], 0.0);
              const bool right_boundary = is_close(face_center[0], 4.0);
              const bool horizontal_boundary =
                is_close(face_center[1], 0.0) ||
                is_close(face_center[1], 1.0);
              if (!left_boundary && !right_boundary && !horizontal_boundary)
                throw std::invalid_argument(
                  "B2 mesh contains an unclassified exterior boundary face");
              const auto boundary_id =
                left_boundary
                  ? chapter6::b2_fixed_boundary_id
                  : right_boundary
                    ? chapter6::b2_outflow_boundary_id
                    : face_center[0] < 1.0
                      ? chapter6::b2_fixed_boundary_id
                      : chapter6::b2_control_boundary_id;
              cell->face(face)->set_boundary_id(boundary_id);
            }
      }
    return std::make_shared<compiler::v1::DealiiCompilationSession<dim>>(
      std::move(mesh), scenario.compile.mesh.mesh_provenance);
  }

  inline compiler::v1::DealiiDiscretisationPolicy
  make_b2_discretisation_policy(const CompileOptions &options)
  {
    validate_common_compile_options(options);
    if (options.execution != ExecutionSelection::assembled)
      throw std::invalid_argument(
        "B2 runtime execution is limited to assembled operators");

    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = options.state_degree;
    policy.execution =
      compiler::v1::DealiiDiscretisationPolicy::Execution::assembled;
    policy.state_solve = {options.state_solve.maximum_iterations,
                          options.state_solve.relative_tolerance,
                          options.state_solve.absolute_tolerance};
    policy.adjoint_solve = {options.adjoint_solve.maximum_iterations,
                            options.adjoint_solve.relative_tolerance,
                            options.adjoint_solve.absolute_tolerance};
    policy.control_metric_solve = {
      options.control_metric_solve.maximum_iterations,
      options.control_metric_solve.relative_tolerance,
      options.control_metric_solve.absolute_tolerance};
    return policy;
  }

  inline compiler::v1::CompilationProduct
  make_b2_compilation_product(const CompileOptions &options)
  {
    if (options.product != ProductSelection::reduced_dto)
      throw std::invalid_argument(
        "B2 runtime execution is limited to the reduced DTO product");
    return compiler::v1::CompilationProduct::reduced_dto;
  }

  inline std::string
  b2_number(const double value)
  {
    std::ostringstream output;
    output << std::setprecision(17) << value;
    return output.str();
  }

  inline const char *
  b2_observation_region(const GraetzCase graetz_case)
  {
    switch (graetz_case)
      {
        case GraetzCase::observation_wings_constant_target:
        case GraetzCase::observation_wings_parabolic_target:
          return "wings";
        case GraetzCase::observation_full_constant_target:
        case GraetzCase::observation_full_parabolic_target:
          return "full";
        default:
          throw std::invalid_argument("B2 has an unknown observation region");
      }
  }

  inline const char *
  b2_target_profile(const GraetzCase graetz_case)
  {
    switch (graetz_case)
      {
        case GraetzCase::observation_wings_constant_target:
        case GraetzCase::observation_full_constant_target:
          return "constant";
        case GraetzCase::observation_wings_parabolic_target:
        case GraetzCase::observation_full_parabolic_target:
          return "parabolic";
        default:
          throw std::invalid_argument("B2 has an unknown target profile");
      }
  }

  inline double
  b2_relative_reduction(const double initial, const double final)
  {
    return (initial - final) /
           std::max(std::abs(initial), std::numeric_limits<double>::epsilon());
  }

  inline std::string
  b2_history(const std::vector<double> &history)
  {
    std::ostringstream output;
    output << std::setprecision(17);
    for (std::size_t index = 0; index < history.size(); ++index)
      {
        if (index != 0)
          output << ',';
        output << history[index];
      }
    return output.str();
  }

  template <typename Backend>
  double
  b2_block_values_l2_norm(const contract::BlockValuesT<Backend> &values)
  {
    double squared_norm = 0.0;
    for (std::size_t block = 0; block < values.n_blocks(); ++block)
      squared_norm += Backend::dot(values.block(block), values.block(block));
    return std::sqrt(std::max(0.0, squared_norm));
  }

  template <typename Backend>
  std::pair<double, double>
  b2_block_extrema(const contract::BlockValuesT<Backend> &values,
                   const std::size_t                      block = 0)
  {
    contract::require(block < values.n_blocks(),
                      "B2 extrema requested an absent block");
    contract::require(values.block(block).size() > 0,
                      "B2 extrema need a nonempty coefficient block");
    double minimum = values.block(block)[0];
    double maximum = minimum;
    for (std::size_t index = 1; index < values.block(block).size(); ++index)
      {
        minimum = std::min(minimum, values.block(block)[index]);
        maximum = std::max(maximum, values.block(block)[index]);
      }
    contract::require(std::isfinite(minimum) && std::isfinite(maximum),
                      "B2 extrema contain a non-finite coefficient");
    return {minimum, maximum};
  }

  template <typename Backend>
  struct B2DerivativeEvidenceT
  {
    double residual_jvp_error;
    double residual_vjp_error;
    double reduced_gradient_finite_difference_error;
    double reduced_taylor_step;
    double reduced_taylor_error;
    double reduced_taylor_order;
    double initial_objective;
    double initial_gradient_norm;
  };

  template <typename Backend>
  contract::PrimalBlockT<Backend>
  b2_shifted_point(const contract::PrimalBlockT<Backend> &point,
                   const contract::PrimalBlockT<Backend> &direction,
                   const double                            step)
  {
    contract::require(point.layout()->compatible_with(*direction.layout()),
                      "B2 evidence uses incompatible full-point directions");
    contract::PrimalBlockT<Backend> shifted = point;
    for (std::size_t block = 0; block < shifted.n_blocks(); ++block)
      shifted.add_scaled_block(block, step, direction.block(block));
    return shifted;
  }

  template <typename Backend>
  contract::PrimalBlockT<Backend>
  b2_control_shifted(const contract::PrimalBlockT<Backend> &control,
                     const contract::PrimalBlockT<Backend> &direction,
                     const double                            step)
  {
    contract::require(control.layout()->compatible_with(*direction.layout()),
                      "B2 evidence uses incompatible control directions");
    contract::PrimalBlockT<Backend> shifted = control;
    shifted.add_scaled_block(0, step, direction.block(0));
    return shifted;
  }

  template <typename Backend>
  B2DerivativeEvidenceT<Backend>
  make_b2_derivative_evidence(
    const compiler::v1::CompiledProblemT<Backend> &problem,
    const contract::ReducedDTOT<Backend> &         reduced,
    const contract::StateControlPartitionT<Backend> &partition,
    const contract::PrimalBlockT<Backend> &        initial_control)
  {
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;
    using Vector = typename Backend::Vector;

    const auto &model = problem.executable_model();
    const Primal zero_point = Primal::zeros(model.variable_layout());
    Vector state_direction_values(model.variable_layout()->dimension(0));
    Vector control_direction_values(model.variable_layout()->dimension(1));
    for (std::size_t index = 0;
         index < state_direction_values.size();
         ++index)
      state_direction_values[index] =
        (index % 2 == 0 ? 0.02 : -0.015) * static_cast<double>(index + 1);
    for (std::size_t index = 0;
         index < control_direction_values.size();
         ++index)
      control_direction_values[index] =
        (index % 2 == 0 ? -0.03 : 0.025) * static_cast<double>(index + 1);
    const Primal full_direction(model.variable_layout(),
                                {std::move(state_direction_values),
                                 std::move(control_direction_values)});

    Vector test_seed_values(model.test_layout()->dimension(0));
    for (std::size_t index = 0; index < test_seed_values.size(); ++index)
      test_seed_values[index] =
        (index % 3 == 0 ? 0.01 : -0.02) * static_cast<double>(index + 1);
    const Primal test_seed(model.test_layout(), {std::move(test_seed_values)});

    constexpr double residual_step = 1.0e-6;
    const Covector residual_zero = model.residual(zero_point);
    const Covector residual_plus = model.residual(
      b2_shifted_point(zero_point, full_direction, residual_step));
    const Covector residual_jvp =
      model.residual_jvp(zero_point, full_direction);
    Covector residual_finite_difference = residual_plus;
    residual_finite_difference.add_scaled_block(
      0, -1.0, residual_zero.block(0));
    residual_finite_difference.scale_block(0, 1.0 / residual_step);
    residual_finite_difference.add_scaled_block(
      0, -1.0, residual_jvp.block(0));
    const double residual_jvp_error =
      b2_block_values_l2_norm(residual_finite_difference);

    const Covector residual_transpose =
      model.residual_vjp(zero_point, test_seed);
    const double residual_vjp_error = std::abs(
      contract::pair(residual_jvp, test_seed) -
      contract::pair(residual_transpose, full_direction));

    Vector reduced_direction_values(partition.control_layout()->dimension(0));
    for (std::size_t index = 0; index < reduced_direction_values.size(); ++index)
      reduced_direction_values[index] =
        (index % 2 == 0 ? 0.03 : -0.02) * static_cast<double>(index + 1);
    const Primal reduced_direction(partition.control_layout(),
                                   {std::move(reduced_direction_values)});

    constexpr double taylor_step = 1.0e-4;
    const auto base_evaluation = reduced.evaluate(initial_control);
    const auto plus_evaluation = reduced.evaluate(b2_control_shifted(
      initial_control, reduced_direction, taylor_step));
    const auto minus_evaluation = reduced.evaluate(b2_control_shifted(
      initial_control, reduced_direction, -taylor_step));
    const auto half_plus_evaluation = reduced.evaluate(b2_control_shifted(
      initial_control, reduced_direction, 0.5 * taylor_step));
    const double directional_derivative = contract::pair(
      base_evaluation.reduced_derivative, reduced_direction);
    const double central_derivative =
      (plus_evaluation.objective_value - minus_evaluation.objective_value) /
      (2.0 * taylor_step);
    const double reduced_gradient_finite_difference_error =
      std::abs(central_derivative - directional_derivative);
    const double taylor_error = std::abs(
      plus_evaluation.objective_value - base_evaluation.objective_value -
      taylor_step * directional_derivative);
    const double half_taylor_error = std::abs(
      half_plus_evaluation.objective_value - base_evaluation.objective_value -
      0.5 * taylor_step * directional_derivative);
    const double taylor_order =
      taylor_error > std::numeric_limits<double>::epsilon() &&
          half_taylor_error > std::numeric_limits<double>::epsilon()
        ? std::log(taylor_error / half_taylor_error) / std::log(2.0)
        : std::numeric_limits<double>::infinity();

    contract::require(std::isfinite(residual_jvp_error) &&
                        std::isfinite(residual_vjp_error) &&
                        std::isfinite(reduced_gradient_finite_difference_error) &&
                        std::isfinite(taylor_error) &&
                        std::isfinite(taylor_order),
                      "B2 derivative evidence produced a non-finite value");
    contract::require(residual_jvp_error <= 1.0e-7,
                      "B2 residual JVP finite-difference evidence failed");
    contract::require(residual_vjp_error <= 1.0e-7,
                      "B2 residual VJP transpose evidence failed");
    contract::require(reduced_gradient_finite_difference_error <= 1.0e-7,
                      "B2 reduced gradient finite-difference evidence failed");
    contract::require(taylor_order >= 1.5,
                      "B2 reduced Taylor evidence failed");

    return {residual_jvp_error,
            residual_vjp_error,
            reduced_gradient_finite_difference_error,
            taylor_step,
            taylor_error,
            taylor_order,
            base_evaluation.objective_value,
            b2_block_values_l2_norm(base_evaluation.reduced_derivative)};
  }

  template <int dim>
  class B2ReducedExecutionAdapterT final
  {
  public:
    B2ReducedExecutionAdapterT(
      const B2RuntimeDataT<dim> &runtime,
      std::shared_ptr<compiler::v1::DealiiCompilationSession<dim>> session,
      experiment::RunEnvironmentRecord environment,
      std::filesystem::path          native_output_directory = {})
      : runtime_(&runtime)
      , session_(std::move(session))
      , environment_(std::move(environment))
      , native_output_directory_(std::move(native_output_directory))
    {
      if (runtime_ == nullptr || !session_)
        throw std::invalid_argument(
          "B2 execution needs runtime data and an owned compilation session");
    }

    Evidence
    operator()(const semantic::v1::ProblemSpec &specification,
               const B2Scenario &                 scenario) const
    {
      validate_b2(scenario);
      if (specification.id != "scalar_convection_neumann_subdomain_control")
        throw std::invalid_argument(
          "B2 execution received a non-Graetz scalar ProblemSpec");
      if (scenario.solver.method != ReducedMethod::bfgs)
        throw std::invalid_argument("B2 execution requires full BFGS");

      const auto policy = make_b2_discretisation_policy(scenario.compile);
      const auto product = make_b2_compilation_product(scenario.compile);
      const auto bindings = make_b2_data_bindings(scenario.problem, *runtime_);
      compiler::v1::DealiiCompiler compiler;
      auto compilation = compiler.compile(specification,
                                          session_,
                                          bindings,
                                          policy,
                                          std::nullopt,
                                          std::nullopt,
                                          product);
      if (!compilation.succeeded() || !compilation.problem)
        {
          std::ostringstream message;
          message << "B2 compilation failed for the selected recipe and bindings";
          for (const auto &diagnostic : compilation.diagnostics.diagnostics())
            message << " [" << diagnostic.component_id << ":"
                    << diagnostic.capability << ":" << diagnostic.remedy << "]";
          throw std::runtime_error(message.str());
        }

      using Model = compiler::v1::detail::NeumannBoundaryControlModel<dim>;
      const auto *model =
        dynamic_cast<const Model *>(&compilation.problem->executable_model());
      if (model == nullptr)
        throw std::runtime_error(
          "B2 execution needs the Neumann boundary model");

      const auto reduced = compilation.problem->make_reduced_dto();
      const contract::StateControlPartitionT<Backend> partition(
        compilation.problem->executable_model(), 0, 1);
      const auto initial_control =
        contract::PrimalBlockT<Backend>::zeros(partition.control_layout());
      const auto initial_evaluation = reduced.evaluate(initial_control);
      const auto derivative_evidence = make_b2_derivative_evidence(
        *compilation.problem, reduced, partition, initial_control);
      const auto report =
        solvers::ReducedFullBfgsSolverT<Backend>(
          reduced,
          compilation.problem->metric(),
          scenario.solver.parameters)
          .solve(initial_control);
      if (scenario.experiment.retain_fields &&
          !native_output_directory_.empty())
        {
          model->write_native_output(native_output_directory_,
                                     report.final_evaluation.state,
                                     report.control,
                                     report.final_evaluation.adjoint,
                                     &initial_evaluation.state,
                                     &runtime_->forcing,
                                     &runtime_->desired_state);
        }
      contract::require(!report.objective_history.empty() &&
                          !report.gradient_norm_history.empty(),
                        "B2 solver did not retain objective and gradient histories");
      const double initial_objective = report.objective_history.front();
      const double final_objective = report.objective_history.back();
      const double initial_gradient_norm =
        report.gradient_norm_history.front();
      const double final_gradient_norm = report.gradient_norm_history.back();

      const auto initial_objective_components =
        model->objective_components(initial_evaluation.full_point);
      const auto final_objective_components =
        model->objective_components(report.final_evaluation.full_point);
      contract::require(
        std::abs(initial_objective_components.state_tracking +
                   initial_objective_components.control_regularisation -
                 initial_objective) <=
          1.0e-12 * std::max(1.0, std::abs(initial_objective)) &&
          std::abs(final_objective_components.state_tracking +
                     final_objective_components.control_regularisation -
                   final_objective) <=
            1.0e-12 * std::max(1.0, std::abs(final_objective)),
        "B2 objective components do not sum to the reported objective");

      std::vector<double> coefficient_derivative_norm_history{
        b2_block_values_l2_norm(initial_evaluation.reduced_derivative)};
      coefficient_derivative_norm_history.reserve(
        report.iteration_records.size() + 1);
      for (const auto &record : report.iteration_records)
        coefficient_derivative_norm_history.push_back(
          b2_block_values_l2_norm(
            record.common.accepted_evaluation.reduced_derivative));
      contract::require(
        coefficient_derivative_norm_history.size() ==
          report.objective_history.size(),
        "B2 coefficient-derivative and objective histories disagree");

      const auto uncontrolled_state_extrema =
        b2_block_extrema(initial_evaluation.state);
      const auto optimized_state_extrema =
        b2_block_extrema(report.final_evaluation.state);
      const auto control_extrema = b2_block_extrema(report.control);
      const auto adjoint_extrema =
        b2_block_extrema(report.final_evaluation.adjoint);

      const auto &triangulation = session_->triangulation();
      std::size_t boundary_face_count = 0;
      std::size_t fixed_boundary_face_count = 0;
      std::size_t control_boundary_face_count = 0;
      std::size_t outflow_boundary_face_count = 0;
      double      observation_measure = 0.0;
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->material_id() ==
              scenario.problem.recipe.observed_material_id)
            observation_measure += cell->measure();
          for (unsigned int face = 0; face < cell->n_faces(); ++face)
            if (cell->face(face)->at_boundary())
              {
                ++boundary_face_count;
                switch (cell->face(face)->boundary_id())
                  {
                    case chapter6::b2_fixed_boundary_id:
                      ++fixed_boundary_face_count;
                      break;
                    case chapter6::b2_control_boundary_id:
                      ++control_boundary_face_count;
                      break;
                    case chapter6::b2_outflow_boundary_id:
                      ++outflow_boundary_face_count;
                      break;
                    default:
                      throw std::runtime_error(
                        "B2 evidence found an unclassified boundary face");
                  }
              }
        }
      contract::require(
        boundary_face_count == fixed_boundary_face_count +
                                 control_boundary_face_count +
                                 outflow_boundary_face_count &&
          observation_measure > 0.0,
        "B2 structural evidence is incomplete");
      const auto solver_policy =
        experiment::make_reduced_search_policy_snapshot(report);

      std::vector<std::string> selected_fields{
        "objective_history",
        "gradient_norm_history",
        "step_length_history",
        "objective_change_history",
        "solve_counts"};
      if (scenario.experiment.retain_fields)
        {
          selected_fields.insert(selected_fields.end(),
                                 {"state",
                                  "state_uncontrolled",
                                  "control",
                                  "adjoint",
                                  "negative_adjoint",
                                  "target",
                                  "forcing",
                                  "observation_region"});
        }

      std::vector<benchmark::ArtifactField> fields{
        {"b2.graetz_case", graetz_case_name(scenario.problem.graetz_case)},
        {"b2.observation_region",
         b2_observation_region(scenario.problem.graetz_case)},
        {"b2.target_profile", b2_target_profile(scenario.problem.graetz_case)},
        {"b2.fixed_temperature", b2_number(scenario.problem.fixed_temperature)},
        {"b2.regularisation_weight",
         b2_number(scenario.problem.data.regularisation_weight)},
        {"benchmark.mesh_vertices",
         std::to_string(triangulation.n_vertices())},
        {"benchmark.mesh_active_cells",
         std::to_string(triangulation.n_active_cells())},
        {"benchmark.boundary_face_count",
         std::to_string(boundary_face_count)},
        {"benchmark.fixed_boundary_face_count",
         std::to_string(fixed_boundary_face_count)},
        {"benchmark.control_boundary_face_count",
         std::to_string(control_boundary_face_count)},
        {"benchmark.outflow_boundary_face_count",
         std::to_string(outflow_boundary_face_count)},
        {"benchmark.observation_measure", b2_number(observation_measure)},
        {"benchmark.state_dimension",
         std::to_string(report.final_evaluation.state.block(0).size())},
        {"benchmark.state_physical_dimension",
         std::to_string(model->physical_state_dimension())},
        {"benchmark.state_independent_dimension",
         std::to_string(model->independent_state_dimension())},
        {"benchmark.control_dimension",
         std::to_string(report.control.block(0).size())},
        {"benchmark.control_physical_dimension",
         std::to_string(model->physical_control_dimension())},
        {"benchmark.control_independent_dimension",
         std::to_string(model->independent_control_dimension())},
        {"benchmark.adjoint_dimension",
         std::to_string(report.final_evaluation.adjoint.block(0).size())},
        {"benchmark.adjoint_physical_dimension",
         std::to_string(model->physical_state_dimension())},
        {"benchmark.adjoint_independent_dimension",
         std::to_string(model->independent_state_dimension())},
        {"b2.derivative_evidence", "finite_difference_and_taylor"},
        {"b2.residual_jvp_error",
         b2_number(derivative_evidence.residual_jvp_error)},
        {"b2.residual_vjp_error",
         b2_number(derivative_evidence.residual_vjp_error)},
        {"b2.reduced_gradient_finite_difference_error",
         b2_number(derivative_evidence.reduced_gradient_finite_difference_error)},
        {"b2.reduced_taylor_step",
         b2_number(derivative_evidence.reduced_taylor_step)},
        {"b2.reduced_taylor_error",
         b2_number(derivative_evidence.reduced_taylor_error)},
        {"b2.reduced_taylor_order",
         b2_number(derivative_evidence.reduced_taylor_order)},
        {"b2.derivative_evidence_passed", "true"},
        {"b2.initial_objective", b2_number(initial_objective)},
        {"b2.initial_tracking_objective",
         b2_number(initial_objective_components.state_tracking)},
        {"b2.initial_control_regularisation_objective",
         b2_number(initial_objective_components.control_regularisation)},
        {"b2.final_objective", b2_number(final_objective)},
        {"b2.final_tracking_objective",
         b2_number(final_objective_components.state_tracking)},
        {"b2.final_control_regularisation_objective",
         b2_number(final_objective_components.control_regularisation)},
        {"b2.relative_objective_reduction",
         b2_number(b2_relative_reduction(initial_objective, final_objective))},
        {"b2.initial_gradient_norm", b2_number(initial_gradient_norm)},
        {"b2.initial_metric_gradient_norm",
         b2_number(initial_gradient_norm)},
        {"b2.final_gradient_norm", b2_number(final_gradient_norm)},
        {"b2.final_metric_gradient_norm", b2_number(final_gradient_norm)},
        {"b2.relative_gradient_reduction",
         b2_number(b2_relative_reduction(initial_gradient_norm,
                                         final_gradient_norm))},
        {"b2.initial_coefficient_derivative_norm",
         b2_number(coefficient_derivative_norm_history.front())},
        {"b2.final_coefficient_derivative_norm",
         b2_number(coefficient_derivative_norm_history.back())},
        {"b2.relative_coefficient_derivative_reduction",
         b2_number(b2_relative_reduction(
           coefficient_derivative_norm_history.front(),
           coefficient_derivative_norm_history.back()))},
        {"b2.state_l2_norm",
         b2_number(b2_block_values_l2_norm(report.final_evaluation.state))},
        {"b2.control_l2_norm", b2_number(b2_block_values_l2_norm(report.control))},
        {"b2.uncontrolled_state_min",
         b2_number(uncontrolled_state_extrema.first)},
        {"b2.uncontrolled_state_max",
         b2_number(uncontrolled_state_extrema.second)},
        {"b2.optimized_state_min", b2_number(optimized_state_extrema.first)},
        {"b2.optimized_state_max", b2_number(optimized_state_extrema.second)},
        {"b2.control_min", b2_number(control_extrema.first)},
        {"b2.control_max", b2_number(control_extrema.second)},
        {"b2.adjoint_min", b2_number(adjoint_extrema.first)},
        {"b2.adjoint_max", b2_number(adjoint_extrema.second)},
        {"solver.method", "bfgs"},
        {"solver.initial_control", scenario.solver.initial_control},
        {"solver.policy", report.policy_name},
        {"solver.stopping_reason",
         solvers::reduced_stopping_reason_name(report.stopping_reason)},
        {"solver.accepted_iterations",
         std::to_string(report.accepted_iterations)},
        {"solver.line_search_trial_count",
         std::to_string(report.line_search_trial_count)},
        {"solver.state_solve_count", std::to_string(report.state_solve_count)},
        {"solver.adjoint_solve_count",
         std::to_string(report.adjoint_solve_count)},
        {"solver.metric_solve_count",
         std::to_string(report.metric_solve_count)},
        {"solver.hessian_action_count",
         std::to_string(report.hessian_action_count)},
        {"solver.direction_reset_count",
         std::to_string(report.direction_reset_count)},
        {"solver.objective_history", b2_history(report.objective_history)},
        {"solver.gradient_norm_history",
         b2_history(report.gradient_norm_history)},
        {"solver.metric_gradient_norm_history",
         b2_history(report.gradient_norm_history)},
        {"solver.coefficient_derivative_norm_history",
         b2_history(coefficient_derivative_norm_history)},
        {"solver.step_length_history", b2_history(report.step_length_history)},
        {"solver.step_norm_history", b2_history(report.step_norm_history)},
        {"solver.objective_change_history",
         b2_history(report.objective_change_history)}};

      Envelope envelope{compilation.problem->manifest(),
                        solver_policy,
              report,
              environment_};
      return {std::move(envelope),
              std::move(compilation.diagnostics),
              {},
              std::move(selected_fields),
              std::move(fields)};
    }

  private:
    const B2RuntimeDataT<dim> *runtime_;
    std::shared_ptr<compiler::v1::DealiiCompilationSession<dim>> session_;
    experiment::RunEnvironmentRecord environment_;
    std::filesystem::path native_output_directory_;
  };
} // namespace nmopt::application::chapter6::dealii
