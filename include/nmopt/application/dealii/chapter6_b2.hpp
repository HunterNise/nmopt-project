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
      const double fixed_temperature = 1.0)
      : desired_state_(graetz_case)
      , fixed_temperature_(fixed_temperature)
      , fixed_temperature_value_(fixed_temperature)
      , graetz_case_(graetz_case)
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

  private:
    ::dealii::Functions::ZeroFunction<dim> forcing_;
    B2DesiredStateFunction<dim>             desired_state_;
    ::dealii::Functions::ConstantFunction<dim> fixed_temperature_;
    B2ConservativeTransportFunction<dim>    conservative_transport_;
    double                                  fixed_temperature_value_;
    GraetzCase                              graetz_case_;
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
    ::dealii::GridGenerator::hyper_rectangle(*mesh, lower, upper);
    mesh->refine_global(scenario.compile.mesh.refinement);

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
        for (unsigned int face = 0;
             face < ::dealii::GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary())
            {
              const auto face_center = cell->face(face)->center();
              const auto boundary_id =
                is_close(face_center[0], 0.0)
                  ? 0
                  : 1;
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

  template <int dim>
  class B2ReducedExecutionAdapterT final
  {
  public:
    B2ReducedExecutionAdapterT(
      const B2RuntimeDataT<dim> &runtime,
      std::shared_ptr<compiler::v1::DealiiCompilationSession<dim>> session,
      experiment::RunEnvironmentRecord environment,
      std::filesystem::path          field_output_directory = {})
      : runtime_(&runtime)
      , session_(std::move(session))
      , environment_(std::move(environment))
      , field_output_directory_(std::move(field_output_directory))
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

      const auto reduced = compilation.problem->make_reduced_dto();
      const contract::StateControlPartitionT<Backend> partition(
        compilation.problem->executable_model(), 0, 1);
      const auto initial_control =
        contract::PrimalBlockT<Backend>::zeros(partition.control_layout());
      const auto report =
        solvers::ReducedFullBfgsSolverT<Backend>(
          reduced,
          compilation.problem->metric(),
          scenario.solver.parameters)
          .solve(initial_control);
      if (scenario.experiment.retain_fields &&
          !field_output_directory_.empty())
        {
          const auto *model = dynamic_cast<const
            compiler::v1::detail::NeumannBoundaryControlModel<dim> *>(
            &compilation.problem->executable_model());
          if (model == nullptr)
            throw std::runtime_error(
              "B2 field output needs the Neumann boundary model");
          model->write_field_output(field_output_directory_,
                                    report.final_evaluation.state,
                                    report.control,
                                    report.final_evaluation.adjoint);
        }
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
                                 {"state", "control", "adjoint"});
        }

      std::vector<benchmark::ArtifactField> fields{
        {"b2.graetz_case", graetz_case_name(scenario.problem.graetz_case)},
        {"b2.fixed_temperature", b2_number(scenario.problem.fixed_temperature)},
        {"b2.regularisation_weight",
         b2_number(scenario.problem.data.regularisation_weight)},
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
    std::filesystem::path field_output_directory_;
  };
} // namespace nmopt::application::chapter6::dealii
