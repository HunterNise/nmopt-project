#pragma once

#include "nmopt/application/chapter6.hpp"
#include "nmopt/application/runner.hpp"
#include "nmopt/compiler/v1/dealii_continuous_control.hpp"
#include "nmopt/compiler/v1/dealii_compiler.hpp"
#include "nmopt/experiment/reduced_envelope.hpp"

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
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
  class B1DesiredStateFunction final : public ::dealii::Function<dim>
  {
  public:
    static_assert(dim >= 2,
                  "The Chapter 6 B1 target requires at least two coordinates");

    double
    value(const ::dealii::Point<dim> &point,
          const unsigned int          component = 0) const override
    {
      (void)component;
      return 10.0 * point[0] * (1.0 - point[0]) * point[1] *
             (1.0 - point[1]);
    }
  };

  template <int dim>
  struct B1RuntimeDataT
  {
    B1RuntimeDataT(const ::dealii::Function<dim> &forcing_function,
                   const ::dealii::Function<dim> &desired_state_function)
      : forcing(forcing_function)
      , desired_state(desired_state_function)
    {}

    const ::dealii::Function<dim> &forcing;
    const ::dealii::Function<dim> &desired_state;
  };

  template <int dim>
  class B1SelectedDataT final
  {
  public:
    explicit B1SelectedDataT(const B1ForcingSelection forcing_selection =
                               B1ForcingSelection::manufactured_zero)
      : forcing_selection_(forcing_selection)
      , forcing_(forcing_value(forcing_selection))
    {}

    B1ForcingSelection
    forcing_selection() const
    {
      return forcing_selection_;
    }

    B1RuntimeDataT<dim>
    runtime_data() const
    {
      return {forcing_, desired_state_};
    }

  private:
    static double
    forcing_value(const B1ForcingSelection selection)
    {
      switch (selection)
        {
          case B1ForcingSelection::manufactured_zero:
            return 0.0;
          case B1ForcingSelection::figure_inferred_constant_one:
            return 1.0;
        }
      throw std::invalid_argument("unknown B1 forcing selection");
    }

    B1ForcingSelection                        forcing_selection_;
    ::dealii::Functions::ConstantFunction<dim> forcing_;
    B1DesiredStateFunction<dim>                desired_state_;
  };

  template <int dim>
  B1RuntimeDataT<dim>
  make_b1_runtime_data(
    const B1Scenario &scenario,
    const B1SelectedDataT<dim> &data)
  {
    validate_b1(scenario);
    if (scenario.problem.forcing_selection != data.forcing_selection())
      throw std::invalid_argument(
        "B1 runtime-data forcing does not match the scenario selection");
    if (scenario.compile.mesh.dimension != dim)
      throw std::invalid_argument(
        "B1 runtime data dimension does not match the scenario mesh");
    return data.runtime_data();
  }

  template <int dim>
  compiler::v1::DealiiDataBindings<dim>
  make_b1_data_bindings(const B1ProblemParameters &parameters,
                        const double               regularisation,
                        const B1RuntimeDataT<dim> & runtime)
  {
    if (!std::isfinite(regularisation) || regularisation <= 0.0)
      throw std::invalid_argument(
        "B1 runtime regularisation must be positive and finite");
    if (!std::isfinite(parameters.data.diffusion) ||
        parameters.data.diffusion <= 0.0)
      throw std::invalid_argument("B1 runtime diffusion must be positive");
    if (!std::isfinite(parameters.data.reaction) ||
        parameters.data.reaction < 0.0)
      throw std::invalid_argument("B1 runtime reaction must be nonnegative");
    if (parameters.data.forcing_provenance.empty() ||
        parameters.data.desired_state_provenance.empty())
      throw std::invalid_argument(
        "B1 runtime bindings need forcing and target provenance");

    return {runtime.forcing,
            runtime.desired_state,
            std::optional<double>{parameters.data.diffusion},
            parameters.data.reaction,
            regularisation,
            {parameters.data.forcing_provenance,
             parameters.data.desired_state_provenance,
             ""}};
  }

  template <int dim>
  std::shared_ptr<compiler::v1::DealiiCompilationSession<dim>>
  make_b1_compilation_session(const B1Scenario &scenario)
  {
    validate_b1(scenario);
    if (scenario.compile.mesh.dimension != dim)
      throw std::invalid_argument(
        "B1 compilation session dimension does not match the scenario mesh");

    auto mesh = std::make_unique<::dealii::Triangulation<dim>>();
    ::dealii::GridGenerator::hyper_cube(*mesh, 0.0, 1.0);
    mesh->refine_global(scenario.compile.mesh.refinement);
    return std::make_shared<compiler::v1::DealiiCompilationSession<dim>>(
      std::move(mesh), scenario.compile.mesh.mesh_provenance);
  }

  inline compiler::v1::DealiiDiscretisationPolicy
  make_b1_discretisation_policy(const CompileOptions &options)
  {
    validate_common_compile_options(options);
    if (options.execution != ExecutionSelection::assembled)
      throw std::invalid_argument(
        "B1 runtime execution is limited to assembled operators");

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
  make_b1_compilation_product(const CompileOptions &options)
  {
    if (options.product != ProductSelection::reduced_dto)
      throw std::invalid_argument(
        "B1 runtime execution is limited to the reduced DTO product");
    return compiler::v1::CompilationProduct::reduced_dto;
  }

  inline const char *
  b1_method_name(const ReducedMethod method)
  {
    switch (method)
      {
        case ReducedMethod::steepest_descent:
          return "steepest_descent";
        case ReducedMethod::limited_memory_bfgs:
          return "limited_memory_bfgs";
        case ReducedMethod::bfgs:
          return "bfgs";
      }
    return "unknown";
  }

  inline std::string
  b1_number(const double value)
  {
    std::ostringstream output;
    output << std::setprecision(17) << value;
    return output.str();
  }

  inline std::string
  b1_history(const std::vector<double> &history)
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
  struct B1HessianEvidenceT
  {
    double direction_norm;
    double symmetry_error;
    double finite_difference_step;
    double finite_difference_error;
  };

  template <typename Backend>
  contract::PrimalBlockT<Backend>
  b1_shifted_control(const contract::PrimalBlockT<Backend> &control,
                     const contract::PrimalBlockT<Backend> &direction,
                     const double                            step)
  {
    contract::require(control.layout()->compatible_with(*direction.layout()),
                      "B1 Hessian evidence uses incompatible control directions");
    contract::PrimalBlockT<Backend> shifted = control;
    shifted.add_scaled_block(0, step, direction.block(0));
    return shifted;
  }

  template <typename Backend>
  B1HessianEvidenceT<Backend>
  make_b1_hessian_evidence(
    const compiler::v1::CompiledProblemT<Backend> &problem,
    const contract::PrimalBlockT<Backend> &        control)
  {
    const auto *hessian = problem.reduced_hessian();
    contract::require(hessian != nullptr,
                      "B1 compiled target omitted its reduced Hessian");

    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;
    using Vector = typename Backend::Vector;

    Vector direction_values(control.layout()->dimension(0));
    Vector second_direction_values(control.layout()->dimension(0));
    for (std::size_t index = 0; index < direction_values.size(); ++index)
      {
        direction_values[index] =
          (index % 2 == 0 ? 0.04 : -0.03) * static_cast<double>(index + 1);
        second_direction_values[index] =
          (index % 3 == 0 ? -0.02 : 0.03) * static_cast<double>(index + 1);
      }
    const Primal direction(control.layout(), {std::move(direction_values)});
    const Primal second_direction(
      control.layout(), {std::move(second_direction_values)});

    const Covector hessian_action = hessian->apply(control, direction);
    const Covector second_hessian_action =
      hessian->apply(control, second_direction);
    const double symmetry_error = std::abs(
      contract::pair(hessian_action, second_direction) -
      contract::pair(second_hessian_action, direction));

    constexpr double finite_difference_step = 1.0e-5;
    const auto reduced = problem.make_reduced_dto();
    const Covector reduced_derivative_plus =
      reduced.evaluate(b1_shifted_control(
                         control, direction, finite_difference_step))
        .reduced_derivative;
    const Covector reduced_derivative_minus =
      reduced.evaluate(b1_shifted_control(
                         control, direction, -finite_difference_step))
        .reduced_derivative;
    Covector finite_difference = reduced_derivative_plus;
    finite_difference.add_scaled_block(
      0, -1.0, reduced_derivative_minus.block(0));
    finite_difference.scale_block(0, 1.0 / (2.0 * finite_difference_step));
    finite_difference.add_scaled_block(0, -1.0, hessian_action.block(0));
    const double finite_difference_error =
      finite_difference.block(0).l2_norm();

    contract::require(std::isfinite(symmetry_error) && symmetry_error <= 1.0e-10,
                      "B1 reduced Hessian symmetry evidence failed");
    contract::require(std::isfinite(finite_difference_error) &&
                        finite_difference_error <= 2.0e-8,
                      "B1 reduced Hessian finite-difference evidence failed");

    return {direction.block(0).l2_norm(),
            symmetry_error,
            finite_difference_step,
            finite_difference_error};
  }

  template <int dim>
  class B1ReducedExecutionAdapterT final
  {
  public:
    B1ReducedExecutionAdapterT(
      const double                  regularisation,
      const B1RuntimeDataT<dim> &runtime,
      std::shared_ptr<compiler::v1::DealiiCompilationSession<dim>> session,
      experiment::RunEnvironmentRecord environment,
      std::filesystem::path          native_output_directory = {})
      : regularisation_(regularisation)
      , runtime_(&runtime)
      , session_(std::move(session))
      , environment_(std::move(environment))
      , native_output_directory_(std::move(native_output_directory))
    {
      if (runtime_ == nullptr || !session_)
        throw std::invalid_argument(
          "B1 execution needs runtime data and an owned compilation session");
      if (!std::isfinite(regularisation_) || regularisation_ <= 0.0)
        throw std::invalid_argument(
          "B1 execution regularisation must be positive and finite");
    }

    Evidence
    operator()(const semantic::v1::ProblemSpec &specification,
               const B1Scenario &                 scenario) const
    {
      validate_b1(scenario);
      const auto expected_specification_id = [&scenario]() -> const char * {
        switch (scenario.problem.recipe.discretisation)
          {
            case chapter5::DistributedControlDiscretisation::cellwise_constant:
              return "scalar_diffusion_reaction_volume_control";
            case chapter5::DistributedControlDiscretisation::
              homogeneous_dirichlet_continuous:
              return "scalar_diffusion_reaction_l2_state_tracking_continuous_control";
          }
        throw std::invalid_argument(
          "B1 execution has an unknown control discretisation");
      }();
      if (specification.id != expected_specification_id)
        throw std::invalid_argument(
          "B1 execution received a non-distributed scalar ProblemSpec");
      if (std::find(scenario.problem.regularisation_sweep.begin(),
                    scenario.problem.regularisation_sweep.end(),
                    regularisation_) ==
          scenario.problem.regularisation_sweep.end())
        throw std::invalid_argument(
          "B1 execution regularisation is not in the frozen scenario sweep");

      const auto policy = make_b1_discretisation_policy(scenario.compile);
      const auto product = make_b1_compilation_product(scenario.compile);
      const auto bindings = make_b1_data_bindings(
        scenario.problem, regularisation_, *runtime_);
      compiler::v1::DealiiCompiler compiler;
      auto compilation = compiler.compile(specification,
                                          session_,
                                          bindings,
                                          policy,
                                          std::nullopt,
                                          std::nullopt,
                                          product);
      if (!compilation.succeeded() || !compilation.problem)
        throw std::runtime_error(
          "B1 compilation failed for the selected recipe and bindings");

      const auto reduced = compilation.problem->make_reduced_dto();
      const contract::StateControlPartitionT<Backend> partition(
        compilation.problem->executable_model(), 0, 1);
      const auto initial_control =
        contract::PrimalBlockT<Backend>::zeros(partition.control_layout());
      const auto hessian_evidence =
        make_b1_hessian_evidence(*compilation.problem, initial_control);

      std::optional<solvers::ReducedSolverResultT<Backend>> report;
      switch (scenario.solver.method)
        {
          case ReducedMethod::steepest_descent:
            report.emplace(
              solvers::ReducedGradientSolverT<Backend>(
                reduced,
                compilation.problem->metric(),
                scenario.solver.parameters)
                .solve(initial_control));
            break;
          case ReducedMethod::limited_memory_bfgs:
            report.emplace(
              solvers::ReducedLimitedMemoryBfgsSolverT<Backend>(
                reduced,
                compilation.problem->metric(),
                scenario.solver.parameters,
                solvers::LimitedMemoryBfgsDirectionPolicyT<Backend>(
                  scenario.solver.limited_memory_bfgs))
                .solve(initial_control));
            break;
          case ReducedMethod::bfgs:
            throw std::invalid_argument(
              "B1 execution does not implement the full BFGS selection");
        }

      const auto &report_value = *report;
      struct DiscreteDimensions
      {
        std::size_t state_physical;
        std::size_t state_independent;
        std::size_t control_physical;
        std::size_t control_independent;
      };
      std::optional<DiscreteDimensions> dimensions;
      const auto &executable = compilation.problem->executable_model();
      const auto record_model = [this, &report_value, &scenario, &dimensions](
                                  const auto &model) {
        dimensions = DiscreteDimensions{model.physical_state_dimension(),
                                        model.independent_state_dimension(),
                                        model.physical_control_dimension(),
                                        model.independent_control_dimension()};
        if (scenario.experiment.retain_fields &&
            !native_output_directory_.empty())
          model.write_native_output(native_output_directory_,
                                    report_value.final_evaluation.state,
                                    report_value.control,
                                    report_value.final_evaluation.adjoint,
                                    &runtime_->forcing,
                                    &runtime_->desired_state);
      };
      switch (scenario.problem.recipe.discretisation)
        {
          case chapter5::DistributedControlDiscretisation::cellwise_constant:
            {
              const auto *model = dynamic_cast<const
                nmopt::dealii_backend::ScalarDiffusionReactionModel<dim> *>(
                &executable);
              if (model == nullptr)
                throw std::runtime_error(
                  "B1 cellwise execution needs the scalar diffusion model");
              record_model(*model);
              break;
            }
          case chapter5::DistributedControlDiscretisation::
            homogeneous_dirichlet_continuous:
            {
              using ContinuousModel =
                compiler::v1::detail::ContinuousControlModel<dim>;
              const auto *model =
                dynamic_cast<const ContinuousModel *>(&executable);
              if (model == nullptr)
                throw std::runtime_error(
                  "B1 continuous execution needs the continuous-control model");
              record_model(*model);
              break;
            }
        }
      contract::require(dimensions.has_value(),
                        "B1 execution omitted discrete dimension evidence");
      const auto solver_policy =
        experiment::make_reduced_search_policy_snapshot(report_value);

      std::vector<std::string> selected_fields{
        "objective_history",
        "gradient_norm_history",
        "step_length_history",
        "objective_change_history",
        "solve_counts",
        "hessian_evidence"};
      if (scenario.experiment.retain_fields)
        {
          selected_fields.insert(selected_fields.end(),
                                 {"state", "control", "adjoint"});
          selected_fields.insert(selected_fields.end(),
                                 {"target", "forcing", "negative_adjoint"});
        }

      std::vector<benchmark::ArtifactField> fields{
        {"b1.forcing_selection", b1_forcing_selection(scenario)},
        {"b1.control_discretisation",
         chapter5::distributed_control_discretisation_name(
           scenario.problem.recipe.discretisation)},
        {"benchmark.state_dimension",
         std::to_string(report_value.final_evaluation.state.block(0).size())},
        {"benchmark.state_physical_dimension",
         std::to_string(dimensions->state_physical)},
        {"benchmark.state_independent_dimension",
         std::to_string(dimensions->state_independent)},
        {"benchmark.control_dimension",
         std::to_string(report_value.control.block(0).size())},
        {"benchmark.control_physical_dimension",
         std::to_string(dimensions->control_physical)},
        {"benchmark.control_independent_dimension",
         std::to_string(dimensions->control_independent)},
        {"benchmark.adjoint_dimension",
         std::to_string(report_value.final_evaluation.adjoint.block(0).size())},
        {"benchmark.adjoint_physical_dimension",
         std::to_string(dimensions->state_physical)},
        {"benchmark.adjoint_independent_dimension",
         std::to_string(dimensions->state_independent)},
        {"b1.regularisation_weight", b1_number(regularisation_)},
        {"b1.hessian_evidence", "centered_finite_difference"},
        {"b1.hessian_direction_norm",
         b1_number(hessian_evidence.direction_norm)},
        {"b1.hessian_symmetry_error",
         b1_number(hessian_evidence.symmetry_error)},
        {"b1.hessian_finite_difference_step",
         b1_number(hessian_evidence.finite_difference_step)},
        {"b1.hessian_finite_difference_error",
         b1_number(hessian_evidence.finite_difference_error)},
        {"b1.hessian_finite_difference_passed", "true"},
        {"solver.method", b1_method_name(scenario.solver.method)},
        {"solver.initial_control", scenario.solver.initial_control},
        {"solver.policy", report_value.policy_name},
        {"solver.stopping_reason",
         solvers::reduced_stopping_reason_name(report_value.stopping_reason)},
        {"solver.accepted_iterations",
         std::to_string(report_value.accepted_iterations)},
        {"solver.line_search_trial_count",
         std::to_string(report_value.line_search_trial_count)},
        {"solver.state_solve_count",
         std::to_string(report_value.state_solve_count)},
        {"solver.adjoint_solve_count",
         std::to_string(report_value.adjoint_solve_count)},
        {"solver.metric_solve_count",
         std::to_string(report_value.metric_solve_count)},
        {"solver.hessian_action_count",
         std::to_string(report_value.hessian_action_count)},
        {"solver.direction_reset_count",
         std::to_string(report_value.direction_reset_count)},
        {"solver.objective_history",
         b1_history(report_value.objective_history)},
        {"solver.gradient_norm_history",
         b1_history(report_value.gradient_norm_history)},
        {"solver.step_length_history",
         b1_history(report_value.step_length_history)},
        {"solver.step_norm_history",
         b1_history(report_value.step_norm_history)},
        {"solver.objective_change_history",
         b1_history(report_value.objective_change_history)}};
      if (scenario.solver.declared_minimum_step_length.has_value())
        fields.push_back({"solver.declared_minimum_step_length",
                          b1_number(
                            *scenario.solver.declared_minimum_step_length)});
      if (scenario.solver.method == ReducedMethod::limited_memory_bfgs)
        {
          fields.push_back(
            {"solver.l_bfgs_memory",
             std::to_string(
               scenario.solver.limited_memory_bfgs.memory_size)});
          fields.push_back(
            {"solver.l_bfgs_curvature_tolerance",
             b1_number(
               scenario.solver.limited_memory_bfgs.curvature_tolerance)});
          fields.push_back(
            {"solver.l_bfgs_initial_scaling",
             solvers::limited_memory_bfgs_initial_scaling_name(
               scenario.solver.limited_memory_bfgs.
                 initial_inverse_hessian_scaling)});
        }
      Envelope envelope{compilation.problem->manifest(),
                        solver_policy,
              std::move(*report),
              environment_};
      return {std::move(envelope),
              std::move(compilation.diagnostics),
              {},
              std::move(selected_fields),
              std::move(fields)};
    }

  private:
    static const char *
    b1_forcing_selection(const B1Scenario &scenario)
    {
      return chapter6::b1_forcing_selection_name(
        scenario.problem.forcing_selection);
    }

    double                                                    regularisation_;
    const B1RuntimeDataT<dim> *                               runtime_;
    std::shared_ptr<compiler::v1::DealiiCompilationSession<dim>> session_;
    experiment::RunEnvironmentRecord                             environment_;
    std::filesystem::path                                         native_output_directory_;
  };
} // namespace nmopt::application::chapter6::dealii
