#include "problem_b_binding.hpp"

#include "nmopt/solvers/reduced_gradient.hpp"

#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  using Backend   = nmopt::dealii_backend::SerialBackend;
  using Binding   = external_dealii_step4::minimal::ProblemBBinding;
  using Parameters = nmopt::solvers::ReducedSolverParameters;
  using Problem   = Binding::Problem;
  using Solver    = nmopt::solvers::ReducedSearchSolverT<Backend>;
  using ProblemVector = Problem::Vector;

  Parameters
  frozen_problem_b_parameters()
  {
    Parameters parameters;
    parameters.maximum_iterations = 5000;
    parameters.maximum_line_search_trials = 30;
    parameters.gradient_tolerance = 1.0e-6;
    parameters.stopping_criterion =
      nmopt::solvers::ReducedStoppingCriterion::gradient_norm;
    parameters.relative_gradient_tolerance = 0.0;
    parameters.objective_change_tolerance = 0.0;
    parameters.step_tolerance = 0.0;
    parameters.objective_target = std::nullopt;
    parameters.initial_step_length = 1.0;
    parameters.minimum_step_length = 0.0;
    parameters.armijo_fraction = 1.0e-4;
    parameters.backtracking_factor = 0.5;
    return parameters;
  }

  std::filesystem::path
  default_output_path()
  {
    const auto root = std::filesystem::current_path() /
                      "runs/external-dealii/step-4/minimal/manual/problem-b";
    std::filesystem::create_directories(root);
    const auto timestamp = std::chrono::duration_cast<std::chrono::microseconds>(
      std::chrono::steady_clock::now().time_since_epoch()).count();
    for (unsigned int attempt = 0;; ++attempt)
      {
        const auto directory = root /
                               (std::to_string(timestamp) + "-" +
                                std::to_string(attempt));
        if (std::filesystem::create_directory(directory))
          return directory / "solution.vtk";
      }
  }
}

int
main(const int argc, char **argv)
{
  if (argc > 2)
    {
      std::cerr << "usage: " << argv[0] << " [output.vtk]\n";
      return 2;
    }

  try
    {
      const std::filesystem::path output =
        argc == 2 ? std::filesystem::path(argv[1]) :
                    default_output_path();
      if (output.has_parent_path())
        std::filesystem::create_directories(output.parent_path());

      Binding::Application application;
      application.prepare_for_external_use();
      Problem problem(application);
      // nmopt integration: bind the native problem, wrap the initial control,
      // and solve the reduced problem; physical-state reconstruction remains
      // application-owned below.
      Binding binding(problem);

      ProblemVector initial_control(problem.control_dimension());
      initial_control = 0.0;
      const auto initial =
        nmopt::contract::PrimalBlockT<Backend>(
          binding.control_layout(), std::vector<Backend::Vector>{
                                      initial_control});

      Solver solver(binding.reduced(),
                    binding.metric(),
                    frozen_problem_b_parameters());
      const auto result = solver.solve(initial);
      if (result.stopping_reason !=
          nmopt::solvers::ReducedStoppingReason::gradient_tolerance)
        throw std::runtime_error(
          "Problem B optimization did not stop by gradient tolerance");

      const auto full_state = problem.coordinates().reconstruct(
        result.final_evaluation.state.block(0));
      application.output_results(full_state, output);

      std::cout << std::setprecision(17)
                << "minimal Problem B consumer completed\n"
                << "stopping_reason gradient_tolerance\n"
                << "free_state_dimension " << problem.state_dimension() << '\n'
                << "control_dimension " << problem.control_dimension() << '\n'
                << "accepted_iterations " << result.accepted_iterations << '\n'
                << "line_search_trials " << result.line_search_trial_count << '\n'
                << "final_objective "
                << result.final_evaluation.objective_value << '\n'
                << "final_gradient_norm "
                << result.gradient_norm_history.back() << '\n'
                << "output " << output.string() << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "minimal Problem B consumer failed: " << exception.what()
                << '\n';
      return 1;
    }
}
