#include "problem_a_binding.hpp"

#include "nmopt/solvers/reduced_gradient.hpp"

#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

namespace
{
  using Backend = nmopt::dealii_backend::SerialBackend;
  using Binding = external_dealii_step4::minimal::ProblemABinding;
  using Parameters = nmopt::solvers::ReducedSolverParameters;
  using ProblemA = external_dealii_step4::ProblemA;
  using Solver = nmopt::solvers::ReducedSearchSolverT<Backend>;
  using ProblemVector = ProblemA::Vector;

  Parameters
  frozen_problem_a_parameters()
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
                    std::filesystem::path("solution.vtk");
      if (output.has_parent_path())
        std::filesystem::create_directories(output.parent_path());

      ProblemA problem;
      Binding binding(problem);
      ProblemVector initial_control(problem.control_dimension());
      initial_control = 0.0;

      const auto initial =
        nmopt::contract::PrimalBlockT<Backend>(binding.control_layout(),
                                               std::vector<Backend::Vector>{
                                                 initial_control});
      Solver solver(binding.reduced(),
                    binding.metric(),
                    frozen_problem_a_parameters());
      const auto result = solver.solve(initial);
      if (result.stopping_reason !=
          nmopt::solvers::ReducedStoppingReason::gradient_tolerance)
        throw std::runtime_error(
          "Problem A optimization did not stop by gradient tolerance");

      const auto &state = result.final_evaluation.state.block(0);
      binding.problem().output_results(state, output);

      std::cout << "minimal Problem A consumer completed\n"
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
      std::cerr << "minimal Problem A consumer failed: " << exception.what()
                << '\n';
      return 1;
    }
}
