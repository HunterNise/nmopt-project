#include "../../apps/external-dealii/step-4/verification/verification.hpp"
#include "../../apps/external-dealii/step-4/evaluation/native_optimization.hpp"

#include "external_step4_evidence.hpp"
#include "../support/scenario_dispatch.hpp"

#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
  using Matrix          = external_dealii_step4::ProblemA::Matrix;
  using ProblemA        = external_dealii_step4::ProblemA;
  using NativeReduced   = external_dealii_step4::NativeReduced;
  using NativeArmijoSolver = external_dealii_step4::NativeArmijoSolver;
  using NativeOptimizationResult =
    external_dealii_step4::NativeOptimizationResult;
  using NativeOptimizationStoppingReason =
    external_dealii_step4::NativeOptimizationStoppingReason;
  using OptimizationPolicy = external_dealii_step4::OptimizationPolicy;
  using Instrumentation = external_dealii_step4::Instrumentation;
  using Vector          = ProblemA::Vector;
  namespace verification = external_dealii_step4::verification;

  void
  require(const bool condition, const std::string &message)
  {
    if (!condition)
      {
        external_dealii_step4_test::note_failure(message);
        throw std::runtime_error(message);
      }
  }

  void
  require_close(const double       actual,
                const double       expected,
                const double       tolerance,
                const std::string &message)
  {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        !std::isfinite(tolerance) || std::abs(actual - expected) > tolerance)
      {
        std::ostringstream detail;
        detail << message << ": actual=" << actual
               << ", expected=" << expected
               << ", tolerance=" << tolerance;
        external_dealii_step4_test::note_failure(detail.str());
        throw std::runtime_error(detail.str());
      }
  }

  template <typename Callable>
  void
  require_rejected(Callable &&callable, const std::string &message)
  {
    bool rejected = false;
    try
      {
        callable();
      }
    catch (const std::exception &)
      {
        rejected = true;
      }
    require(rejected, message);
  }

  double
  matrix_difference(const std::vector<double> &left, const Matrix &right)
  {
    require(left.size() ==
              static_cast<std::size_t>(right.m()) * right.n(),
            "Step4 matrix dimensions changed");

    double squared_difference = 0.0;
    for (unsigned int row = 0; row < right.m(); ++row)
      for (unsigned int column = 0; column < right.n(); ++column)
        {
          const double difference =
            left[static_cast<std::size_t>(row) * right.n() + column] -
            right.el(row, column);
          squared_difference += difference * difference;
        }
    return std::sqrt(squared_difference);
  }

  std::vector<double>
  matrix_values(const Matrix &matrix)
  {
    std::vector<double> values(static_cast<std::size_t>(matrix.m()) *
                               matrix.n());
    for (unsigned int row = 0; row < matrix.m(); ++row)
      for (unsigned int column = 0; column < matrix.n(); ++column)
        values[static_cast<std::size_t>(row) * matrix.n() + column] =
          matrix.el(row, column);
    return values;
  }

  double
  matrix_symmetry_error(const Matrix &matrix)
  {
    double squared_difference = 0.0;
    double squared_norm       = 0.0;
    for (unsigned int row = 0; row < matrix.m(); ++row)
      for (unsigned int column = 0; column < matrix.n(); ++column)
        {
          const double value = matrix.el(row, column);
          squared_norm += value * value;
          const double difference = value - matrix.el(column, row);
          squared_difference += difference * difference;
        }
    return std::sqrt(squared_difference) /
           std::max(1.0, std::sqrt(squared_norm));
  }

  double
  vector_difference(const Vector &left, const Vector &right)
  {
    require(left.size() == right.size(), "Step4 vector dimensions changed");
    Vector difference = left;
    difference.add(-1.0, right);
    return difference.l2_norm();
  }

  std::filesystem::path
  native_reference_artifact_root(const char *const scenario)
  {
    auto directory = std::filesystem::current_path();
    while (true)
      {
        if (std::filesystem::exists(
              directory / "apps/external-dealii/step-4/source/upstream/step-4.cc"))
          {
            return external_dealii_step4_test::create_unique_artifact_root(
              directory / "runs/external-dealii/step-4/working/"
                        "native-reference",
              scenario);
          }

        const auto parent = directory.parent_path();
        if (parent == directory)
          break;
        directory = parent;
      }

    throw std::runtime_error("could not locate the repository artifact root");
  }

  void
  require_original_assembly(const std::vector<double> &matrix_before,
                            const Vector &             rhs_before,
                            const Step4<2> &           tutorial)
  {
    require_close(matrix_difference(matrix_before,
                                    tutorial.system_matrix_view()),
                  0.0,
                  0.0,
                  "Step4 assembled matrix changed during solve");
    require_close(vector_difference(rhs_before, tutorial.system_rhs_view()),
                  0.0,
                  0.0,
                  "Step4 assembled RHS changed during solve");
  }

  void
  run_matrix_rhs_solve_contract()
  {
    Step4<2> tutorial;
    tutorial.prepare_for_external_use();

    const auto &matrix = tutorial.system_matrix_view();
    const auto &rhs    = tutorial.system_rhs_view();
    require(matrix.m() == 289 && matrix.n() == 289,
            "Step4 2D assembled matrix has the wrong dimension");
    require(rhs.size() == 289, "Step4 2D assembled RHS has the wrong dimension");
    require(matrix_symmetry_error(matrix) <= 1e-14,
            "Step4 2D assembled matrix is not symmetric");

    const auto matrix_before = matrix_values(matrix);
    const Vector rhs_before  = rhs;

    Vector state(rhs.size());
    state = 0.0;
    const auto evidence = tutorial.solve(rhs, state);
    require(evidence.converged, "Step4 supplied-RHS solve did not converge");
    require(evidence.iterations > 0,
            "Step4 supplied-RHS solve forced no work for a nonzero RHS");
    require(evidence.final_residual <= 1e-10,
            "Step4 supplied-RHS solve has a large residual");
    require_original_assembly(matrix_before, rhs_before, tutorial);

    Vector second_state(rhs.size());
    second_state = 0.5;
    const auto second_evidence = tutorial.solve(rhs, second_state);
    require(second_evidence.converged,
            "Step4 nonzero-initial supplied-RHS solve did not converge");
    require(second_evidence.final_residual <= 1e-10,
            "Step4 nonzero-initial solve has a large residual");
    require_original_assembly(matrix_before, rhs_before, tutorial);

    Vector zero_rhs(rhs.size());
    zero_rhs = 0.0;
    Vector zero_state(rhs.size());
    zero_state = 0.0;
    const auto zero_evidence = tutorial.solve(zero_rhs, zero_state);
    require(zero_evidence.converged,
            "Step4 zero-RHS zero-initial solve did not converge");
    require(zero_evidence.iterations == 0,
            "Step4 zero-RHS zero-initial solve forced iterations");
    require(zero_evidence.final_residual == 0.0,
            "Step4 zero-RHS zero-initial solve changed the solution");
    require_original_assembly(matrix_before, rhs_before, tutorial);
  }

  void
  run_fe_access_contract()
  {
    Step4<2> tutorial;
    tutorial.prepare_for_external_use();

    const auto &dof_handler     = tutorial.dof_handler_view();
    const auto &finite_element  = dof_handler.get_fe();
    const auto &boundary_values = tutorial.boundary_values_view();
    require(dof_handler.n_dofs() == 289,
            "Step4 FE access view has the wrong DoF dimension");
    require(finite_element.degree == 1 &&
              finite_element.n_dofs_per_cell() == 4,
            "Step4 FE access view has the wrong Q1 element");
    require(boundary_values.size() == 64,
            "Step4 boundary access view has the wrong boundary count");
    for (const auto &[index, value] : boundary_values)
      {
        require(index < dof_handler.n_dofs(),
                "Step4 boundary access view contains an invalid DoF index");
        require(std::isfinite(value),
                "Step4 boundary access view contains a non-finite value");
      }
  }

  void
  run_supplied_state_output_contract()
  {
    Step4<2> tutorial;
    tutorial.prepare_for_external_use();

    Vector state(tutorial.system_rhs_view().size());
    state = 7.0;

    const auto output_directory =
      std::filesystem::temp_directory_path() / "nmopt-external-step4-output";
    std::filesystem::remove_all(output_directory);
    std::filesystem::create_directories(output_directory);
    const auto output_file = output_directory / "supplied-state.vtk";
    tutorial.output_results(state, output_file);

    std::ifstream input(output_file);
    require(static_cast<bool>(input),
            "Step4 supplied-state output could not be opened");
    const std::string contents((std::istreambuf_iterator<char>(input)), {});
    const auto lookup_table = contents.find("LOOKUP_TABLE default");
    require(contents.find("# vtk DataFile Version") != std::string::npos,
            "Step4 supplied-state output is not legacy VTK");
    require(contents.find("SCALARS solution double 1") != std::string::npos,
            "Step4 supplied-state output changed the field identity");
    require(lookup_table != std::string::npos &&
              contents.find("7", lookup_table) != std::string::npos,
            "Step4 supplied-state output did not use the supplied vector");

    std::filesystem::remove_all(output_directory);
  }

  void
  run_problem_a_operations_contract()
  {
    Instrumentation instrumentation;
    ProblemA        problem(instrumentation);
    require(problem.state_dimension() ==
              external_dealii_step4::scenario::dimension,
            "Problem A state dimension does not match the frozen scenario");
    require(problem.control_dimension() ==
              external_dealii_step4::scenario::dimension,
            "Problem A control dimension does not match the frozen scenario");
    require(instrumentation.assembly_calls == 1,
            "Problem A did not assemble exactly once during construction");

    const auto controls = external_dealii_step4::scenario::reduced_controls();
    require(controls.size() == 4, "Problem A scenario controls are incomplete");
    for (const auto &control : controls)
      require(control.size() == problem.control_dimension(),
              "Problem A scenario control has the wrong dimension");

    const auto state_result = problem.solve_state(controls.front());
    require(state_result.evidence.converged,
            "Problem A zero-control state solve did not converge");
    require(state_result.evidence.final_residual <= 1e-10,
            "Problem A zero-control state solve has a large residual");

    const auto objective =
      problem.objective(state_result.solution, controls.front());
    require(objective >= 0.0, "Problem A objective is negative");
    const auto objective_derivative =
      problem.objective_derivative(state_result.solution, controls.front());
    require(vector_difference(objective_derivative.state, state_result.solution) ==
              0.0,
            "Problem A state objective derivative has the wrong value");
    require(vector_difference(objective_derivative.control, controls.front()) ==
              0.0,
            "Problem A control objective derivative has the wrong value");

    const auto adjoint_result =
      problem.solve_adjoint(objective_derivative.state);
    require(adjoint_result.evidence.converged,
            "Problem A adjoint solve did not converge");
    require(adjoint_result.evidence.final_residual <= 1e-10,
            "Problem A adjoint solve has a large residual");

    Vector seed(problem.state_dimension());
    seed = 0.25;
    const auto control_pullback = problem.control_vjp(seed);
    seed *= -1.0;
    require(vector_difference(control_pullback, seed) == 0.0,
            "Problem A control pullback has the wrong sign");
  }

  void
  run_native_reduced_contract()
  {
    Instrumentation instrumentation;
    ProblemA        problem(instrumentation);
    NativeReduced  reduced(problem, instrumentation);
    const auto     controls = external_dealii_step4::scenario::reduced_controls();

    const auto value = reduced.evaluate_value(controls.front());
    require(instrumentation.value_evaluations == 1,
            "native reduced value stage was not counted");
    require(instrumentation.derivative_augmentations == 0,
            "native reduced derivative ran during value evaluation");
    require(instrumentation.state_solve_calls == 1,
            "native reduced value stage did not perform one state solve");
    require(instrumentation.adjoint_solve_calls == 0,
            "native reduced value stage performed an adjoint solve");
    require(instrumentation.objective_calls == 1,
            "native reduced value stage did not evaluate the objective once");

    const auto derivative = reduced.augment_derivative(value);
    require(instrumentation.derivative_augmentations == 1,
            "native reduced derivative stage was not counted");
    require(instrumentation.state_solve_calls == 1,
            "native reduced derivative stage repeated the state solve");
    require(instrumentation.adjoint_solve_calls == 1,
            "native reduced derivative stage did not perform one adjoint solve");
    require(instrumentation.objective_derivative_calls == 1,
            "native reduced derivative stage did not evaluate objective partials once");
    require(instrumentation.control_vjp_calls == 1,
            "native reduced derivative stage did not perform one control pullback");

    Vector expected_gradient = derivative.control_derivative;
    expected_gradient += derivative.adjoint;
    require(vector_difference(derivative.reduced_derivative, expected_gradient) ==
              0.0,
            "native reduced derivative has the wrong control sign");
    require(vector_difference(value.state, derivative.state_derivative) == 0.0,
            "native reduced derivative changed the retained state");

    const auto repeated_before = reduced.evaluate_value(controls[2]);
    const auto intervening      = reduced.evaluate_value(controls[3]);
    const auto repeated_after   = reduced.evaluate_value(controls[2]);
    (void)intervening;
    require(vector_difference(repeated_before.state, repeated_after.state) == 0.0,
            "native reduced repeated control changed the state");
    require_close(repeated_before.objective,
                  repeated_after.objective,
                  0.0,
                  "native reduced repeated control changed the objective");

    const auto output_directory =
      std::filesystem::temp_directory_path() /
      "nmopt-external-step4-native-reduced-output";
    std::filesystem::remove_all(output_directory);
    std::filesystem::create_directories(output_directory);
    const auto output_file = output_directory / "retained-state.vtk";
    problem.output_results(repeated_after.state, output_file);
    require(std::filesystem::exists(output_file),
            "native reduced retained-state output is missing");
    require(instrumentation.output_calls == 1,
            "native reduced output was not counted");
    std::filesystem::remove_all(output_directory);
  }

  void
  run_native_derivative_verification()
  {
    Instrumentation instrumentation;
    const auto artifact_root = native_reference_artifact_root("derivatives");
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact_root, "native_derivatives", {{"native", &instrumentation}});
    try
      {
    std::ofstream finite_difference_output(
      artifact_root / "reduced-finite-differences.csv");
    require(static_cast<bool>(finite_difference_output),
            "could not open the reduced finite-difference trace");
    finite_difference_output
      << "control,direction,step,analytic,centered,error,bound\n";

    std::ofstream taylor_output(artifact_root / "reduced-taylor.csv");
    require(static_cast<bool>(taylor_output),
            "could not open the reduced Taylor trace");
    taylor_output << "direction,step,remainder,ratio,repeat_variation\n";

    std::ofstream summary(
      artifact_root / "native-derivative-verification.txt");
    require(static_cast<bool>(summary),
            "could not open the native derivative verification summary");
    summary << "status running\n";

    ProblemA      problem(instrumentation);
    NativeReduced reduced(problem, instrumentation);
    const auto    point = verification::make_off_solution_point();

    const auto centered_jvp =
      verification::centered_residual_jvp(problem, point, 1.0e-6);
    const auto analytic_jvp =
      problem.residual_jvp(point.state_tangent, point.control_tangent);
    const double residual_jvp_error =
      verification::scaled_vector_error(centered_jvp, analytic_jvp);
    summary << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "residual_jvp_scaled_error " << residual_jvp_error << '\n';
    require(residual_jvp_error <= 1.0e-8,
            "residual JVP centered difference does not agree");

    const auto full_vjp = problem.residual_vjp(point.test_seed);
    const auto zero = [&point]() {
      Vector value(point.state.size());
      value = 0.0;
      return value;
    }();
    const auto state_only_jvp =
      problem.residual_jvp(point.state_tangent, zero);
    const auto control_only_jvp =
      problem.residual_jvp(zero, point.control_tangent);
    const double full_left = analytic_jvp * point.test_seed;
    const double full_right =
      full_vjp.state * point.state_tangent +
      full_vjp.control * point.control_tangent;
    const double state_left = state_only_jvp * point.test_seed;
    const double state_right = full_vjp.state * point.state_tangent;
    const double control_left = control_only_jvp * point.test_seed;
    const double control_right = full_vjp.control * point.control_tangent;
    const double full_pairing_error =
      verification::scaled_scalar_error(full_left, full_right);
    const double state_pairing_error =
      verification::scaled_scalar_error(state_left, state_right);
    const double control_pairing_error =
      verification::scaled_scalar_error(control_left, control_right);
    summary << "full_pairing_scaled_error " << full_pairing_error << '\n'
            << "state_pairing_scaled_error " << state_pairing_error << '\n'
            << "control_pairing_scaled_error " << control_pairing_error
            << '\n';
    require(full_pairing_error <= 1.0e-12,
            "full residual JVP/VJP pairing failed");
    require(state_pairing_error <= 1.0e-12,
            "state-only residual JVP/VJP pairing failed");
    require(control_pairing_error <= 1.0e-12,
            "control-only residual JVP/VJP pairing failed");

    const auto objective_derivative =
      problem.objective_derivative(point.state, point.control);
    const double analytic_objective_derivative =
      objective_derivative.state * point.state_tangent +
      objective_derivative.control * point.control_tangent;
    const double centered_objective_derivative =
      verification::centered_objective_derivative(problem, point, 1.0e-6);
    const double objective_derivative_error = verification::scaled_scalar_error(
      analytic_objective_derivative, centered_objective_derivative);
    summary << "objective_derivative_scaled_error "
            << objective_derivative_error << '\n';
    require(objective_derivative_error <= 1.0e-8,
            "objective directional derivative does not agree");

    const auto controls = external_dealii_step4::scenario::reduced_controls();
    const std::vector<std::pair<const char *, Vector>> directions{
      {"r", external_dealii_step4::scenario::normalized_ramp_direction()},
      {"a",
       external_dealii_step4::scenario::normalized_alternating_direction()}};
    const std::vector<double> finite_difference_steps{
      1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6};

    for (std::size_t control_index = 0; control_index < controls.size();
         ++control_index)
      {
        const auto value = reduced.evaluate_value(controls[control_index]);
        const auto derivative = reduced.augment_derivative(value);

        for (const auto &[direction_name, direction] : directions)
          {
            const double analytic = derivative.reduced_derivative * direction;
            std::vector<bool> passes;
            passes.reserve(finite_difference_steps.size());

            for (const double step : finite_difference_steps)
              {
                Vector control_plus  = controls[control_index];
                Vector control_minus = controls[control_index];
                control_plus.add(step, direction);
                control_minus.add(-step, direction);
                const auto plus  = reduced.evaluate_value(control_plus);
                const auto minus = reduced.evaluate_value(control_minus);
                const double centered =
                  (plus.objective - minus.objective) / (2.0 * step);
                const double error = std::abs(centered - analytic);
                const double bound =
                  1.0e-7 * std::max(1.0, std::abs(analytic));
                passes.push_back(error <= bound);
                finite_difference_output
                  << control_index << ',' << direction_name << ','
                  << std::setprecision(
                       std::numeric_limits<double>::max_digits10)
                  << step << ',' << analytic << ',' << centered << ','
                  << error << ',' << bound << '\n';
              }

            bool adjacent_passes = false;
            for (std::size_t index = 1; index < passes.size(); ++index)
              adjacent_passes = adjacent_passes ||
                                (passes[index - 1] && passes[index]);
            require(adjacent_passes,
                    "reduced centered finite differences lack adjacent usable steps");
          }
      }

    const Vector zero_control = controls.front();
    const auto base_value = reduced.evaluate_value(zero_control);
    const auto base_derivative = reduced.augment_derivative(base_value);
    const auto repeated_value = reduced.evaluate_value(zero_control);
    const double repeat_variation =
      std::abs(repeated_value.objective - base_value.objective);
    const std::vector<double> taylor_steps{0.1, 0.05, 0.025};

    for (const auto &[direction_name, direction] : directions)
      {
        const double slope = base_derivative.reduced_derivative * direction;
        std::vector<double> remainders;
        remainders.reserve(taylor_steps.size());
        for (std::size_t index = 0; index < taylor_steps.size(); ++index)
          {
            const double step = taylor_steps[index];
            Vector control = zero_control;
            control.add(step, direction);
            const auto trial = reduced.evaluate_value(control);
            const double remainder =
              trial.objective - base_value.objective - step * slope;
            const double ratio = index == 0 ?
                                   std::numeric_limits<double>::quiet_NaN() :
                                   remainders[index - 1] / remainder;
            remainders.push_back(remainder);
            taylor_output
              << direction_name << ','
              << std::setprecision(std::numeric_limits<double>::max_digits10)
              << step << ',' << remainder << ',' << ratio << ','
              << repeat_variation << '\n';
            require(std::isfinite(remainder) &&
                      remainder > 10.0 * repeat_variation,
                    "reduced Taylor remainder is not positive and resolved");
          }

        for (std::size_t index = 1; index < remainders.size(); ++index)
          {
            const double ratio = remainders[index - 1] / remainders[index];
            require(ratio >= 3.5 && ratio <= 4.5,
                    "reduced Taylor remainder does not have quadratic halving");
          }
      }

    summary << "value_evaluations " << instrumentation.value_evaluations << '\n'
            << "derivative_augmentations "
            << instrumentation.derivative_augmentations << '\n'
            << "state_solve_calls " << instrumentation.state_solve_calls << '\n'
            << "adjoint_solve_calls " << instrumentation.adjoint_solve_calls
            << '\n';
    summary.flush();
    evidence.complete();
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
  }

  void
  run_native_oracle_contract()
  {
    Instrumentation instrumentation;
    const auto artifact_root = native_reference_artifact_root("oracle");
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact_root, "native_oracle", {{"native", &instrumentation}});
    try
      {
    std::ofstream output(artifact_root / "native-oracle.txt");
    require(static_cast<bool>(output),
            "could not open the native oracle trace");
    output << "status running\n";

    ProblemA problem(instrumentation);
    const double symmetry_error =
      verification::matrix_symmetry_error(problem.system_matrix());
    output << std::setprecision(std::numeric_limits<double>::max_digits10)
           << "matrix_symmetry_error " << symmetry_error << '\n';
    require(symmetry_error <= 1.0e-14,
            "Problem A assembled matrix failed the symmetry check");

    const auto oracle = verification::optimum_oracle(problem);
    output << "system_residual " << oracle.system_residual << '\n'
           << "stationarity_residual " << oracle.stationarity_residual << '\n';
    require(oracle.system_residual <= 1.0e-10,
            "dense optimum oracle has a large system residual");
    require(oracle.stationarity_residual <= 1.0e-10,
            "dense optimum oracle has a large stationarity residual");

    NativeReduced reduced(problem, instrumentation);
    const auto value = reduced.evaluate_value(oracle.control);
    const auto derivative = reduced.augment_derivative(value);
    const double state_error =
      verification::scaled_vector_error(value.state, oracle.state);
    const double gradient_norm = derivative.reduced_derivative.l2_norm();
    output << "native_state_scaled_error " << state_error << '\n'
           << "native_oracle_gradient_norm " << gradient_norm << '\n'
           << "native_state_solve_iterations " << value.state_solve.iterations
           << '\n'
           << "native_adjoint_solve_iterations "
           << derivative.adjoint_solve.iterations << '\n';
    verification::require_vector_close(value.state,
                                       oracle.state,
                                       1.0e-11,
                                       1.0e-10,
                                       "native state differs from dense oracle");
    verification::require_vector_close(
      value.control,
      oracle.control,
      0.0,
      0.0,
      "oracle control changed during evaluation");
    require(gradient_norm <= 1.0e-8,
            "native reduced gradient is not zero at the dense oracle");
    output.flush();
    evidence.complete();
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
  }

  std::filesystem::path
  native_optimization_artifact_root()
  {
    auto directory = std::filesystem::current_path();
    while (true)
      {
        if (std::filesystem::exists(
              directory / "apps/external-dealii/step-4/source/upstream/step-4.cc"))
          {
            const auto root =
              external_dealii_step4_test::create_unique_artifact_root(
                directory / "runs/external-dealii/step-4/optimization",
                "native");
            const auto native_root = root / "native";
            std::filesystem::create_directories(native_root);
            return native_root;
          }

        const auto parent = directory.parent_path();
        if (parent == directory)
          break;
        directory = parent;
      }

    throw std::runtime_error("could not locate the optimization artifact root");
  }

  void
  write_native_optimization_trace(const std::filesystem::path &root,
                                  const NativeOptimizationResult &result)
  {
    std::ofstream trace(root / "trace.csv");
    require(static_cast<bool>(trace), "could not open the native optimization trace");
    trace << "record,iteration,trial,step_length,objective,actual_slope,"
             "armijo_bound,objective_finite,slope_negative,accepted,"
             "objective_before,objective_after,objective_change,"
             "actual_step_norm,gradient_norm\n";
    trace << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const auto &trial : result.trial_records)
      trace << "trial," << trial.iteration << ',' << trial.trial << ','
            << trial.step_length << ',' << trial.objective_value << ','
            << trial.actual_slope << ',' << trial.sufficient_decrease_bound
            << ',' << trial.objective_finite << ',' << trial.slope_negative
            << ',' << trial.accepted << ",,,,,\n";
    for (const auto &iteration : result.accepted_iterations)
      trace << "accepted," << iteration.iteration << ",,"
            << iteration.requested_step_length << ','
            << iteration.objective_after << ',' << iteration.actual_slope
            << ",,,,1," << iteration.objective_before << ','
            << iteration.objective_after << ',' << iteration.objective_change
            << ',' << iteration.actual_step_norm << ','
            << iteration.gradient_norm << '\n';

  }

  void
  write_native_optimization_summary(const std::filesystem::path &root,
                                    const NativeOptimizationResult &result,
                                    const Instrumentation &instrumentation,
                                    const verification::OracleResult &oracle)
  {
    std::ofstream summary(root / "summary.txt");
    require(static_cast<bool>(summary),
            "could not open the native optimization summary");
    summary << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "stopping_reason "
            << external_dealii_step4::native_optimization_stopping_reason_name(
                 result.stopping_reason)
            << '\n'
            << "accepted_iterations " << result.accepted_iteration_count << '\n'
            << "line_search_trials " << result.line_search_trial_count << '\n'
            << "final_objective " << result.value.objective << '\n'
            << "final_gradient_norm "
            << result.derivative.reduced_derivative.l2_norm() << '\n'
            << "oracle_system_residual " << oracle.system_residual << '\n'
            << "oracle_stationarity_residual " << oracle.stationarity_residual
            << '\n'
            << "oracle_control_distance_bound " << 2.0e-6 << '\n'
            << "final_control_oracle_error "
            << vector_difference(result.value.control, oracle.control) << '\n'
            << "oracle_objective "
            << (0.5 * (oracle.state * oracle.state) +
                0.5 * (oracle.control * oracle.control)) << '\n'
            << "final_objective_gap "
            << (result.value.objective -
                (0.5 * (oracle.state * oracle.state) +
                 0.5 * (oracle.control * oracle.control))) << '\n'
            << "assembly_calls " << instrumentation.assembly_calls << '\n'
            << "state_solve_calls " << instrumentation.state_solve_calls << '\n'
            << "adjoint_solve_calls " << instrumentation.adjoint_solve_calls
            << '\n'
            << "value_evaluations " << instrumentation.value_evaluations << '\n'
            << "derivative_augmentations "
            << instrumentation.derivative_augmentations << '\n';
  }

  void
  run_native_optimization_contract(
    const bool fail_after_trace = false,
    std::filesystem::path *const failure_artifact = nullptr)
  {
    Instrumentation instrumentation;
    Instrumentation verification_instrumentation;
    const auto artifact_root = native_optimization_artifact_root();
    if (failure_artifact != nullptr)
      *failure_artifact = artifact_root;
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact_root,
      "native_optimization",
      {{"native", &instrumentation},
       {"verification", &verification_instrumentation}});
    try
      {
    ProblemA      problem(instrumentation);
    NativeReduced reduced(problem, instrumentation);
    Vector        initial_control(problem.control_dimension());
    initial_control = 0.0;

    const OptimizationPolicy policy =
      external_dealii_step4::frozen_optimization_policy();
    NativeArmijoSolver solver(reduced, policy);
    const auto result = solver.solve(initial_control);

    write_native_optimization_trace(artifact_root, result);
    if (fail_after_trace)
      throw std::runtime_error("failure after completed native optimization trace");

    verification::require_finite(result.value.control,
                                  "native optimization final control");
    verification::require_finite(result.value.state,
                                  "native optimization final state");
    verification::require_finite(result.value.objective,
                                  "native optimization final objective");
    verification::require_finite(result.derivative.reduced_derivative,
                                  "native optimization final gradient");

    require(result.stopping_reason ==
              NativeOptimizationStoppingReason::gradient_tolerance,
            "native optimization did not stop by gradient tolerance");
    require(result.accepted_iteration_count > 0,
            "native optimization accepted no iterations");
    const double final_gradient_norm =
      result.derivative.reduced_derivative.l2_norm();
    require(std::isfinite(final_gradient_norm) && final_gradient_norm <= 1.1e-6,
            "native optimization final gradient exceeds the audit bound");
    require(result.objective_history.size() ==
              result.accepted_iteration_count + 1,
            "native optimization objective history has the wrong size");
    require(result.accepted_iterations.size() ==
              result.accepted_iteration_count,
            "native optimization accepted trace has the wrong size");
    require(instrumentation.assembly_calls == 1,
            "native optimization assembled more than once");
    require(instrumentation.state_solve_calls ==
              1 + result.line_search_trial_count,
            "native optimization state count violates the staged schedule");
    require(instrumentation.adjoint_solve_calls ==
              1 + result.accepted_iteration_count,
            "native optimization adjoint count violates state reuse");
    require(instrumentation.value_evaluations ==
              1 + result.line_search_trial_count,
            "native optimization value count violates the trial schedule");
    require(instrumentation.derivative_augmentations ==
              1 + result.accepted_iteration_count,
            "native optimization derivative count includes rejected trials");
    require(instrumentation.objective_calls ==
              1 + result.line_search_trial_count &&
              instrumentation.objective_derivative_calls ==
                1 + result.accepted_iteration_count,
            "native optimization objective callback counts are inconsistent");
    require(instrumentation.residual_calls == 0 &&
              instrumentation.residual_jvp_calls == 0 &&
              instrumentation.residual_vjp_calls == 0,
            "native optimization unexpectedly used residual callbacks");
    require(instrumentation.control_vjp_calls ==
              1 + result.accepted_iteration_count,
            "native optimization control pullback count is inconsistent");
    require(instrumentation.explicit_matrix_vmult_calls == 0 &&
              instrumentation.explicit_matrix_tvmult_calls == 0,
            "native optimization unexpectedly used explicit matrix callbacks");
    require(instrumentation.metric_apply_calls == 0 &&
              instrumentation.metric_inverse_apply_calls == 0,
            "native optimization unexpectedly used metric callbacks");
    require(instrumentation.solve_failures == 0,
            "native optimization encountered a solve failure");

    const auto oracle = verification::optimum_oracle(problem);
    Vector nonfinite_control = oracle.control;
    nonfinite_control[0] = std::numeric_limits<double>::quiet_NaN();
    require_rejected(
      [&] {
        verification::require_vector_close(nonfinite_control,
                                           oracle.control,
                                           2.0e-6,
                                           0.0,
                                           "non-finite oracle control probe");
      },
      "oracle control acceptance did not reject a non-finite value");

    Vector outside_oracle_bound = oracle.control;
    outside_oracle_bound[0] += 3.0e-6;
    require_rejected(
      [&] {
        verification::require_vector_close(outside_oracle_bound,
                                           oracle.control,
                                           2.0e-6,
                                           0.0,
                                           "oracle control bound probe");
      },
      "oracle control acceptance did not reject an out-of-bound value");

    Vector nonfinite_gradient = result.derivative.reduced_derivative;
    nonfinite_gradient[0] = std::numeric_limits<double>::quiet_NaN();
    require_rejected(
      [&] {
        verification::require_vector_close(nonfinite_gradient,
                                           result.derivative.reduced_derivative,
                                           1.0e-11,
                                           1.0e-10,
                                           "non-finite gradient probe");
      },
      "gradient acceptance did not reject a non-finite value");

    Vector nonfinite_residual = problem.system_rhs();
    nonfinite_residual[0] = std::numeric_limits<double>::quiet_NaN();
    require_rejected(
      [&] {
        verification::normalized_equation_residual(nonfinite_residual,
                                                   problem.system_rhs());
      },
      "residual acceptance did not reject a non-finite value");

    require(oracle.system_residual <= 1.0e-10 &&
              oracle.stationarity_residual <= 1.0e-10,
            "native optimization oracle audit failed");
    verification::require_vector_close(result.value.control,
                                       oracle.control,
                                       2.0e-6,
                                       0.0,
                                       "native optimization control differs from oracle");

    problem.output_results(result.value.state, artifact_root / "solution.vtk");
    require(instrumentation.output_calls == 1,
            "native optimization output was not written once");

    const auto runtime_counts =
      external_dealii_step4::runtime_counter_snapshot(instrumentation);
    ProblemA      verification_problem(verification_instrumentation);
    NativeReduced verification_reduced(verification_problem,
                                       verification_instrumentation);
    const auto fresh_value =
      verification_reduced.evaluate_value(result.value.control);
    const auto fresh_derivative =
      verification_reduced.augment_derivative(fresh_value);
    const double fresh_gradient_norm =
      fresh_derivative.reduced_derivative.l2_norm();
    const double gradient_difference = vector_difference(
      fresh_derivative.reduced_derivative,
      result.derivative.reduced_derivative);
    require(fresh_gradient_norm <= 1.1e-6,
            "fresh native final gradient exceeds the audit bound");
    verification::require_vector_close(
      fresh_derivative.reduced_derivative,
      result.derivative.reduced_derivative,
      1.0e-11,
      1.0e-10,
      "fresh native final gradient differs from the returned gradient");
    require(verification_instrumentation.assembly_calls == 1 &&
              verification_instrumentation.state_solve_calls == 1 &&
              verification_instrumentation.adjoint_solve_calls == 1 &&
              verification_instrumentation.value_evaluations == 1 &&
              verification_instrumentation.derivative_augmentations == 1,
            "fresh native gradient verification did not recompute one state and adjoint");
    require(external_dealii_step4::runtime_counter_snapshot(instrumentation) ==
              runtime_counts,
            "fresh native gradient verification changed runtime counters");

    std::ofstream gradient_audit(artifact_root / "gradient-audit.csv");
    require(static_cast<bool>(gradient_audit),
            "could not open the native gradient audit artifact");
    gradient_audit
      << "path,returned_gradient_norm,fresh_gradient_norm,"
         "gradient_difference,fresh_state_monitored_residual,"
         "fresh_adjoint_monitored_residual\n"
      << std::setprecision(std::numeric_limits<double>::max_digits10)
      << "native," << result.derivative.reduced_derivative.l2_norm() << ','
      << fresh_gradient_norm << ',' << gradient_difference << ','
      << fresh_value.state_solve.final_residual << ','
      << fresh_derivative.adjoint_solve.final_residual << '\n';
    gradient_audit.flush();
    write_native_optimization_summary(artifact_root,
                                      result,
                                      instrumentation,
                                      oracle);
    evidence.complete();
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
  }

  void
  run_native_trace_failure_contract()
  {
    std::filesystem::path previous_artifact;
    for (unsigned int attempt = 0; attempt < 2; ++attempt)
      {
        std::filesystem::path artifact;
        bool thrown = false;
        try
          {
            run_native_optimization_contract(true, &artifact);
          }
        catch (const std::exception &exception)
          {
            thrown = true;
            require(exception.what() ==
                      std::string("failure after completed native optimization trace"),
                    "native trace failure did not propagate its diagnostic");
          }
        require(thrown, "native trace failure probe did not throw");
        require(artifact != previous_artifact,
                "native trace failure retry reused the run directory");
        const auto read = [&artifact](const char *const filename) {
          std::ifstream input(artifact / filename);
          require(static_cast<bool>(input),
                  std::string("native trace failure lost ") + filename);
          return std::string{std::istreambuf_iterator<char>(input),
                             std::istreambuf_iterator<char>()};
        };
        require(read("status.txt").find("status failed") != std::string::npos,
                "native trace failure did not retain failed status");
        require(read("failure.txt").find(
                  "failure after completed native optimization trace") !=
                  std::string::npos,
                "native trace failure lost its original diagnostic");
        const auto trace = read("trace.csv");
        require(trace.find("\ntrial,") != std::string::npos &&
                  trace.find("\naccepted,") != std::string::npos,
                "native trace failure lost completed trial or accepted records");
        require(read("counters.csv").find("native,state_solve_calls,") !=
                  std::string::npos &&
                  read("solve-records.csv").find("native,success,state") !=
                    std::string::npos,
                "native trace failure lost counters or solve records");
        require(!std::filesystem::exists(artifact / "summary.txt") &&
                  !std::filesystem::exists(artifact / "gradient-audit.csv"),
                "native trace failure ran the later acceptance audits");
        previous_artifact = artifact;
      }
  }

  void
  run_native_optimization_limit_contract()
  {
    Instrumentation instrumentation;
    ProblemA        problem(instrumentation);
    NativeReduced   reduced(problem, instrumentation);
    Vector          initial_control(problem.control_dimension());
    initial_control = 0.0;
    OptimizationPolicy policy =
      external_dealii_step4::frozen_optimization_policy();
    policy.maximum_iterations = 1;

    NativeArmijoSolver solver(reduced, policy);
    const auto result = solver.solve(initial_control);
    require(result.stopping_reason ==
              NativeOptimizationStoppingReason::maximum_iterations,
            "native optimization did not honor iteration-limit precedence");
    require(result.accepted_iteration_count == 1,
            "native optimization iteration-limit probe did not accept one step");
    require(result.gradient_norm_history.size() == 2,
            "native optimization did not check the gradient after acceptance");
    require(std::isfinite(result.derivative.reduced_derivative.l2_norm()) &&
              result.derivative.reduced_derivative.l2_norm() >
                policy.gradient_tolerance,
            "native iteration-limit probe was falsely accepted as converged");
    require(instrumentation.state_solve_calls ==
              1 + result.line_search_trial_count &&
              instrumentation.adjoint_solve_calls ==
                1 + result.accepted_iteration_count,
            "native iteration-limit probe violated staged counts");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"matrix_rhs_solve",
         "nmopt.external_tutorial_step_4.native_matrix_rhs_solve",
         {"dealii", "application", "external", "tutorial", "reuse"},
         60,
         run_matrix_rhs_solve_contract},
        {"fe_access",
         "nmopt.external_tutorial_step_4.native_fe_access",
         {"dealii", "application", "external", "tutorial", "native",
          "problem_b"},
         60,
         run_fe_access_contract},
        {"supplied_state_output",
         "nmopt.external_tutorial_step_4.native_supplied_state_output",
         {"dealii", "application", "external", "tutorial", "reuse"},
         30,
         run_supplied_state_output_contract},
        {"problem_a_operations",
         "nmopt.external_tutorial_step_4.native_problem_a_operations",
         {"dealii", "application", "external", "tutorial", "native", "control"},
         60,
         run_problem_a_operations_contract},
        {"native_reduced",
         "nmopt.external_tutorial_step_4.native_reduced",
         {"dealii", "application", "external", "tutorial", "native", "control"},
         60,
         run_native_reduced_contract},
        {"native_derivatives",
         "nmopt.external_tutorial_step_4.native_derivatives",
         {"dealii", "application", "external", "tutorial", "native", "verification"},
         180,
         run_native_derivative_verification},
        {"native_oracle",
         "nmopt.external_tutorial_step_4.native_oracle",
         {"dealii", "application", "external", "tutorial", "native", "verification"},
         60,
         run_native_oracle_contract},
        {"native_optimization",
         "nmopt.external_tutorial_step_4.native_optimization",
         {"dealii", "application", "external", "tutorial", "native", "optimization"},
         300,
         [] { run_native_optimization_contract(); }},
        {"native_trace_failure",
         "nmopt.external_tutorial_step_4.native_trace_failure",
         {"dealii", "application", "external", "tutorial", "native", "diagnostics"},
         300,
         run_native_trace_failure_contract},
        {"native_optimization_limit",
         "nmopt.external_tutorial_step_4.native_optimization_limit",
         {"dealii", "application", "external", "tutorial", "native", "optimization"},
         120,
         run_native_optimization_limit_contract}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "external Step4 native reuse scenarios passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "external Step4 native reuse test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
