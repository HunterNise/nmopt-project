#include "../../apps/external-dealii/step-4/verification/verification.hpp"
#include "../../apps/external-dealii/step-4/evaluation/native_optimization.hpp"
#include "../../apps/external-dealii/step-4/evaluation/native_problem_b_optimization.hpp"
#include "../../apps/external-dealii/step-4/evaluation/native_problem_b_reduced.hpp"
#include "../../apps/external-dealii/step-4/integration/problem_b_coordinates.hpp"
#include "../../apps/external-dealii/step-4/integration/problem_b_mass.hpp"
#include "../../apps/external-dealii/step-4/integration/problem_b.hpp"
#include "../../apps/external-dealii/step-4/integration/problem_b_metric.hpp"
#include "../../apps/external-dealii/step-4/verification/problem_b_verification.hpp"

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
  using NativeProblemBArmijoSolver =
    external_dealii_step4::NativeProblemBArmijoSolver<2, Step4<2>>;
  using NativeProblemBOptimizationResult =
    external_dealii_step4::NativeProblemBOptimizationResult<2, Step4<2>>;
  using NativeProblemBOptimizationStoppingReason =
    external_dealii_step4::NativeProblemBOptimizationStoppingReason;
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

  std::filesystem::path
  native_problem_b_artifact_root(const char *const scenario)
  {
    auto directory = std::filesystem::current_path();
    while (true)
      {
        if (std::filesystem::exists(
              directory / "apps/external-dealii/step-4/source/upstream/step-4.cc"))
          return external_dealii_step4_test::create_unique_artifact_root(
            directory / "runs/external-dealii/step-4/problem-b/"
                      "native-verification",
            scenario);

        const auto parent = directory.parent_path();
        if (parent == directory)
          break;
        directory = parent;
      }

    throw std::runtime_error(
      "could not locate the Problem B verification artifact root");
  }

  std::filesystem::path
  native_problem_b_optimization_artifact_root(
    const char *const prefix = "native")
  {
    auto directory = std::filesystem::current_path();
    while (true)
      {
        if (std::filesystem::exists(
              directory / "apps/external-dealii/step-4/source/upstream/step-4.cc"))
          return external_dealii_step4_test::create_unique_artifact_root(
            directory / "runs/external-dealii/step-4/problem-b/optimization",
            prefix);

        const auto parent = directory.parent_path();
        if (parent == directory)
          break;
        directory = parent;
      }

    throw std::runtime_error(
      "could not locate the Problem B optimization artifact root");
  }

  void
  write_native_problem_b_optimization_trace(
    const std::filesystem::path &root,
    const NativeProblemBOptimizationResult &result)
  {
    std::ofstream trace(root / "trace.csv");
    require(static_cast<bool>(trace),
            "could not open the native Problem B optimization trace");
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
  write_native_problem_b_metric_trace(
    const std::filesystem::path &root,
    const NativeProblemBOptimizationResult &result)
  {
    std::ofstream trace(root / "metric-solves.csv");
    require(static_cast<bool>(trace),
            "could not open the native Problem B metric trace");
    trace << "check,converged,iterations,initial_residual,final_residual\n"
          << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t index = 0;
         index < result.metric_inverse_evidence.size();
         ++index)
      {
        const auto &evidence = result.metric_inverse_evidence[index];
        trace << index << ',' << evidence.converged << ','
              << evidence.iterations << ',' << evidence.initial_residual
              << ',' << evidence.final_residual << '\n';
      }
  }

  void
  write_native_problem_b_optimization_summary(
    const std::filesystem::path &root,
    const NativeProblemBOptimizationResult &result,
    const Instrumentation &             instrumentation,
    const external_dealii_step4::problem_b_verification::DenseOracle &oracle,
    const double final_gradient_norm,
    const double oracle_control_mass_distance,
    const double oracle_objective,
    const double oracle_control_gradient_norm)
  {
    std::ofstream summary(root / "summary.txt");
    require(static_cast<bool>(summary),
            "could not open the native Problem B optimization summary");
    unsigned int metric_inverse_iterations = 0;
    for (const auto &evidence : result.metric_inverse_evidence)
      metric_inverse_iterations += evidence.iterations;

    summary << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "stopping_reason "
            << external_dealii_step4::
                 native_problem_b_optimization_stopping_reason_name(
                   result.stopping_reason)
            << '\n'
            << "accepted_iterations " << result.accepted_iteration_count << '\n'
            << "line_search_trials " << result.line_search_trial_count << '\n'
            << "final_objective " << result.value.objective << '\n'
            << "final_gradient_norm " << final_gradient_norm << '\n'
            << "oracle_system_residual "
            << oracle.residuals.state_stationarity << '\n'
            << "oracle_control_residual "
            << oracle.residuals.control_stationarity << '\n'
            << "oracle_feasibility_residual "
            << oracle.residuals.feasibility << '\n'
            << "oracle_control_mass_distance "
            << oracle_control_mass_distance << '\n'
            << "oracle_control_gradient_norm "
            << oracle_control_gradient_norm << '\n'
            << "oracle_objective " << oracle_objective << '\n'
            << "final_objective_gap "
            << (result.value.objective - oracle_objective) << '\n'
            << "final_state_solve_iterations "
            << result.value.state_solve.iterations << '\n'
            << "final_state_solve_initial_residual "
            << result.value.state_solve.initial_residual << '\n'
            << "final_state_solve_final_residual "
            << result.value.state_solve.final_residual << '\n'
            << "final_adjoint_solve_iterations "
            << result.derivative.adjoint_solve.iterations << '\n'
            << "final_adjoint_solve_initial_residual "
            << result.derivative.adjoint_solve.initial_residual << '\n'
            << "final_adjoint_solve_final_residual "
            << result.derivative.adjoint_solve.final_residual << '\n'
            << "metric_inverse_iterations " << metric_inverse_iterations << '\n'
            << "state_solve_calls " << instrumentation.state_solve_calls << '\n'
            << "adjoint_solve_calls " << instrumentation.adjoint_solve_calls
            << '\n'
            << "value_evaluations " << instrumentation.value_evaluations << '\n'
            << "derivative_augmentations "
            << instrumentation.derivative_augmentations << '\n'
            << "metric_apply_calls " << instrumentation.metric_apply_calls
            << '\n'
            << "metric_inverse_apply_calls "
            << instrumentation.metric_inverse_apply_calls << '\n';
  }

  void
  write_native_problem_b_failure_summary(
    const std::filesystem::path &root,
    const char *const             outcome,
    const std::string &           detail,
    const std::size_t             state_solve_calls,
    const std::size_t             metric_inverse_apply_calls)
  {
    std::ofstream summary(root / "failure-summary.txt");
    require(static_cast<bool>(summary),
            "could not open the native Problem B failure summary");
    summary << "outcome " << outcome << '\n'
            << "detail " << detail << '\n'
            << "state_solve_calls " << state_solve_calls << '\n'
            << "metric_inverse_apply_calls " << metric_inverse_apply_calls
            << '\n';
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
  run_problem_b_coordinates_contract()
  {
    Step4<2> tutorial;
    tutorial.prepare_for_external_use();

    using Coordinates =
      external_dealii_step4::ProblemBCoordinates<2>;
    const Coordinates coordinates(tutorial.dof_handler_view(),
                                  tutorial.boundary_values_view());
    require(coordinates.full_dimension() == 289,
            "Problem B full coordinate dimension is wrong");
    require(coordinates.free_dimension() == 225,
            "Problem B free coordinate dimension is wrong");
    require(coordinates.boundary_dimension() == 64,
            "Problem B boundary coordinate dimension is wrong");

    const auto &free_indices = coordinates.free_indices();
    require(free_indices.size() == coordinates.free_dimension(),
            "Problem B free index map has the wrong dimension");
    for (std::size_t position = 0; position < free_indices.size(); ++position)
      {
        require(free_indices[position] < coordinates.full_dimension(),
                "Problem B free index map contains an invalid index");
        require(tutorial.boundary_values_view().find(free_indices[position]) ==
                  tutorial.boundary_values_view().end(),
                "Problem B free index map contains a boundary index");
        if (position > 0)
          require(free_indices[position - 1] < free_indices[position],
                  "Problem B free index map is not in native order");
      }

    const auto &lifting = coordinates.lifting();
    require(lifting.size() == coordinates.full_dimension(),
            "Problem B lifting has the wrong dimension");
    for (std::size_t index = 0; index < coordinates.full_dimension(); ++index)
      {
        const auto boundary = tutorial.boundary_values_view().find(index);
        const double expected = boundary == tutorial.boundary_values_view().end() ?
                                  0.0 : boundary->second;
        require_close(lifting[index],
                      expected,
                      0.0,
                      "Problem B lifting has the wrong coefficient");
      }

    Vector free_vector(coordinates.free_dimension());
    for (std::size_t index = 0; index < coordinates.free_dimension(); ++index)
      free_vector[index] =
        0.25 + static_cast<double>(index) /
                  static_cast<double>(coordinates.free_dimension());

    const auto reconstructed = coordinates.reconstruct(free_vector);
    const auto restricted    = coordinates.restrict(reconstructed);
    require(vector_difference(restricted, free_vector) == 0.0,
            "Problem B reconstruction/restriction changed free coordinates");
    for (const auto &[index, value] : tutorial.boundary_values_view())
      require_close(reconstructed[index],
                    value,
                    0.0,
                    "Problem B reconstruction changed a boundary coefficient");

    const auto embedded = coordinates.embed_free(free_vector);
    require(vector_difference(coordinates.restrict(embedded), free_vector) ==
              0.0,
            "Problem B free embedding/restriction changed free coordinates");
    for (const auto &[index, value] : tutorial.boundary_values_view())
      {
        (void)value;
        require_close(embedded[index],
                      0.0,
                      0.0,
                      "Problem B free embedding changed a boundary coordinate");
      }

    Vector zero_free(coordinates.free_dimension());
    zero_free = 0.0;
    require(vector_difference(coordinates.reconstruct(zero_free), lifting) ==
              0.0,
            "Problem B zero reconstruction is not the lifting");
    require_rejected(
      [&coordinates] {
        Vector wrong(coordinates.free_dimension() - 1);
        coordinates.reconstruct(wrong);
      },
      "Problem B reconstruction accepted a wrong free dimension");
  }

  void
  run_problem_b_mass_contract()
  {
    Step4<2> tutorial;
    tutorial.prepare_for_external_use();

    using Coordinates =
      external_dealii_step4::ProblemBCoordinates<2>;
    using Mass = external_dealii_step4::ProblemBMass<2>;
    const Coordinates coordinates(tutorial.dof_handler_view(),
                                  tutorial.boundary_values_view());
    const Mass mass(tutorial.dof_handler_view(), coordinates);

    require(mass.full_dimension() == 289,
            "Problem B mass full dimension is wrong");
    require(mass.free_dimension() == 225,
            "Problem B mass free dimension is wrong");

    const auto &mass_matrix = mass.mass_matrix();
    const auto &coupling    = mass.coupling_matrix();
    require(mass_matrix.m() == 289 && mass_matrix.n() == 289,
            "Problem B mass matrix has the wrong dimensions");
    require(coupling.m() == 225 && coupling.n() == 289,
            "Problem B coupling matrix has the wrong dimensions");
    require(matrix_symmetry_error(mass_matrix) <= 1.0e-14,
            "Problem B mass matrix is not symmetric");

    bool has_positive_diagonal = true;
    bool has_positive_off_diagonal = false;
    for (unsigned int index = 0; index < mass_matrix.m(); ++index)
      {
        has_positive_diagonal = has_positive_diagonal &&
                                mass_matrix.el(index, index) > 0.0;
        for (unsigned int column = 0; column < mass_matrix.n(); ++column)
          if (index != column && mass_matrix.el(index, column) > 0.0)
            has_positive_off_diagonal = true;
      }
    require(has_positive_diagonal,
            "Problem B mass matrix has a non-positive diagonal");
    require(has_positive_off_diagonal,
            "Problem B mass matrix has no consistent off-diagonal coupling");

    Vector ones(mass.full_dimension());
    ones = 1.0;
    const auto row_sums = mass.mass_apply(ones);
    double total_row_sum = 0.0;
    for (unsigned int index = 0; index < row_sums.size(); ++index)
      {
        require(row_sums[index] > 0.0,
                "Problem B mass row sum is not positive");
        total_row_sum += row_sums[index];
      }
    require_close(total_row_sum,
                  4.0,
                  1.0e-12,
                  "Problem B mass row sums do not integrate one");

    for (std::size_t row = 0; row < mass.free_dimension(); ++row)
      for (std::size_t column = 0; column < mass.full_dimension(); ++column)
        require_close(coupling.el(row, column),
                      mass_matrix.el(coordinates.free_indices()[row], column),
                      0.0,
                      "Problem B coupling is not P-transposed mass");

    std::size_t boundary_columns_with_coupling = 0;
    for (const auto &[index, value] : tutorial.boundary_values_view())
      {
        (void)value;
        bool has_coupling = false;
        for (unsigned int row = 0; row < coupling.m(); ++row)
          has_coupling = has_coupling || coupling.el(row, index) != 0.0;
        if (has_coupling)
          ++boundary_columns_with_coupling;
      }
    require(boundary_columns_with_coupling ==
              tutorial.boundary_values_view().size(),
            "Problem B coupling dropped boundary control columns");

    Vector full_vector(mass.full_dimension());
    for (unsigned int index = 0; index < full_vector.size(); ++index)
      full_vector[index] = 0.125 + 0.001 * static_cast<double>(index);
    const auto expected_coupling =
      coordinates.restrict(mass.mass_apply(full_vector));
    require_close(vector_difference(mass.coupling_apply(full_vector),
                                    expected_coupling),
                  0.0,
                  1.0e-14,
                  "Problem B coupling action disagrees with restricted mass action");

    Vector free_vector(mass.free_dimension());
    for (unsigned int index = 0; index < free_vector.size(); ++index)
      free_vector[index] = 0.25 - 0.002 * static_cast<double>(index);
    const auto expected_transpose =
      mass.mass_apply(coordinates.embed_free(free_vector));
    require(vector_difference(mass.coupling_transpose_apply(free_vector),
                              expected_transpose) <= 1.0e-14,
            "Problem B transpose coupling disagrees with mass embedding");

    require_rejected(
      [&mass] {
        Vector wrong(mass.full_dimension() - 1);
        mass.mass_apply(wrong);
      },
      "Problem B mass action accepted a wrong full dimension");
    require_rejected(
      [&mass] {
        Vector wrong(mass.free_dimension() - 1);
        mass.coupling_transpose_apply(wrong);
      },
      "Problem B transpose coupling accepted a wrong free dimension");
  }

  void
  run_problem_b_operations_contract()
  {
    Step4<2> tutorial;
    tutorial.prepare_for_external_use();

    using Problem = external_dealii_step4::ProblemB<2, Step4<2>>;
    using Metric  = external_dealii_step4::ProblemBMetric<2>;
    Problem problem(tutorial);
    Metric  metric(problem.mass());

    require(problem.state_dimension() == 225,
            "Problem B state dimension is wrong");
    require(problem.control_dimension() == 289,
            "Problem B control dimension is wrong");
    require(problem.alpha() == 1.0, "Problem B regularization is not one");
    require(problem.free_system_rhs().size() == problem.state_dimension(),
            "Problem B free system RHS has the wrong dimension");

    const auto &coordinates = problem.coordinates();
    const auto &mass         = problem.mass();

    const auto check_state_solution =
      [&problem, &tutorial, &coordinates, &mass](const Vector &control,
                                                  const char *const name) {
        const auto result = problem.solve_state(control);
        require(result.evidence.converged,
                std::string("Problem B ") + name + " solve did not converge");

        Vector rhs = tutorial.system_rhs_view();
        rhs.add(1.0,
                coordinates.embed_free(mass.coupling_apply(control)));
        Vector lhs(rhs.size());
        tutorial.system_matrix_view().vmult(lhs, result.full_solution);
        require_close(vector_difference(lhs, rhs) /
                        std::max(1.0, rhs.l2_norm()),
                      0.0,
                      1.0e-10,
                      std::string("Problem B ") + name +
                        " full equation residual is too large");
        require_close(vector_difference(
                        coordinates.reconstruct(result.solution),
                        result.full_solution),
                      0.0,
                      1.0e-11,
                      std::string("Problem B ") + name +
                        " full solution disagrees with reconstruction");
        for (const auto &[index, value] : tutorial.boundary_values_view())
          require_close(result.full_solution[index],
                        value,
                        1.0e-11,
                        std::string("Problem B ") + name +
                          " solve changed a prescribed boundary value");
        require(problem.residual(result.solution, control).l2_norm() <= 1.0e-10,
                std::string("Problem B ") + name +
                  " reduced residual is too large");
        return result;
      };

    Vector zero_control(problem.control_dimension());
    zero_control = 0.0;
    const auto zero_state =
      check_state_solution(zero_control, "zero-control state");

    const auto controls = external_dealii_step4::scenario::reduced_controls();
    const auto nonzero_state =
      check_state_solution(controls[1], "nonzero-control state");

    const Vector &state = nonzero_state.solution;
    const Vector &control = controls[1];
    const double objective = problem.objective(state, control);
    require(std::isfinite(objective) && objective >= 0.0,
            "Problem B objective is invalid");
    const auto objective_derivative = problem.objective_derivative(state,
                                                                    control);

    Vector state_tangent(problem.state_dimension());
    for (unsigned int index = 0; index < state_tangent.size(); ++index)
      state_tangent[index] = 0.05 + 0.001 * static_cast<double>(index);
    Vector control_tangent(problem.control_dimension());
    for (unsigned int index = 0; index < control_tangent.size(); ++index)
      control_tangent[index] = -0.025 + 0.0005 * static_cast<double>(index);

    const double step = 1.0e-6;
    Vector state_plus = state;
    Vector state_minus = state;
    state_plus.add(step, state_tangent);
    state_minus.add(-step, state_tangent);
    Vector control_plus = control;
    Vector control_minus = control;
    control_plus.add(step, control_tangent);
    control_minus.add(-step, control_tangent);
    const double centered_objective =
      (problem.objective(state_plus, control_plus) -
       problem.objective(state_minus, control_minus)) /
      (2.0 * step);
    const double analytic_objective =
      objective_derivative.state * state_tangent +
      objective_derivative.control * control_tangent;
    require_close(std::abs(centered_objective - analytic_objective) /
                    std::max(1.0, std::abs(analytic_objective)),
                  0.0,
                  1.0e-8,
                  "Problem B objective derivative disagrees with centered difference");

    Vector test_seed(problem.state_dimension());
    for (unsigned int index = 0; index < test_seed.size(); ++index)
      test_seed[index] = 0.1 - 0.00075 * static_cast<double>(index);
    const auto residual_jvp =
      problem.residual_jvp(state_tangent, control_tangent);
    const auto residual_vjp = problem.residual_vjp(test_seed);
    const double jvp_pairing = residual_jvp * test_seed;
    const double vjp_pairing = residual_vjp.state * state_tangent +
                               residual_vjp.control * control_tangent;
    require_close(std::abs(jvp_pairing - vjp_pairing) /
                    std::max(1.0,
                             std::max(std::abs(jvp_pairing),
                                      std::abs(vjp_pairing))),
                  0.0,
                  1.0e-12,
                  "Problem B residual JVP/VJP pairing failed");

    Vector zero_state_tangent(problem.state_dimension());
    zero_state_tangent = 0.0;
    Vector zero_control_tangent(problem.control_dimension());
    zero_control_tangent = 0.0;
    const auto state_jvp =
      problem.residual_jvp(state_tangent, zero_control_tangent);
    const auto state_vjp = problem.residual_vjp(test_seed);
    require_close(std::abs(state_jvp * test_seed -
                            state_vjp.state * state_tangent) /
                    std::max(1.0, std::abs(state_jvp * test_seed)),
                  0.0,
                  1.0e-12,
                  "Problem B state residual pairing failed");
    const auto control_jvp =
      problem.residual_jvp(zero_state_tangent, control_tangent);
    require_close(std::abs(control_jvp * test_seed -
                            state_vjp.control * control_tangent) /
                    std::max(1.0, std::abs(control_jvp * test_seed)),
                  0.0,
                  1.0e-12,
                  "Problem B control residual pairing failed");

    const auto adjoint =
      problem.solve_adjoint(objective_derivative.state);
    require(adjoint.evidence.converged,
            "Problem B adjoint solve did not converge");
    Vector adjoint_rhs = coordinates.embed_free(objective_derivative.state);
    Vector adjoint_lhs(adjoint_rhs.size());
    tutorial.system_matrix_view().Tvmult(adjoint_lhs,
                                         adjoint.full_solution);
    require_close(vector_difference(adjoint_lhs, adjoint_rhs) /
                    std::max(1.0, adjoint_rhs.l2_norm()),
                  0.0,
                  1.0e-10,
                  "Problem B adjoint equation residual is too large");
    for (const auto &[index, value] : tutorial.boundary_values_view())
      {
        (void)value;
        require_close(adjoint.full_solution[index],
                      0.0,
                      1.0e-11,
                      "Problem B adjoint solve has a boundary value");
      }

    const auto reduced_pullback = problem.residual_vjp(adjoint.solution);
    Vector reduced_gradient = objective_derivative.control;
    reduced_gradient.add(1.0, reduced_pullback.control);
    require(reduced_gradient.size() == problem.control_dimension(),
            "Problem B reduced gradient has the wrong dimension");

    const Vector ramp = external_dealii_step4::scenario::ramp_vector();
    const Vector alternating =
      external_dealii_step4::scenario::alternating_vector();
    const Vector mass_ramp = metric.apply(ramp);
    const Vector mass_alternating = metric.apply(alternating);
    require(ramp * mass_ramp > 0.0 && alternating * mass_alternating > 0.0,
            "Problem B mass metric is not positive");
    require_close(std::abs(ramp * mass_alternating -
                           alternating * mass_ramp),
                  0.0,
                  1.0e-12,
                  "Problem B mass metric is not symmetric");

    const auto recovered_ramp = metric.inverse_apply(mass_ramp);
    const auto recovered_alternating = metric.inverse_apply(mass_alternating);
    require(recovered_ramp.evidence.converged &&
              recovered_alternating.evidence.converged,
            "Problem B mass metric solve did not converge");
    require_close(vector_difference(recovered_ramp.solution, ramp),
                  0.0,
                  1.0e-10,
                  "Problem B mass metric did not recover the ramp");
    require_close(
      vector_difference(recovered_alternating.solution, alternating),
      0.0,
      1.0e-10,
      "Problem B mass metric did not recover the alternating vector");
    require(recovered_ramp.evidence.iterations > 0 &&
              recovered_alternating.evidence.iterations > 0,
            "Problem B mass metric solve did no work");
    require_rejected(
      [&metric] {
        Vector wrong(1);
        metric.apply(wrong);
      },
      "Problem B metric accepted a wrong vector dimension");

    (void)zero_state;
  }

  void
  run_native_problem_b_reduced_contract()
  {
    Instrumentation instrumentation;
    Step4<2>        tutorial;
    tutorial.prepare_for_external_use();

    using Problem = external_dealii_step4::ProblemB<2, Step4<2>>;
    using Reduced = external_dealii_step4::NativeProblemBReduced<2, Step4<2>>;
    Problem problem(tutorial);
    Reduced reduced(problem, instrumentation);

    const auto controls = external_dealii_step4::scenario::reduced_controls();
    const auto value = reduced.evaluate_value(controls.front());
    require(instrumentation.value_evaluations == 1,
            "native Problem B value stage was not counted");
    require(instrumentation.derivative_augmentations == 0,
            "native Problem B derivative ran during value evaluation");
    require(instrumentation.state_solve_calls == 1,
            "native Problem B value stage did not perform one state solve");
    require(instrumentation.adjoint_solve_calls == 0,
            "native Problem B value stage performed an adjoint solve");
    require(instrumentation.objective_calls == 1,
            "native Problem B value stage did not evaluate the objective once");
    require(instrumentation.objective_derivative_calls == 0,
            "native Problem B value stage evaluated objective derivatives");
    require(instrumentation.control_vjp_calls == 0,
            "native Problem B value stage performed a control pullback");
    require(value.state.size() == problem.state_dimension() &&
              value.full_state.size() == problem.control_dimension(),
            "native Problem B value stage returned wrong state dimensions");
    require(std::isfinite(value.objective),
            "native Problem B value stage returned a nonfinite objective");

    const Vector retained_state      = value.state;
    const Vector retained_full_state = value.full_state;
    const auto derivative = reduced.augment_derivative(value);
    require(instrumentation.derivative_augmentations == 1,
            "native Problem B derivative stage was not counted");
    require(instrumentation.state_solve_calls == 1,
            "native Problem B derivative stage repeated the state solve");
    require(instrumentation.adjoint_solve_calls == 1,
            "native Problem B derivative stage did not perform one adjoint solve");
    require(instrumentation.objective_derivative_calls == 1,
            "native Problem B derivative stage did not evaluate objective partials once");
    require(instrumentation.control_vjp_calls == 1,
            "native Problem B derivative stage did not perform one control pullback");
    require(derivative.state_derivative.size() == problem.state_dimension() &&
              derivative.control_derivative.size() == problem.control_dimension() &&
              derivative.adjoint.size() == problem.state_dimension() &&
              derivative.full_adjoint.size() == problem.control_dimension() &&
              derivative.reduced_derivative.size() == problem.control_dimension(),
            "native Problem B derivative stage returned wrong dimensions");
    require(vector_difference(value.state, retained_state) == 0.0 &&
              vector_difference(value.full_state, retained_full_state) == 0.0,
            "native Problem B derivative stage changed retained state data");

    Vector expected_reduced = derivative.control_derivative;
    expected_reduced.add(-1.0, problem.control_vjp(derivative.adjoint));
    require(vector_difference(derivative.reduced_derivative, expected_reduced) ==
              0.0,
            "native Problem B reduced derivative has the wrong composition");

    const auto repeated_before = reduced.evaluate_value(controls[2]);
    const auto intervening      = reduced.evaluate_value(controls[3]);
    const auto repeated_after   = reduced.evaluate_value(controls[2]);
    (void)intervening;
    require(vector_difference(repeated_before.state, repeated_after.state) == 0.0,
            "native Problem B repeated control changed the state");
    require_close(repeated_before.objective,
                  repeated_after.objective,
                  0.0,
                  "native Problem B repeated control changed the objective");
    require(instrumentation.state_solve_calls == 4 &&
              instrumentation.objective_calls == 4 &&
              instrumentation.adjoint_solve_calls == 1,
            "native Problem B repeated values used the wrong staged work");

    require_rejected(
      [&reduced] {
        Vector wrong(1);
        reduced.evaluate_value(wrong);
      },
      "native Problem B evaluator accepted a wrong control dimension");
  }

  void
  run_native_problem_b_oracle_contract()
  {
    Step4<2> tutorial;
    tutorial.prepare_for_external_use();

    using Problem = external_dealii_step4::ProblemB<2, Step4<2>>;
    namespace verification =
      external_dealii_step4::problem_b_verification;
    Problem problem(tutorial);
    const auto operators = verification::make_dense_operators(problem, tutorial);

    require(operators.K.m() == 225 && operators.K.n() == 225,
            "Problem B dense K has the wrong dimensions");
    require(operators.B.m() == 225 && operators.B.n() == 289,
            "Problem B dense B has the wrong dimensions");
    require(operators.M.m() == 289 && operators.M.n() == 289,
            "Problem B dense M has the wrong dimensions");
    require(operators.Q.m() == 225 && operators.Q.n() == 225,
            "Problem B dense Q has the wrong dimensions");

    for (std::size_t index = 0; index < problem.state_dimension(); ++index)
      {
        Vector basis(problem.state_dimension());
        basis = 0.0;
        basis[index] = 1.0;
        const auto round_trip = problem.coordinates().restrict(
          problem.coordinates().embed_free(basis));
        require(verification::vector_difference(round_trip, basis) == 0.0,
                "Problem B dense oracle found a non-identity P-transpose P");
      }

    require(verification::relative_symmetry_error(
              tutorial.system_matrix_view()) <= 1.0e-14,
            "Problem B full system matrix is not symmetric");
    require(verification::relative_symmetry_error(operators.K) <= 1.0e-14,
            "Problem B dense K is not symmetric");
    require(verification::relative_symmetry_error(operators.M) <= 1.0e-14,
            "Problem B dense M is not symmetric");
    verification::require_dense_cholesky(operators.K, "Problem B dense K");
    verification::require_dense_cholesky(operators.M, "Problem B dense M");

    Vector ones(problem.control_dimension());
    ones = 1.0;
    Vector mass_ones(problem.control_dimension());
    operators.M.vmult(mass_ones, ones);
    require_close(ones * mass_ones,
                  4.0,
                  1.0e-12,
                  "Problem B dense mass does not integrate one");

    const auto oracle = verification::dense_kkt_oracle(operators);
    require(oracle.residuals.state_stationarity <= 1.0e-10 &&
              oracle.residuals.control_stationarity <= 1.0e-10 &&
              oracle.residuals.feasibility <= 1.0e-10,
            "Problem B dense KKT oracle has a large residual");

    const auto audit_state_and_adjoint =
      [&problem, &tutorial](const Vector &control, const char *const name) {
        const auto state = problem.solve_state(control);
        require(state.evidence.converged,
                std::string("Problem B ") + name +
                  " state solve did not converge");
        const auto state_audit = verification::audit_state_solution(
          problem, tutorial, state.solution, state.full_solution, control);
        require(state_audit.full_equation.normalized <= 1.0e-10,
                std::string("Problem B ") + name +
                  " full state equation residual is too large");
        require(state_audit.reduced_equation.normalized <= 1.0e-10,
                std::string("Problem B ") + name +
                  " reduced state equation residual is too large");
        require(state_audit.boundary_error <= 1.0e-11,
                std::string("Problem B ") + name +
                  " state boundary error is too large");
        require(state_audit.reconstruction_error <=
                  1.0e-11 + 1.0e-10 * state.full_solution.l2_norm(),
                std::string("Problem B ") + name +
                  " state reconstruction error is too large");

        const auto state_derivative =
          problem.objective_derivative(state.solution, control).state;
        const auto adjoint = problem.solve_adjoint(state_derivative);
        require(adjoint.evidence.converged,
                std::string("Problem B ") + name +
                  " adjoint solve did not converge");
        const auto adjoint_audit = verification::audit_adjoint_solution(
          problem,
          tutorial,
          adjoint.solution,
          adjoint.full_solution,
          state_derivative);
        require(adjoint_audit.full_equation.normalized <= 1.0e-10,
                std::string("Problem B ") + name +
                  " full adjoint equation residual is too large");
        require(adjoint_audit.boundary_error <= 1.0e-11,
                std::string("Problem B ") + name +
                  " adjoint boundary error is too large");
        require(adjoint_audit.homogeneous_reconstruction_error <=
                  1.0e-11 + 1.0e-10 * adjoint.full_solution.l2_norm(),
                std::string("Problem B ") + name +
                  " adjoint reconstruction error is too large");
        return std::pair<decltype(state), decltype(adjoint)>(state, adjoint);
      };

    const auto controls = external_dealii_step4::scenario::reduced_controls();
    for (std::size_t index = 0; index < controls.size(); ++index)
      (void)audit_state_and_adjoint(controls[index], "sample");

    const auto oracle_state_and_adjoint =
      audit_state_and_adjoint(oracle.control, "oracle");
    verification::require_vector_close(oracle_state_and_adjoint.first.solution,
                                       oracle.state,
                                       1.0e-11,
                                       1.0e-10,
                                       "native state differs from Problem B dense oracle");
    verification::require_vector_close(
      oracle_state_and_adjoint.second.solution,
      oracle.adjoint,
      1.0e-11,
      1.0e-10,
                                       "native adjoint differs from Problem B dense oracle");
  }

  void
  run_native_problem_b_derivative_metric_contract()
  {
    Instrumentation instrumentation;
    const auto artifact_root = native_problem_b_artifact_root("derivative-metric");
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact_root,
      "native_problem_b_derivative_metric",
      {{"native", &instrumentation}});
    try
      {
        Step4<2> tutorial;
        tutorial.prepare_for_external_use();

        using Problem = external_dealii_step4::ProblemB<2, Step4<2>>;
        using Metric  = external_dealii_step4::ProblemBMetric<2>;
        using Reduced =
          external_dealii_step4::NativeProblemBReduced<2, Step4<2>>;
        using ValueEvaluation = typename Reduced::ValueEvaluation;
        using DerivativeEvaluation = typename Reduced::DerivativeEvaluation;
        namespace problem_b_verification =
          external_dealii_step4::problem_b_verification;

        Problem problem(tutorial);
        Metric  metric(problem.mass());
        Reduced reduced(problem, instrumentation);
        std::ofstream output(artifact_root / "derivative-metric.csv");
        require(static_cast<bool>(output),
                "could not open the Problem B derivative/metric trace");
        output << "record,index,direction,step,analytic,value,error,bound\n"
               << std::setprecision(std::numeric_limits<double>::max_digits10);

        Vector state_tangent(problem.state_dimension());
        Vector control_tangent(problem.control_dimension());
        Vector test_seed(problem.state_dimension());
        for (unsigned int index = 0; index < state_tangent.size(); ++index)
          {
            state_tangent[index] = 0.05 + 0.001 * static_cast<double>(index);
            test_seed[index] = 0.1 - 0.00075 * static_cast<double>(index);
          }
        for (unsigned int index = 0; index < control_tangent.size(); ++index)
          control_tangent[index] = -0.025 + 0.0005 * static_cast<double>(index);

        Vector ramp = external_dealii_step4::scenario::ramp_vector();
        Vector alternating =
          external_dealii_step4::scenario::alternating_vector();
        const Vector ramp_raw = ramp;
        const Vector alternating_raw = alternating;
        ramp /= problem_b_verification::mass_norm(metric, ramp);
        alternating /= problem_b_verification::mass_norm(metric, alternating);

        const double ramp_metric_pairing =
          ramp_raw * metric.apply(ramp_raw);
        const double alternating_metric_pairing =
          alternating_raw * metric.apply(alternating_raw);
        require(ramp_metric_pairing > 0.0 &&
                  alternating_metric_pairing > 0.0,
                "Problem B mass metric is not positive");
        require_close(std::abs(ramp_raw * metric.apply(alternating_raw) -
                               alternating_raw * metric.apply(ramp_raw)),
                      0.0,
                      1.0e-12,
                      "Problem B mass metric is not symmetric");

        const auto audit_metric_vector =
          [&metric, &output](const Vector &vector, const char *const name) {
            const auto rhs = metric.apply(vector);
            const auto inverse = metric.inverse_apply(rhs);
            require(inverse.evidence.converged,
                    std::string("Problem B ") + name +
                      " metric inverse did not converge");
            const auto residual = problem_b_verification::normalized_equation_residual(
              metric.apply(inverse.solution), rhs);
            require(residual.normalized <= 1.0e-10,
                    std::string("Problem B ") + name +
                      " metric equation residual is too large");
            problem_b_verification::require_vector_close(
              inverse.solution,
              vector,
              1.0e-11,
              1.0e-10,
              std::string("Problem B ") + name +
                " metric inverse changed the vector");
            output << "metric_inverse,," << name << ','
                   << inverse.evidence.iterations << ",,"
                   << residual.normalized << ",,1e-10\n";
          };
        audit_metric_vector(ramp_raw, "ramp");
        audit_metric_vector(alternating_raw, "alternating");

        const auto controls =
          external_dealii_step4::scenario::reduced_controls();
        std::vector<std::pair<ValueEvaluation, DerivativeEvaluation>> samples;
        samples.reserve(controls.size());

        const auto audit_sample =
          [&problem,
           &metric,
           &reduced,
           &tutorial,
           &output,
           &state_tangent,
           &control_tangent,
           &test_seed](const Vector &control, const std::size_t index) {
            const auto value = reduced.evaluate_value(control);
            const auto derivative = reduced.augment_derivative(value);
            const auto objective_derivative =
              problem.objective_derivative(value.state, control);
            problem_b_verification::require_vector_close(
              derivative.state_derivative,
              objective_derivative.state,
              0.0,
              0.0,
              "Problem B state objective derivative changed");
            problem_b_verification::require_vector_close(
              derivative.control_derivative,
              objective_derivative.control,
              0.0,
              0.0,
              "Problem B control objective derivative changed");

            const auto centered_jvp =
              problem_b_verification::centered_residual_jvp(
                problem,
                value.state,
                control,
                state_tangent,
                control_tangent,
                1.0e-6);
            const auto analytic_jvp =
              problem.residual_jvp(state_tangent, control_tangent);
            const double residual_jvp_error =
              problem_b_verification::scaled_vector_error(centered_jvp,
                                                          analytic_jvp);
            const double analytic_objective =
              objective_derivative.state * state_tangent +
              objective_derivative.control * control_tangent;
            const double centered_objective =
              problem_b_verification::centered_objective_directional_derivative(
                problem,
                value.state,
                control,
                state_tangent,
                control_tangent,
                1.0e-6);
            const double objective_error =
              problem_b_verification::scaled_scalar_error(
                centered_objective, analytic_objective);
            require(residual_jvp_error <= 1.0e-8,
                    "Problem B residual JVP centered difference failed");
            require(objective_error <= 1.0e-8,
                    "Problem B objective derivative centered difference failed");

            const auto jvp =
              problem.residual_jvp(state_tangent, control_tangent);
            const auto vjp = problem.residual_vjp(test_seed);
            const double full_left = jvp * test_seed;
            const double full_right =
              vjp.state * state_tangent + vjp.control * control_tangent;
            Vector zero_control(problem.control_dimension());
            zero_control = 0.0;
            const double state_left =
              problem.residual_jvp(state_tangent, zero_control) * test_seed;
            const double state_right = vjp.state * state_tangent;
            Vector zero_state(problem.state_dimension());
            zero_state = 0.0;
            const double control_left =
              problem.residual_jvp(zero_state, control_tangent) * test_seed;
            const double control_right = vjp.control * control_tangent;
            const double full_pairing_error =
              problem_b_verification::scaled_scalar_error(full_left,
                                                          full_right);
            const double state_pairing_error =
              problem_b_verification::scaled_scalar_error(state_left,
                                                          state_right);
            const double control_pairing_error =
              problem_b_verification::scaled_scalar_error(control_left,
                                                          control_right);
            require(full_pairing_error <= 1.0e-12 &&
                      state_pairing_error <= 1.0e-12 &&
                      control_pairing_error <= 1.0e-12,
                    "Problem B residual JVP/VJP pairing failed");

            const auto physical_state =
              problem.coordinates().reconstruct(value.state);
            const double quadrature_objective =
              0.5 * problem_b_verification::cell_quadrature_pairing<2>(
                      tutorial, physical_state, physical_state) +
              0.5 * problem_b_verification::cell_quadrature_pairing<2>(
                      tutorial, control, control);
            require_close(problem.objective(value.state, control),
                          quadrature_objective,
                          1.0e-12,
                          "Problem B objective disagrees with cell quadrature");

            const auto metric_gradient =
              metric.inverse_apply(derivative.reduced_derivative);
            require(metric_gradient.evidence.converged,
                    "Problem B metric gradient solve did not converge");
            const auto metric_residual =
              problem_b_verification::normalized_equation_residual(
                metric.apply(metric_gradient.solution),
                derivative.reduced_derivative);
            require(metric_residual.normalized <= 1.0e-10,
                    "Problem B metric gradient residual is too large");

            const double objective_value =
              problem.objective(value.state, control);
            const double objective_quadrature_error =
              std::abs(objective_value - quadrature_objective);
            output << "derivative_jvp," << index << ",,1e-6,,"
                   << residual_jvp_error << ",,1e-8\n";
            output << "objective_derivative," << index << ",,1e-6,"
                   << analytic_objective << ',' << centered_objective << ','
                   << objective_error << ",1e-8\n";
            output << "pairing_full," << index << ",,1e-6," << full_left
                   << ',' << full_right << ',' << full_pairing_error
                   << ",1e-12\n";
            output << "pairing_state," << index << ",,1e-6," << state_left
                   << ',' << state_right << ',' << state_pairing_error
                   << ",1e-12\n";
            output << "pairing_control," << index << ",,1e-6," << control_left
                   << ',' << control_right << ',' << control_pairing_error
                   << ",1e-12\n";
            output << "objective_quadrature," << index << ",,,"
                   << objective_value << ',' << quadrature_objective << ','
                   << objective_quadrature_error << ",1e-12\n";
            output << "metric_gradient," << index << ",,"
                   << metric_gradient.evidence.iterations << ",,"
                   << metric_residual.normalized << ",,1e-10\n";
            return std::make_pair(value, derivative);
          };

        for (std::size_t index = 0; index < controls.size(); ++index)
          samples.push_back(audit_sample(controls[index], index));

        const auto repeated_ramp = reduced.evaluate_value(controls[2]);
        require(problem_b_verification::vector_difference(
                  repeated_ramp.state, samples[2].first.state) == 0.0,
                "Problem B repeated ramp changed the state");
        require_close(repeated_ramp.objective,
                      samples[2].first.objective,
                      0.0,
                      "Problem B repeated ramp changed the objective");

        Vector mass_left(problem.control_dimension());
        Vector mass_right(problem.control_dimension());
        for (unsigned int index = 0; index < mass_left.size(); ++index)
          {
            mass_left[index] = 0.125 + 0.001 * static_cast<double>(index);
            mass_right[index] = -0.05 + 0.0007 * static_cast<double>(index);
          }
        const double assembled_pairing =
          mass_left * problem.mass().mass_apply(mass_right);
        const double quadrature_pairing =
          problem_b_verification::cell_quadrature_pairing<2>(
            tutorial, mass_left, mass_right);
        require_close(assembled_pairing,
                      quadrature_pairing,
                      1.0e-12,
                      "Problem B mass pairing disagrees with cell quadrature");

        Vector free_test(problem.state_dimension());
        for (unsigned int index = 0; index < free_test.size(); ++index)
          free_test[index] = 0.2 - 0.0009 * static_cast<double>(index);
        const double assembled_weak_action =
          free_test * problem.mass().coupling_apply(mass_left);
        const double quadrature_weak_action =
          problem_b_verification::cell_quadrature_pairing<2>(
            tutorial,
            mass_left,
            problem.coordinates().embed_free(free_test));
        require_close(assembled_weak_action,
                      quadrature_weak_action,
                      1.0e-12,
                      "Problem B weak coupling disagrees with cell quadrature");

        const std::vector<double> finite_difference_steps{
          1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6};
        const Vector directions[] = {ramp, alternating};
        const char *const direction_names[] = {"ramp", "alternating"};
        for (std::size_t sample_index = 0; sample_index < samples.size();
             ++sample_index)
          for (std::size_t direction_index = 0; direction_index < 2;
               ++direction_index)
            {
              const double analytic =
                samples[sample_index].second.reduced_derivative *
                directions[direction_index];
              bool passes[5] = {false, false, false, false, false};
              for (std::size_t step_index = 0;
                   step_index < finite_difference_steps.size();
                   ++step_index)
                {
                  const double step = finite_difference_steps[step_index];
                  Vector plus = samples[sample_index].first.control;
                  Vector minus = samples[sample_index].first.control;
                  plus.add(step, directions[direction_index]);
                  minus.add(-step, directions[direction_index]);
                  const auto plus_value = reduced.evaluate_value(plus);
                  const auto minus_value = reduced.evaluate_value(minus);
                  const double centered =
                    (plus_value.objective - minus_value.objective) /
                    (2.0 * step);
                  const double error = std::abs(centered - analytic);
                  const double bound =
                    1.0e-7 * std::max(1.0, std::abs(analytic));
                  passes[step_index] = std::isfinite(error) && error <= bound;
                  output << "reduced_fd," << sample_index << ','
                         << direction_names[direction_index] << ',' << step
                         << ',' << analytic << ',' << centered << ',' << error
                         << ',' << bound
                         << '\n';
                }
              bool adjacent_passes = false;
              for (std::size_t step_index = 1;
                   step_index < finite_difference_steps.size();
                   ++step_index)
                adjacent_passes = adjacent_passes ||
                                  (passes[step_index - 1] && passes[step_index]);
              require(adjacent_passes,
                      "Problem B reduced centered differences lack adjacent usable steps");
            }

        const auto base_value = reduced.evaluate_value(controls.front());
        const auto base_derivative = reduced.augment_derivative(base_value);
        const auto repeated_zero = reduced.evaluate_value(controls.front());
        const double repeat_variation =
          std::abs(repeated_zero.objective - base_value.objective);
        const std::vector<double> taylor_steps{0.1, 0.05, 0.025};
        for (std::size_t direction_index = 0; direction_index < 2;
             ++direction_index)
          {
            const double slope = base_derivative.reduced_derivative *
                                 directions[direction_index];
            std::vector<double> remainders;
            for (const double step : taylor_steps)
              {
                Vector control = controls.front();
                control.add(step, directions[direction_index]);
                const auto trial = reduced.evaluate_value(control);
                const double remainder = trial.objective -
                                         base_value.objective - step * slope;
                const double ratio = remainders.empty() ?
                                       std::numeric_limits<double>::quiet_NaN() :
                                       remainders.back() / remainder;
                remainders.push_back(remainder);
                output << "taylor,0," << direction_names[direction_index]
                       << ',' << step << ',' << remainder << ',' << ratio
                       << ',' << repeat_variation << ",\n";
                require(std::isfinite(remainder) &&
                          remainder > 10.0 * repeat_variation,
                        "Problem B Taylor remainder is not positive and resolved");
              }
            for (std::size_t step_index = 1; step_index < remainders.size();
                 ++step_index)
              {
                const double ratio =
                  remainders[step_index - 1] / remainders[step_index];
                require(ratio >= 3.5 && ratio <= 4.5,
                        "Problem B Taylor remainder does not have quadratic halving");
              }
          }

        output.flush();
        evidence.complete();
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
  }

  void
  run_native_problem_b_optimization_contract()
  {
    Instrumentation instrumentation;
    Step4<2>        tutorial;
    tutorial.prepare_for_external_use();

    using Problem = external_dealii_step4::ProblemB<2, Step4<2>>;
    using Metric  = external_dealii_step4::ProblemBMetric<2>;
    using Reduced = external_dealii_step4::NativeProblemBReduced<2, Step4<2>>;

    Problem problem(tutorial);
    Metric  metric(problem.mass());
    Reduced reduced(problem, instrumentation);
    Vector initial_control(problem.control_dimension());
    initial_control = 0.0;

    OptimizationPolicy policy =
      external_dealii_step4::frozen_optimization_policy();
    policy.maximum_iterations = 1;

    NativeProblemBArmijoSolver solver(reduced, metric, instrumentation, policy);
    const auto result = solver.solve(initial_control);

    require(result.stopping_reason ==
              NativeProblemBOptimizationStoppingReason::maximum_iterations,
            "native Problem B optimization did not honor iteration-limit "
            "precedence");
    require(result.accepted_iteration_count == 1,
            "native Problem B optimization did not accept one step");
    require(result.gradient_norm_history.size() == 2,
            "native Problem B optimization did not check the gradient after "
            "acceptance");
    require(result.metric_inverse_evidence.size() == 2,
            "native Problem B optimization did not retain both metric solves");
    for (const auto &evidence : result.metric_inverse_evidence)
      require(evidence.converged,
              "native Problem B optimization metric solve did not converge");

    require(result.trial_records.size() == result.line_search_trial_count &&
              !result.trial_records.empty() &&
              result.trial_records.back().accepted,
            "native Problem B optimization did not retain an accepted trial");
    require(result.accepted_iterations.size() == 1 &&
              result.accepted_iterations.front().actual_step_norm > 0.0,
            "native Problem B optimization did not record its metric step");
    require(result.gradient_norm_history.back() > policy.gradient_tolerance,
            "native Problem B iteration-limit probe was falsely converged");

    require(instrumentation.state_solve_calls ==
              1 + result.line_search_trial_count &&
              instrumentation.adjoint_solve_calls ==
                1 + result.accepted_iteration_count,
            "native Problem B optimization violated the staged solve schedule");
    require(instrumentation.value_evaluations ==
              1 + result.line_search_trial_count &&
              instrumentation.derivative_augmentations ==
                1 + result.accepted_iteration_count,
            "native Problem B optimization violated the staged evaluation "
            "schedule");
    require(instrumentation.objective_calls ==
              1 + result.line_search_trial_count &&
              instrumentation.objective_derivative_calls ==
                1 + result.accepted_iteration_count &&
              instrumentation.control_vjp_calls ==
                1 + result.accepted_iteration_count,
            "native Problem B optimization objective schedule is inconsistent");
    require(instrumentation.metric_inverse_apply_calls ==
              1 + result.accepted_iteration_count &&
              instrumentation.metric_apply_calls ==
                1 + 2 * result.accepted_iteration_count,
            "native Problem B optimization metric schedule is inconsistent");
  }

  void
  run_native_problem_b_optimization_evidence()
  {
    Instrumentation instrumentation;
    const auto artifact_root =
      native_problem_b_optimization_artifact_root();
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact_root,
      "native_problem_b_optimization",
      {{"native", &instrumentation}});
    try
      {
        Step4<2> tutorial;
        tutorial.prepare_for_external_use();

        using Problem = external_dealii_step4::ProblemB<2, Step4<2>>;
        using Metric  = external_dealii_step4::ProblemBMetric<2>;
        using Reduced =
          external_dealii_step4::NativeProblemBReduced<2, Step4<2>>;

        Problem problem(tutorial);
        Metric  metric(problem.mass());
        Reduced reduced(problem, instrumentation);
        Vector initial_control(problem.control_dimension());
        initial_control = 0.0;

        const auto policy =
          external_dealii_step4::frozen_optimization_policy();
        NativeProblemBArmijoSolver solver(reduced, metric, instrumentation, policy);
        const auto result = solver.solve(initial_control);

        write_native_problem_b_optimization_trace(artifact_root, result);
        write_native_problem_b_metric_trace(artifact_root, result);

        require(result.stopping_reason ==
                  NativeProblemBOptimizationStoppingReason::gradient_tolerance,
                "native Problem B optimization did not stop by gradient "
                "tolerance");
        require(result.accepted_iteration_count > 0,
                "native Problem B optimization accepted no iterations");
        require(result.objective_history.size() ==
                  result.accepted_iteration_count + 1 &&
                  result.accepted_iterations.size() ==
                    result.accepted_iteration_count,
                "native Problem B optimization histories have inconsistent "
                "sizes");
        require(result.trial_records.size() ==
                  result.line_search_trial_count,
                "native Problem B optimization trial history has the wrong "
                "size");
        require(result.metric_inverse_evidence.size() ==
                  result.accepted_iteration_count + 1,
                "native Problem B optimization metric history has the wrong "
                "size");
        for (const auto &metric_evidence : result.metric_inverse_evidence)
          require(metric_evidence.converged &&
                    std::isfinite(metric_evidence.final_residual),
                  "native Problem B optimization had an invalid metric solve");

        std::size_t accepted_trials = 0;
        for (const auto &trial : result.trial_records)
          {
            if (trial.accepted)
              {
                ++accepted_trials;
                require(trial.objective_finite && trial.slope_negative,
                        "native Problem B accepted an invalid trial");
              }
          }
        require(accepted_trials == result.accepted_iteration_count,
                "native Problem B accepted-trial count is inconsistent");
        for (std::size_t index = 1; index < result.objective_history.size();
             ++index)
          require(result.objective_history[index] <=
                    result.objective_history[index - 1],
                  "native Problem B objective increased after acceptance");

        require(instrumentation.state_solve_calls ==
                  1 + result.line_search_trial_count &&
                  instrumentation.adjoint_solve_calls ==
                    1 + result.accepted_iteration_count &&
                  instrumentation.value_evaluations ==
                    1 + result.line_search_trial_count &&
                  instrumentation.derivative_augmentations ==
                    1 + result.accepted_iteration_count,
                "native Problem B optimization violated the staged schedule");
        require(instrumentation.objective_calls ==
                  1 + result.line_search_trial_count &&
                  instrumentation.objective_derivative_calls ==
                    1 + result.accepted_iteration_count &&
                  instrumentation.control_vjp_calls ==
                    1 + result.accepted_iteration_count,
                "native Problem B optimization objective schedule is "
                "inconsistent");
        require(instrumentation.metric_inverse_apply_calls ==
                  1 + result.accepted_iteration_count &&
                  instrumentation.metric_apply_calls ==
                    1 + 2 * result.accepted_iteration_count,
                "native Problem B optimization metric schedule is "
                "inconsistent");
        require(instrumentation.residual_calls == 0 &&
                  instrumentation.residual_jvp_calls == 0 &&
                  instrumentation.residual_vjp_calls == 0 &&
                  instrumentation.explicit_matrix_vmult_calls == 0 &&
                  instrumentation.explicit_matrix_tvmult_calls == 0 &&
                  instrumentation.solve_failures == 0,
                "native Problem B optimization used an unexpected operation");

        tutorial.output_results(result.value.full_state,
                                artifact_root / "solution.vtk");
        require(std::filesystem::exists(artifact_root / "solution.vtk"),
                "native Problem B optimization did not retain its state output");

        Step4<2> verification_tutorial;
        verification_tutorial.prepare_for_external_use();
        Problem verification_problem(verification_tutorial);
        Metric  verification_metric(verification_problem.mass());
        Instrumentation verification_instrumentation;
        Reduced verification_reduced(verification_problem,
                                     verification_instrumentation);
        const auto fresh_value =
          verification_reduced.evaluate_value(result.value.control);
        const auto fresh_derivative =
          verification_reduced.augment_derivative(fresh_value);

        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_value.state,
          result.value.state,
          1.0e-11,
          1.0e-10,
          "fresh native Problem B final state differs");
        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_value.full_state,
          result.value.full_state,
          1.0e-11,
          1.0e-10,
          "fresh native Problem B full state differs");
        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_derivative.adjoint,
          result.derivative.adjoint,
          1.0e-11,
          1.0e-10,
          "fresh native Problem B final adjoint differs");
        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_derivative.full_adjoint,
          result.derivative.full_adjoint,
          1.0e-11,
          1.0e-10,
          "fresh native Problem B full adjoint differs");
        external_dealii_step4::problem_b_verification::require_vector_close(
          fresh_derivative.reduced_derivative,
          result.derivative.reduced_derivative,
          1.0e-11,
          1.0e-10,
          "fresh native Problem B final covector differs");
        require_close(fresh_value.objective,
                      result.value.objective,
                      1.0e-12,
                      "fresh native Problem B final objective differs");
        require(fresh_value.state_solve.converged &&
                  fresh_derivative.adjoint_solve.converged,
                "fresh native Problem B final solves did not converge");

        const auto state_audit =
          external_dealii_step4::problem_b_verification::audit_state_solution(
            verification_problem,
            verification_tutorial,
            fresh_value.state,
            fresh_value.full_state,
            fresh_value.control);
        const auto adjoint_audit =
          external_dealii_step4::problem_b_verification::audit_adjoint_solution(
            verification_problem,
            verification_tutorial,
            fresh_derivative.adjoint,
            fresh_derivative.full_adjoint,
            fresh_derivative.state_derivative);
        require(state_audit.full_equation.normalized <= 1.0e-10 &&
                  state_audit.reduced_equation.normalized <= 1.0e-10 &&
                  state_audit.reconstruction_error <= 1.0e-11 &&
                  state_audit.boundary_error <= 1.0e-11,
                "native Problem B final state audit failed");
        require(adjoint_audit.full_equation.normalized <= 1.0e-10 &&
                  adjoint_audit.homogeneous_reconstruction_error <= 1.0e-11 &&
                  adjoint_audit.boundary_error <= 1.0e-11,
                "native Problem B final adjoint audit failed");

        const auto operators =
          external_dealii_step4::problem_b_verification::make_dense_operators(
            verification_problem,
            verification_tutorial);
        const double final_gradient_norm =
          external_dealii_step4::problem_b_verification::dense_mass_norm(
            operators.M,
            fresh_derivative.reduced_derivative);
        require(final_gradient_norm <= 1.1e-6,
                "native Problem B final metric gradient exceeds the audit "
                "bound");
        require_close(result.gradient_norm_history.back(),
                      final_gradient_norm,
                      1.0e-12,
                      "native Problem B returned and fresh gradient norms "
                      "differ");

        const auto oracle =
          external_dealii_step4::problem_b_verification::dense_kkt_oracle(
            operators);
        require(oracle.residuals.state_stationarity <= 1.0e-10 &&
                  oracle.residuals.control_stationarity <= 1.0e-10 &&
                  oracle.residuals.feasibility <= 1.0e-10,
                "native Problem B optimization dense oracle failed");

        Step4<2> oracle_tutorial;
        oracle_tutorial.prepare_for_external_use();
        Problem oracle_problem(oracle_tutorial);
        Instrumentation oracle_instrumentation;
        Reduced oracle_reduced(oracle_problem, oracle_instrumentation);
        const auto oracle_value =
          oracle_reduced.evaluate_value(oracle.control);
        const auto oracle_derivative =
          oracle_reduced.augment_derivative(oracle_value);
        const double oracle_control_gradient_norm =
          external_dealii_step4::problem_b_verification::dense_mass_norm(
            operators.M,
            oracle_derivative.reduced_derivative);
        require(oracle_control_gradient_norm <= 1.0e-8,
                "native Problem B oracle control gradient is not zero");

        Vector control_difference = fresh_value.control;
        control_difference.add(-1.0, oracle.control);
        const double oracle_control_mass_distance =
          external_dealii_step4::problem_b_verification::mass_norm(
            verification_metric,
            control_difference);
        require(oracle_control_mass_distance <= 2.0e-6,
                "native Problem B final control differs from the oracle");

        const double oracle_objective =
          verification_problem.objective(oracle.state, oracle.control);
        require(std::isfinite(oracle_objective) &&
                  std::isfinite(result.value.objective - oracle_objective),
                "native Problem B objective gap is not finite");

        std::ofstream gradient_audit(artifact_root / "gradient-audit.csv");
        require(static_cast<bool>(gradient_audit),
                "could not open the native Problem B gradient audit");
        gradient_audit
          << "path,returned_gradient_norm,fresh_gradient_norm,"
             "gradient_difference,control_mass_distance,"
             "state_residual,adjoint_residual\n"
          << std::setprecision(std::numeric_limits<double>::max_digits10)
          << "native," << result.gradient_norm_history.back() << ','
          << final_gradient_norm << ','
          << external_dealii_step4::problem_b_verification::vector_difference(
               fresh_derivative.reduced_derivative,
               result.derivative.reduced_derivative)
          << ',' << oracle_control_mass_distance << ','
          << state_audit.full_equation.normalized << ','
          << adjoint_audit.full_equation.normalized << '\n';
        gradient_audit.flush();
        write_native_problem_b_optimization_summary(
          artifact_root,
          result,
          instrumentation,
          oracle,
          final_gradient_norm,
          oracle_control_mass_distance,
          oracle_objective,
          oracle_control_gradient_norm);
        evidence.complete();
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
  }

  void
  run_native_problem_b_line_search_failure_contract()
  {
    Instrumentation instrumentation;
    const auto artifact_root =
      native_problem_b_optimization_artifact_root("line-search-failure");
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact_root,
      "native_problem_b_line_search_failure",
      {{"native", &instrumentation}});
    try
      {
        Step4<2> tutorial;
        tutorial.prepare_for_external_use();

        using Problem = external_dealii_step4::ProblemB<2, Step4<2>>;
        using Metric  = external_dealii_step4::ProblemBMetric<2>;
        using Reduced =
          external_dealii_step4::NativeProblemBReduced<2, Step4<2>>;

        Problem problem(tutorial);
        Metric  metric(problem.mass());
        Reduced reduced(problem, instrumentation);
        Vector initial_control(problem.control_dimension());
        initial_control = 0.0;

        auto policy = external_dealii_step4::frozen_optimization_policy();
        policy.maximum_line_search_trials = 1;
        policy.initial_step_length = 1.0e6;
        NativeProblemBArmijoSolver solver(
          reduced, metric, instrumentation, policy);
        const auto result = solver.solve(initial_control);

        write_native_problem_b_optimization_trace(artifact_root, result);
        write_native_problem_b_metric_trace(artifact_root, result);
        write_native_problem_b_failure_summary(
          artifact_root,
          "line_search_failure",
          "bounded Armijo trial was rejected",
          instrumentation.state_solve_calls,
          instrumentation.metric_inverse_apply_calls);

        require(result.stopping_reason ==
                  NativeProblemBOptimizationStoppingReason::line_search_failure,
                "native Problem B line-search failure was not reported");
        require(result.accepted_iteration_count == 0 &&
                  result.line_search_trial_count == 1,
                "native Problem B line-search failure changed the accepted "
                "iteration or trial count");
        require(result.trial_records.size() == 1 &&
                  !result.trial_records.front().accepted,
                "native Problem B line-search failure lost its rejected trial");
        require(instrumentation.state_solve_calls == 2 &&
                  instrumentation.adjoint_solve_calls == 1 &&
                  instrumentation.value_evaluations == 2 &&
                  instrumentation.derivative_augmentations == 1 &&
                  instrumentation.metric_inverse_apply_calls == 1 &&
                  instrumentation.metric_apply_calls == 1,
                "native Problem B line-search failure violated the staged "
                "failure schedule");

        evidence.fail("native Problem B line search returned line_search_failure");
        const auto read = [&artifact_root](const char *const filename) {
          std::ifstream input(artifact_root / filename);
          require(static_cast<bool>(input),
                  std::string("native Problem B line-search failure lost ") +
                    filename);
          return std::string{std::istreambuf_iterator<char>(input),
                             std::istreambuf_iterator<char>()};
        };
        require(read("status.txt").find("status failed") != std::string::npos,
                "native Problem B line-search failure did not retain failed status");
        require(read("failure.txt").find("line_search_failure") !=
                  std::string::npos,
                "native Problem B line-search failure lost its diagnostic");
        require(read("trace.csv").find("\ntrial,") != std::string::npos &&
                  read("metric-solves.csv").find("\n0,1,") !=
                    std::string::npos &&
                  read("failure-summary.txt").find("outcome line_search_failure") !=
                    std::string::npos,
                "native Problem B line-search failure lost available traces");
        require(read("counters.csv").find("native,state_solve_calls,2") !=
                  std::string::npos &&
                  read("solve-records.csv").find("native,success,state") !=
                    std::string::npos,
                "native Problem B line-search failure lost counters or solves");
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
      }
  }

  void
  run_native_problem_b_nonfinite_contract()
  {
    Instrumentation instrumentation;
    const auto artifact_root =
      native_problem_b_optimization_artifact_root("nonfinite-inputs");
    external_dealii_step4_test::EvidenceGuard evidence(
      artifact_root,
      "native_problem_b_nonfinite",
      {{"native", &instrumentation}});
    try
      {
        Step4<2> tutorial;
        tutorial.prepare_for_external_use();

        using Problem = external_dealii_step4::ProblemB<2, Step4<2>>;
        using Metric  = external_dealii_step4::ProblemBMetric<2>;
        using Reduced =
          external_dealii_step4::NativeProblemBReduced<2, Step4<2>>;

        Problem problem(tutorial);
        Metric  metric(problem.mass());
        Reduced reduced(problem, instrumentation);

        Vector nonfinite_metric_rhs(problem.control_dimension());
        nonfinite_metric_rhs = 0.0;
        nonfinite_metric_rhs[0] = std::numeric_limits<double>::quiet_NaN();
        std::string metric_message;
        try
          {
            (void)metric.inverse_apply(nonfinite_metric_rhs);
          }
        catch (const std::exception &exception)
          {
            metric_message = exception.what();
          }
        require(metric_message.find("non-finite") != std::string::npos,
                "native Problem B metric did not reject a non-finite RHS");

        Vector nonfinite_control(problem.control_dimension());
        nonfinite_control = 0.0;
        nonfinite_control[0] = std::numeric_limits<double>::quiet_NaN();
        auto policy = external_dealii_step4::frozen_optimization_policy();
        NativeProblemBArmijoSolver solver(
          reduced, metric, instrumentation, policy);

        std::string solver_message;
        try
          {
            (void)solver.solve(nonfinite_control);
          }
        catch (const std::exception &exception)
          {
            solver_message = exception.what();
          }
        require(!solver_message.empty(),
                "native Problem B optimization accepted non-finite input");
        require(solver_message.find("non-finite") != std::string::npos,
                "native Problem B optimization lost the non-finite diagnostic");
        require(instrumentation.state_solve_calls == 0 &&
                  instrumentation.adjoint_solve_calls == 0 &&
                  instrumentation.value_evaluations == 0 &&
                  instrumentation.metric_inverse_apply_calls == 0 &&
                  instrumentation.solve_failures == 0,
                "native Problem B non-finite input violated failure counters");

        write_native_problem_b_failure_summary(
          artifact_root,
          "nonfinite_input",
          "metric: " + metric_message + "; solver: " + solver_message,
          instrumentation.state_solve_calls,
          instrumentation.metric_inverse_apply_calls);
        evidence.fail(solver_message);
        const auto read = [&artifact_root](const char *const filename) {
          std::ifstream input(artifact_root / filename);
          require(static_cast<bool>(input),
                  std::string("native Problem B non-finite probe lost ") +
                    filename);
          return std::string{std::istreambuf_iterator<char>(input),
                             std::istreambuf_iterator<char>()};
        };
        require(read("status.txt").find("status failed") != std::string::npos,
                "native Problem B non-finite probe did not retain failed status");
        require(read("failure.txt").find(solver_message) != std::string::npos,
                "native Problem B non-finite probe lost its diagnostic");
        require(read("failure-summary.txt").find("metric: ") !=
                  std::string::npos &&
                  read("failure-summary.txt").find("solver: ") !=
                    std::string::npos,
                "native Problem B non-finite probe lost metric diagnostics");
        require(read("counters.csv").find("native,state_solve_calls,0") !=
                  std::string::npos &&
                  read("solve-records.csv").find(
                    "path,status,role,iterations,initial_residual,final_residual") !=
                    std::string::npos,
                "native Problem B non-finite probe lost counters or solve evidence");
      }
    catch (...)
      {
        evidence.fail_current_exception();
        throw;
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
        {"problem_b_coordinates",
         "nmopt.external_tutorial_step_4.native_problem_b_coordinates",
         {"dealii", "application", "external", "tutorial", "native",
          "problem_b"},
         60,
         run_problem_b_coordinates_contract},
        {"problem_b_mass",
         "nmopt.external_tutorial_step_4.native_problem_b_mass",
         {"dealii", "application", "external", "tutorial", "native",
          "problem_b"},
         60,
         run_problem_b_mass_contract},
        {"problem_b_operations",
         "nmopt.external_tutorial_step_4.native_problem_b_operations",
         {"dealii", "application", "external", "tutorial", "native",
          "problem_b"},
         120,
         run_problem_b_operations_contract},
        {"native_problem_b_reduced",
         "nmopt.external_tutorial_step_4.native_problem_b_reduced",
         {"dealii", "application", "external", "tutorial", "native",
          "problem_b", "verification"},
         120,
         run_native_problem_b_reduced_contract},
        {"native_problem_b_oracle",
         "nmopt.external_tutorial_step_4.native_problem_b_oracle",
         {"dealii", "application", "external", "tutorial", "native",
          "problem_b", "verification"},
         120,
         run_native_problem_b_oracle_contract},
        {"native_problem_b_derivative_metric",
         "nmopt.external_tutorial_step_4.native_problem_b_derivative_metric",
         {"dealii", "application", "external", "tutorial", "native",
          "problem_b", "verification"},
         240,
         run_native_problem_b_derivative_metric_contract},
        {"native_problem_b_optimization",
         "nmopt.external_tutorial_step_4.native_problem_b_optimization",
         {"dealii", "application", "external", "tutorial", "native",
          "problem_b", "optimization"},
         180,
         run_native_problem_b_optimization_contract},
        {"native_problem_b_optimization_evidence",
         "nmopt.external_tutorial_step_4.native_problem_b_optimization_evidence",
         {"dealii", "application", "external", "tutorial", "native",
          "problem_b", "optimization", "verification"},
         900,
         run_native_problem_b_optimization_evidence},
        {"native_problem_b_line_search_failure",
         "nmopt.external_tutorial_step_4.native_problem_b_line_search_failure",
         {"dealii", "application", "external", "tutorial", "native",
          "problem_b", "optimization", "diagnostics"},
         180,
         run_native_problem_b_line_search_failure_contract},
        {"native_problem_b_nonfinite",
         "nmopt.external_tutorial_step_4.native_problem_b_nonfinite",
         {"dealii", "application", "external", "tutorial", "native",
          "problem_b", "optimization", "diagnostics"},
         180,
         run_native_problem_b_nonfinite_contract},
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
