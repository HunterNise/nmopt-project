#pragma once

#include "native_reduced.hpp"
#include "scenario.hpp"

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

namespace external_dealii_step4
{
  namespace verification
  {
    using Matrix  = ProblemA::Matrix;
    using Vector  = ProblemA::Vector;
    using Problem = ProblemA;

    inline double
    vector_difference(const Vector &left, const Vector &right)
    {
      if (left.size() != right.size())
        throw std::invalid_argument(
          "verification vectors have incompatible dimensions");
      Vector difference = left;
      difference.add(-1.0, right);
      return difference.l2_norm();
    }

    inline double
    scaled_vector_error(const Vector &left, const Vector &right)
    {
      return vector_difference(left, right) /
             std::max(1.0, std::max(left.l2_norm(), right.l2_norm()));
    }

    inline double
    scaled_scalar_error(const double left, const double right)
    {
      return std::abs(left - right) /
             std::max(1.0, std::max(std::abs(left), std::abs(right)));
    }

    inline void
    require_vector_close(const Vector &      actual,
                         const Vector &      expected,
                         const double        absolute_tolerance,
                         const double        relative_tolerance,
                         const std::string &message)
    {
      const double error = vector_difference(actual, expected);
      const double bound = absolute_tolerance +
                           relative_tolerance *
                             std::max(actual.l2_norm(), expected.l2_norm());
      if (error > bound)
        throw std::runtime_error(message + ": error=" +
                                 std::to_string(error) +
                                 ", bound=" + std::to_string(bound));
    }

    inline void
    require_scalar_close(const double        actual,
                         const double        expected,
                         const double        absolute_tolerance,
                         const double        relative_tolerance,
                         const std::string &message)
    {
      const double bound = absolute_tolerance +
                           relative_tolerance *
                             std::max(std::abs(actual), std::abs(expected));
      if (std::abs(actual - expected) > bound)
        throw std::runtime_error(message + ": actual=" +
                                 std::to_string(actual) +
                                 ", expected=" + std::to_string(expected) +
                                 ", bound=" + std::to_string(bound));
    }

    inline double
    matrix_symmetry_error(const Matrix &matrix)
    {
      if (matrix.m() != matrix.n())
        throw std::invalid_argument("verification matrix is not square");

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

    inline dealii::FullMatrix<double>
    dense_matrix(const Matrix &sparse)
    {
      dealii::FullMatrix<double> dense(sparse.m(), sparse.n());
      dense.copy_from(sparse);
      return dense;
    }

    struct OffSolutionPoint
    {
      Vector state;
      Vector control;
      Vector state_tangent;
      Vector control_tangent;
      Vector test_seed;
    };

    inline OffSolutionPoint
    make_off_solution_point()
    {
      const auto constant = scenario::constant_vector();
      const auto ramp     = scenario::ramp_vector();
      const auto alternating = scenario::alternating_vector();

      Vector state = scenario::scaled(constant, 0.2);
      state.add(0.3, ramp);
      Vector control = scenario::scaled(alternating, 0.2);
      control.add(0.1, constant);

      Vector seed = scenario::scaled(constant, 0.4);
      seed.add(0.2, alternating);

      return {std::move(state),
              std::move(control),
              scenario::normalized_ramp_direction(),
              scenario::normalized_alternating_direction(),
              std::move(seed)};
    }

    inline Vector
    centered_residual_jvp(const Problem &       problem,
                          const OffSolutionPoint &point,
                          const double           step)
    {
      Vector state_plus  = point.state;
      Vector state_minus = point.state;
      state_plus.add(step, point.state_tangent);
      state_minus.add(-step, point.state_tangent);

      Vector control_plus  = point.control;
      Vector control_minus = point.control;
      control_plus.add(step, point.control_tangent);
      control_minus.add(-step, point.control_tangent);

      auto residual_plus = problem.residual(state_plus, control_plus);
      auto residual_minus = problem.residual(state_minus, control_minus);
      residual_plus.add(-1.0, residual_minus);
      residual_plus *= 1.0 / (2.0 * step);
      return residual_plus;
    }

    inline double
    centered_objective_derivative(const Problem &       problem,
                                  const OffSolutionPoint &point,
                                  const double           step)
    {
      Vector state_plus  = point.state;
      Vector state_minus = point.state;
      state_plus.add(step, point.state_tangent);
      state_minus.add(-step, point.state_tangent);

      Vector control_plus  = point.control;
      Vector control_minus = point.control;
      control_plus.add(step, point.control_tangent);
      control_minus.add(-step, point.control_tangent);

      return (problem.objective(state_plus, control_plus) -
              problem.objective(state_minus, control_minus)) /
             (2.0 * step);
    }

    struct OracleResult
    {
      Vector state;
      Vector control;
      double system_residual;
      double stationarity_residual;
    };

    inline OracleResult
    optimum_oracle(const Problem &problem)
    {
      const auto dense = dense_matrix(problem.system_matrix());
      const auto size  = problem.state_dimension();

      dealii::FullMatrix<double> normal_system(size, size);
      for (std::size_t row = 0; row < size; ++row)
        for (std::size_t column = 0; column < size; ++column)
          {
            double value = row == column ? 1.0 : 0.0;
            for (std::size_t inner = 0; inner < size; ++inner)
              value += dense(inner, row) * dense(inner, column);
            normal_system(row, column) = value;
          }

      Vector normal_rhs(size);
      dense.Tvmult(normal_rhs, problem.system_rhs());

      dealii::LAPACKFullMatrix<double> factor(size, size);
      factor = normal_system;
      factor.compute_lu_factorization();
      Vector state = normal_rhs;
      factor.solve(state);

      Vector system_residual(size);
      normal_system.vmult(system_residual, state);
      system_residual.add(-1.0, normal_rhs);
      const double system_scale = std::max(1.0, normal_rhs.l2_norm());

      Vector control(size);
      dense.vmult(control, state);
      control.add(-1.0, problem.system_rhs());

      Vector transposed_control(size);
      dense.Tvmult(transposed_control, control);
      Vector stationarity = state;
      stationarity += transposed_control;
      const double stationarity_scale =
        std::max(1.0,
                 std::max(state.l2_norm(), transposed_control.l2_norm()));

      return {std::move(state),
              std::move(control),
              system_residual.l2_norm() / system_scale,
              stationarity.l2_norm() / stationarity_scale};
    }
  } // namespace verification
} // namespace external_dealii_step4
