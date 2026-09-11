#pragma once

#include "../integration/problem_b.hpp"

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>

namespace external_dealii_step4
{
  namespace problem_b_verification
  {
    using Vector = dealii::Vector<double>;

    struct EquationResidualAudit
    {
      double absolute_norm;
      double scale;
      double normalized;
    };

    struct DenseOperators
    {
      dealii::FullMatrix<double> K;
      dealii::FullMatrix<double> B;
      dealii::FullMatrix<double> M;
      dealii::FullMatrix<double> Q;
      Vector                     lifting;
      Vector                     free_rhs;
      Vector                     c;
      double                     alpha;
    };

    struct KKTResiduals
    {
      double state_stationarity;
      double control_stationarity;
      double feasibility;
    };

    struct DenseOracle
    {
      Vector        state;
      Vector        control;
      Vector        adjoint;
      KKTResiduals  residuals;
    };

    struct StateSolutionAudit
    {
      EquationResidualAudit full_equation;
      EquationResidualAudit reduced_equation;
      double                reconstruction_error;
      double                boundary_error;
    };

    struct AdjointSolutionAudit
    {
      EquationResidualAudit full_equation;
      double                homogeneous_reconstruction_error;
      double                boundary_error;
    };

    inline bool
    vector_is_finite(const Vector &vector)
    {
      for (unsigned int index = 0; index < vector.size(); ++index)
        if (!std::isfinite(vector[index]))
          return false;
      return true;
    }

    inline void
    require_finite(const double        value,
                   const std::string &message)
    {
      if (!std::isfinite(value))
        throw std::runtime_error(message + " is not finite");
    }

    inline void
    require_finite(const Vector &      vector,
                   const std::string &message)
    {
      if (!vector_is_finite(vector))
        throw std::runtime_error(message + " contains a non-finite value");
    }

    inline double
    vector_difference(const Vector &left, const Vector &right)
    {
      if (left.size() != right.size())
        throw std::invalid_argument(
          "Problem B verification vectors have incompatible dimensions");
      require_finite(left, "Problem B verification left vector");
      require_finite(right, "Problem B verification right vector");
      Vector difference = left;
      difference.add(-1.0, right);
      return difference.l2_norm();
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
      require_finite(error, message + " error");
      require_finite(bound, message + " bound");
      if (error > bound)
        throw std::runtime_error(message + ": error=" +
                                 std::to_string(error) +
                                 ", bound=" + std::to_string(bound));
    }

    inline EquationResidualAudit
    normalized_equation_residual(const Vector &lhs, const Vector &rhs)
    {
      if (lhs.size() != rhs.size())
        throw std::invalid_argument(
          "Problem B verification equation dimensions do not match");
      require_finite(lhs, "Problem B verification equation left value");
      require_finite(rhs, "Problem B verification equation right value");

      Vector residual = lhs;
      residual.add(-1.0, rhs);
      const double absolute_norm = residual.l2_norm();
      const double scale         = std::max(1.0, rhs.l2_norm());
      const double normalized    = absolute_norm / scale;
      require_finite(absolute_norm,
                     "Problem B verification equation absolute norm");
      require_finite(scale, "Problem B verification equation scale");
      require_finite(normalized,
                     "Problem B verification equation normalized value");
      return {absolute_norm, scale, normalized};
    }

    inline double
    normalized_block_residual(const Vector &                  residual,
                              const std::initializer_list<double> scales)
    {
      require_finite(residual, "Problem B KKT residual");
      double scale = 1.0;
      for (const double value : scales)
        {
          require_finite(value, "Problem B KKT residual scale");
          scale = std::max(scale, value);
        }
      return residual.l2_norm() / scale;
    }

    inline double
    matrix_value(const dealii::FullMatrix<double> &matrix,
                 const unsigned int                row,
                 const unsigned int                column)
    {
      return matrix(row, column);
    }

    template <typename Matrix>
    double
    matrix_value(const Matrix &      matrix,
                 const unsigned int row,
                 const unsigned int column)
    {
      return matrix.el(row, column);
    }

    template <typename Matrix>
    double
    relative_symmetry_error(const Matrix &matrix)
    {
      if (matrix.m() != matrix.n())
        throw std::invalid_argument(
          "Problem B verification symmetry matrix is not square");

      double squared_difference = 0.0;
      double squared_norm       = 0.0;
      for (unsigned int row = 0; row < matrix.m(); ++row)
        for (unsigned int column = 0; column < matrix.n(); ++column)
          {
            const double value = matrix_value(matrix, row, column);
            const double difference =
              value - matrix_value(matrix, column, row);
            require_finite(value, "Problem B verification matrix value");
            require_finite(difference,
                           "Problem B verification matrix difference");
            squared_norm += value * value;
            squared_difference += difference * difference;
          }
      return std::sqrt(squared_difference) /
             std::max(1.0, std::sqrt(squared_norm));
    }

    template <int dim, typename Application>
    DenseOperators
    make_dense_operators(const ProblemB<dim, Application> &problem,
                         const Application &                application)
    {
      const auto &coordinates = problem.coordinates();
      const auto &mass        = problem.mass();
      const auto &free_indices = coordinates.free_indices();
      const auto &system_matrix = application.system_matrix_view();

      const unsigned int state_size =
        static_cast<unsigned int>(problem.state_dimension());
      const unsigned int control_size =
        static_cast<unsigned int>(problem.control_dimension());
      if (free_indices.size() != state_size ||
          coordinates.lifting().size() != control_size ||
          system_matrix.m() != control_size ||
          system_matrix.n() != control_size)
        throw std::invalid_argument(
          "Problem B verification received incompatible operator dimensions");

      DenseOperators result{dealii::FullMatrix<double>(state_size, state_size),
                            dealii::FullMatrix<double>(state_size, control_size),
                            dealii::FullMatrix<double>(control_size, control_size),
                            dealii::FullMatrix<double>(state_size, state_size),
                            coordinates.lifting(),
                            problem.free_system_rhs(),
                            Vector(state_size),
                            problem.alpha()};

      for (unsigned int row = 0; row < control_size; ++row)
        for (unsigned int column = 0; column < control_size; ++column)
          result.M(row, column) = mass.mass_matrix().el(row, column);

      for (unsigned int row = 0; row < state_size; ++row)
        {
          const auto native_row = free_indices[row];
          for (unsigned int column = 0; column < state_size; ++column)
            {
              const auto native_column = free_indices[column];
              result.K(row, column) =
                system_matrix.el(native_row, native_column);
              result.Q(row, column) = result.M(static_cast<unsigned int>(native_row),
                                               static_cast<unsigned int>(native_column));
            }

          for (unsigned int column = 0; column < control_size; ++column)
            result.B(row, column) = mass.coupling_matrix().el(row, column);
        }

      for (unsigned int row = 0; row < state_size; ++row)
        {
          double value = 0.0;
          for (unsigned int column = 0; column < control_size; ++column)
            value += result.M(static_cast<unsigned int>(free_indices[row]), column) *
                     result.lifting[column];
          result.c[row] = value;
        }

      require_finite(result.lifting, "Problem B dense lifting");
      require_finite(result.free_rhs, "Problem B dense free RHS");
      require_finite(result.c, "Problem B dense objective lifting term");
      return result;
    }

    inline void
    require_dense_cholesky(const dealii::FullMatrix<double> &matrix,
                           const std::string &               name)
    {
      if (matrix.m() != matrix.n())
        throw std::invalid_argument(name + " is not square");

      const unsigned int size = matrix.m();
      dealii::FullMatrix<double> factor(size, size);
      for (unsigned int row = 0; row < size; ++row)
        {
          for (unsigned int column = 0; column <= row; ++column)
            {
              double value = matrix(row, column);
              for (unsigned int inner = 0; inner < column; ++inner)
                value -= factor(row, inner) * factor(column, inner);

              require_finite(value, name + " Cholesky pivot");
              if (row == column)
                {
                  if (!(value > 0.0))
                    throw std::runtime_error(name +
                                             " is not positive definite");
                  factor(row, column) = std::sqrt(value);
                }
              else
                factor(row, column) = value / factor(column, column);
            }
          for (unsigned int column = row + 1; column < size; ++column)
            factor(row, column) = 0.0;
        }
    }

    inline DenseOracle
    dense_kkt_oracle(const DenseOperators &operators)
    {
      const unsigned int state_size   = operators.K.m();
      const unsigned int control_size = operators.M.m();
      const unsigned int adjoint_offset = state_size + control_size;
      const unsigned int total_size     = adjoint_offset + state_size;

      dealii::FullMatrix<double> kkt(total_size, total_size);
      kkt = 0.0;
      Vector rhs(total_size);
      rhs = 0.0;

      for (unsigned int row = 0; row < state_size; ++row)
        {
          rhs[row] = -operators.c[row];
          rhs[adjoint_offset + row] = -operators.free_rhs[row];
          for (unsigned int column = 0; column < state_size; ++column)
            {
              kkt(row, column) = operators.Q(row, column);
              kkt(row, adjoint_offset + column) =
                -operators.K(column, row);
              kkt(adjoint_offset + row, column) =
                -operators.K(row, column);
            }
        }

      for (unsigned int row = 0; row < control_size; ++row)
        for (unsigned int column = 0; column < control_size; ++column)
          kkt(state_size + row, state_size + column) =
            operators.alpha * operators.M(row, column);

      for (unsigned int state_index = 0; state_index < state_size;
           ++state_index)
        for (unsigned int control_index = 0; control_index < control_size;
             ++control_index)
          {
            kkt(state_size + control_index,
                adjoint_offset + state_index) =
              operators.B(state_index, control_index);
            kkt(adjoint_offset + state_index,
                state_size + control_index) =
              operators.B(state_index, control_index);
          }

      dealii::LAPACKFullMatrix<double> factor(total_size, total_size);
      factor = kkt;
      factor.compute_lu_factorization();
      factor.solve(rhs);

      Vector state(state_size);
      Vector control(control_size);
      Vector adjoint(state_size);
      for (unsigned int index = 0; index < state_size; ++index)
        {
          state[index] = rhs[index];
          adjoint[index] = rhs[adjoint_offset + index];
        }
      for (unsigned int index = 0; index < control_size; ++index)
        control[index] = rhs[state_size + index];

      Vector state_block(state_size);
      operators.Q.vmult(state_block, state);
      Vector transposed_adjoint(state_size);
      operators.K.Tvmult(transposed_adjoint, adjoint);
      const double state_term_norm = state_block.l2_norm();
      const double adjoint_term_norm = transposed_adjoint.l2_norm();
      state_block.add(-1.0, transposed_adjoint);
      state_block.add(1.0, operators.c);
      const double lifting_term_norm = operators.c.l2_norm();

      Vector control_block(control_size);
      operators.M.vmult(control_block, control);
      control_block *= operators.alpha;
      Vector transposed_control(control_size);
      operators.B.Tvmult(transposed_control, adjoint);
      const double control_term_norm = control_block.l2_norm();
      const double control_adjoint_norm = transposed_control.l2_norm();
      control_block.add(1.0, transposed_control);

      Vector feasibility(state_size);
      operators.K.vmult(feasibility, state);
      const double state_operator_norm = feasibility.l2_norm();
      feasibility *= -1.0;
      operators.B.vmult(transposed_adjoint, control);
      const double control_operator_norm = transposed_adjoint.l2_norm();
      feasibility.add(1.0, transposed_adjoint);
      feasibility.add(1.0, operators.free_rhs);
      const double rhs_norm = operators.free_rhs.l2_norm();

      const KKTResiduals residuals{
        normalized_block_residual(state_block,
                                  {state_term_norm,
                                   adjoint_term_norm,
                                   lifting_term_norm}),
        normalized_block_residual(control_block,
                                  {control_term_norm,
                                   control_adjoint_norm}),
        normalized_block_residual(feasibility,
                                  {state_operator_norm,
                                   control_operator_norm,
                                   rhs_norm})};

      require_finite(state, "Problem B dense oracle state");
      require_finite(control, "Problem B dense oracle control");
      require_finite(adjoint, "Problem B dense oracle adjoint");
      return {std::move(state),
              std::move(control),
              std::move(adjoint),
              residuals};
    }

    template <int dim, typename Application>
    StateSolutionAudit
    audit_state_solution(const ProblemB<dim, Application> &problem,
                         const Application &                application,
                         const Vector &                     state,
                         const Vector &                     full_state,
                         const Vector &                     control)
    {
      if (state.size() != problem.state_dimension() ||
          full_state.size() != problem.control_dimension() ||
          control.size() != problem.control_dimension())
        throw std::invalid_argument(
          "Problem B state audit received incompatible vector dimensions");

      Vector rhs = application.system_rhs_view();
      rhs.add(1.0,
              problem.coordinates().embed_free(
                problem.mass().coupling_apply(control)));
      Vector lhs(rhs.size());
      application.system_matrix_view().vmult(lhs, full_state);

      const auto full_equation = normalized_equation_residual(lhs, rhs);
      Vector zero_state(problem.state_dimension());
      zero_state = 0.0;
      const auto reduced_equation = normalized_equation_residual(
        problem.residual(state, control), zero_state);
      const auto reconstructed = problem.coordinates().reconstruct(state);
      double boundary_error = 0.0;
      for (const auto &[index, value] : application.boundary_values_view())
        boundary_error = std::max(boundary_error,
                                  std::abs(full_state[index] - value));

      require_finite(boundary_error, "Problem B state boundary error");
      return {full_equation,
              reduced_equation,
              vector_difference(reconstructed, full_state),
              boundary_error};
    }

    template <int dim, typename Application>
    AdjointSolutionAudit
    audit_adjoint_solution(const ProblemB<dim, Application> &problem,
                           const Application &                application,
                           const Vector &                     adjoint,
                           const Vector &                     full_adjoint,
                           const Vector &                     state_derivative)
    {
      if (adjoint.size() != problem.state_dimension() ||
          full_adjoint.size() != problem.control_dimension() ||
          state_derivative.size() != problem.state_dimension())
        throw std::invalid_argument(
          "Problem B adjoint audit received incompatible vector dimensions");

      const Vector rhs = problem.coordinates().embed_free(state_derivative);
      Vector lhs(rhs.size());
      application.system_matrix_view().Tvmult(lhs, full_adjoint);
      const auto full_equation = normalized_equation_residual(lhs, rhs);
      const auto homogeneous =
        problem.coordinates().embed_free(adjoint);
      double boundary_error = 0.0;
      for (const auto &[index, value] : application.boundary_values_view())
        {
          (void)value;
          boundary_error = std::max(boundary_error,
                                    std::abs(full_adjoint[index]));
        }

      require_finite(boundary_error, "Problem B adjoint boundary error");
      return {full_equation,
              vector_difference(homogeneous, full_adjoint),
              boundary_error};
    }
  } // namespace problem_b_verification
} // namespace external_dealii_step4
