#pragma once

#include "nmopt/contract/metric_constraint.hpp"
#include "nmopt/dealii/metric_solve_parameters.hpp"
#include "nmopt/dealii/serial_backend.hpp"
#include "nmopt/dealii/serial_spd_solver.hpp"

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::dealii_backend
{
  // The quotient H1/2 trace metric induced by minimum H1 extension:
  //
  //   G = A_BB - A_BI A_II^{-1} A_IB,
  //
  // where A is a volume H1 mass-plus-stiffness matrix and B is the ordered
  // trace DoF set. Applying G performs one interior extension solve. Applying
  // G^{-1} solves the full volume problem with a trace-supported covector.
  class TraceHhalfMetric final : public contract::MetricT<SerialBackend>
  {
  public:
    using Matrix = dealii::SparseMatrix<double>;
    using Primal = contract::PrimalBlockT<SerialBackend>;
    using Covector = contract::CovectorBlockT<SerialBackend>;
    using Vector = dealii::Vector<double>;

    TraceHhalfMetric(std::string                   id,
                     contract::LayoutPtr           layout,
                     std::shared_ptr<const Matrix> volume_h1_matrix,
                     std::vector<std::size_t>      trace_dofs,
                     MetricSolveParameters         solve_parameters = {})
      : id_(std::move(id))
      , layout_(std::move(layout))
      , volume_h1_matrix_(std::move(volume_h1_matrix))
      , trace_dofs_(std::move(trace_dofs))
      , solve_parameters_(solve_parameters)
    {
      contract::require(!id_.empty(),
                        "H1/2 trace metric identifier must not be empty");
      contract::require(static_cast<bool>(layout_),
                        "H1/2 trace metric needs a layout");
      contract::require(layout_->n_blocks() == 1,
                        "H1/2 trace metric supports exactly one coefficient block");
      contract::require(static_cast<bool>(volume_h1_matrix_),
                        "H1/2 trace metric needs a volume H1 matrix");
      contract::require(volume_h1_matrix_->m() == volume_h1_matrix_->n(),
                        "H1/2 trace metric volume matrix must be square");
      contract::require(trace_dofs_.size() == layout_->dimension(0),
                        "H1/2 trace DoF count does not match its layout");
      contract::require(solve_parameters_.maximum_iterations > 0,
                        "H1/2 trace metric maximum iterations must be positive");
      contract::require(solve_parameters_.relative_tolerance > 0.0,
                        "H1/2 trace metric relative tolerance must be positive");
      contract::require(solve_parameters_.absolute_tolerance > 0.0,
                        "H1/2 trace metric absolute tolerance must be positive");
      build_interior_matrix();
    }

    const std::string &
    id() const override
    {
      return id_;
    }

    const contract::LayoutPtr &
    layout() const override
    {
      return layout_;
    }

    const MetricSolveParameters &
    solve_parameters() const
    {
      return solve_parameters_;
    }

    Vector
    apply_vector(const Vector &trace) const
    {
      contract::require(trace.size() == trace_dofs_.size(),
                        "H1/2 trace metric received an incompatible primal vector");
      Vector extension(volume_h1_matrix_->m());
      for (std::size_t index = 0; index < trace_dofs_.size(); ++index)
        extension[trace_dofs_[index]] = trace[index];

      if (!interior_dofs_.empty())
        {
          Vector interaction(volume_h1_matrix_->m());
          volume_h1_matrix_->vmult(interaction, extension);
          Vector right_hand_side(interior_dofs_.size());
          for (std::size_t index = 0; index < interior_dofs_.size(); ++index)
            right_hand_side[index] = -interaction[interior_dofs_[index]];
          const Vector interior = solve(*interior_matrix_, right_hand_side,
                                        "minimum-extension");
          for (std::size_t index = 0; index < interior_dofs_.size(); ++index)
            extension[interior_dofs_[index]] = interior[index];
        }

      Vector volume_action(volume_h1_matrix_->m());
      volume_h1_matrix_->vmult(volume_action, extension);
      return restrict_to_trace(volume_action);
    }

    Vector
    inverse_apply_vector(const Vector &trace_covector) const
    {
      contract::require(trace_covector.size() == trace_dofs_.size(),
                        "H1/2 trace metric received an incompatible covector vector");
      Vector volume_covector(volume_h1_matrix_->m());
      for (std::size_t index = 0; index < trace_dofs_.size(); ++index)
        volume_covector[trace_dofs_[index]] = trace_covector[index];
      return restrict_to_trace(
        solve(*volume_h1_matrix_, volume_covector, "inverse"));
    }

    Covector
    apply(const Primal &primal) const override
    {
      contract::require(primal.layout()->compatible_with(*layout_),
                        "H1/2 trace metric primal has an incompatible layout");
      return Covector(layout_, {apply_vector(primal.block(0))});
    }

    Primal
    inverse_apply(const Covector &covector) const override
    {
      contract::require(covector.layout()->compatible_with(*layout_),
                        "H1/2 trace metric covector has an incompatible layout");
      return Primal(layout_, {inverse_apply_vector(covector.block(0))});
    }

  private:
    void
    build_interior_matrix()
    {
      const std::size_t invalid = std::numeric_limits<std::size_t>::max();
      std::vector<std::size_t> full_to_interior(volume_h1_matrix_->m(), invalid);
      std::vector<bool> is_trace(volume_h1_matrix_->m(), false);
      for (const auto dof : trace_dofs_)
        {
          contract::require(dof < volume_h1_matrix_->m(),
                            "H1/2 trace DoF lies outside the volume matrix");
          contract::require(!is_trace[dof],
                            "H1/2 trace DoF map contains a duplicate");
          is_trace[dof] = true;
        }
      for (std::size_t dof = 0; dof < is_trace.size(); ++dof)
        if (!is_trace[dof])
          {
            full_to_interior[dof] = interior_dofs_.size();
            interior_dofs_.push_back(dof);
          }

      auto sparsity = std::make_shared<dealii::SparsityPattern>();
      dealii::DynamicSparsityPattern dynamic(interior_dofs_.size(),
                                               interior_dofs_.size());
      for (std::size_t row = 0; row < interior_dofs_.size(); ++row)
        for (auto entry = volume_h1_matrix_->begin(interior_dofs_[row]);
             entry != volume_h1_matrix_->end(interior_dofs_[row]); ++entry)
          if (full_to_interior[entry->column()] != invalid)
            dynamic.add(row, full_to_interior[entry->column()]);
      sparsity->copy_from(dynamic);
      auto matrix = std::make_shared<Matrix>(*sparsity);
      for (std::size_t row = 0; row < interior_dofs_.size(); ++row)
        for (auto entry = volume_h1_matrix_->begin(interior_dofs_[row]);
             entry != volume_h1_matrix_->end(interior_dofs_[row]); ++entry)
          if (full_to_interior[entry->column()] != invalid)
            matrix->add(row,
                        full_to_interior[entry->column()],
                        entry->value());
      interior_sparsity_ = std::move(sparsity);
      interior_matrix_ = std::move(matrix);
    }

    Vector
    restrict_to_trace(const Vector &volume_vector) const
    {
      Vector trace(trace_dofs_.size());
      for (std::size_t index = 0; index < trace_dofs_.size(); ++index)
        trace[index] = volume_vector[trace_dofs_[index]];
      return trace;
    }

    Vector
    solve(const Matrix &matrix,
          const Vector &right_hand_side,
          const char *  operation) const
    {
      Vector solution(right_hand_side.size());
      if (right_hand_side.size() == 0)
        return solution;
      const SPDLinearSolvePolicy policy{solve_parameters_.maximum_iterations,
                                        solve_parameters_.relative_tolerance,
                                        solve_parameters_.absolute_tolerance};
      const auto report =
        solve_serial_spd(matrix, solution, right_hand_side, policy);
      contract::require(report.converged(),
                        std::string("H1/2 trace metric ") + operation +
                          " solve did not converge");
      return solution;
    }

    std::string                                   id_;
    contract::LayoutPtr                           layout_;
    std::shared_ptr<const Matrix>                 volume_h1_matrix_;
    std::vector<std::size_t>                      trace_dofs_;
    std::vector<std::size_t>                      interior_dofs_;
    std::shared_ptr<const dealii::SparsityPattern> interior_sparsity_;
    std::shared_ptr<const Matrix>                  interior_matrix_;
    MetricSolveParameters                         solve_parameters_;
  };
} // namespace nmopt::dealii_backend
