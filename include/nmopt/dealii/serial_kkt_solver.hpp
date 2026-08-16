#pragma once

#include "nmopt/contract/quadratic_kkt_solver.hpp"
#include "nmopt/dealii/serial_backend.hpp"

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace nmopt::dealii_backend
{
  using SerialKKTVector = dealii::Vector<double>;

  using SerialQuadraticKKTProduct =
    contract::EqualityConstrainedQuadraticKKTProductT<SerialBackend>;

  /**
   * Flat deal.II solver view of a typed quadratic KKT product.
   *
   * The typed product remains the source of truth. This adapter only packs
   * its primal/multiplier blocks into the scalar vector expected by the
   * serial MINRES/GMRES implementations and unpacks their actions again.
   */
  class SerialQuadraticKKTOperator final
  {
  public:
    using Vector = dealii::Vector<double>;
    using Product = SerialQuadraticKKTProduct;

    explicit SerialQuadraticKKTOperator(const Product &product)
      : product_(product)
      , domain_dimension_(layout_dimension(product.layout().primal) +
                          layout_dimension(product.layout().multiplier))
      , range_dimension_(layout_dimension(product.layout().stationarity) +
                         layout_dimension(product.layout().equality))
    {
      contract::require(domain_dimension_ == range_dimension_,
                        "Serial KKT operator has incompatible domain and range dimensions");
      contract::require(
        domain_dimension_ <= SerialBackend::maximum_native_size(),
        "Serial KKT operator dimension exceeds the native vector range");
    }

    dealii::types::global_dof_index
    m() const
    {
      return static_cast<dealii::types::global_dof_index>(range_dimension_);
    }

    dealii::types::global_dof_index
    n() const
    {
      return static_cast<dealii::types::global_dof_index>(domain_dimension_);
    }

    void
    vmult(Vector &destination, const Vector &source) const
    {
      contract::require(source.size() == domain_dimension_,
                        "Serial KKT operator source has the wrong dimension");
      contract::require(destination.size() == range_dimension_,
                        "Serial KKT operator destination has the wrong dimension");
      const auto action = product_.apply_kkt(point(source));
      flatten(destination, action, product_);
    }

    void
    right_hand_side(Vector &value) const
    {
      value.reinit(SerialBackend::checked_native_size(range_dimension_));
      const auto zero_point = Product::Point{
        Product::Primal::zeros(product_.layout().primal),
        Product::Primal::zeros(product_.layout().multiplier)};
      const auto zero_residual = product_.residual(zero_point);
      flatten(value, zero_residual, product_, -1.0);
    }

    Product::Point
    point(const Vector &value) const
    {
      contract::require(value.size() == domain_dimension_,
                        "Serial KKT solution has the wrong dimension");
      std::size_t offset = 0;
      auto primal_blocks = read_blocks(value,
                                       product_.layout().primal,
                                       offset);
      auto multiplier_blocks = read_blocks(value,
                                            product_.layout().multiplier,
                                            offset);
      contract::require(offset == domain_dimension_,
                        "Serial KKT solution packing did not consume its input");
      return Product::Point{
        Product::Primal(product_.layout().primal, std::move(primal_blocks)),
        Product::Primal(product_.layout().multiplier,
                        std::move(multiplier_blocks))};
    }

  private:
    static std::size_t
    layout_dimension(const contract::LayoutPtr &layout)
    {
      std::size_t dimension = 0;
      for (std::size_t block = 0; block < layout->n_blocks(); ++block)
        dimension += layout->dimension(block);
      return dimension;
    }

    static std::vector<Vector>
    read_blocks(const Vector &                 source,
                const contract::LayoutPtr &   layout,
                std::size_t &                 offset)
    {
      std::vector<Vector> blocks;
      blocks.reserve(layout->n_blocks());
      for (std::size_t block = 0; block < layout->n_blocks(); ++block)
        {
          Vector value(SerialBackend::checked_native_size(
            layout->dimension(block)));
          for (std::size_t index = 0; index < layout->dimension(block); ++index)
            value[SerialBackend::checked_native_size(index)] =
              source[SerialBackend::checked_native_size(offset + index)];
          offset += layout->dimension(block);
          blocks.emplace_back(std::move(value));
        }
      return blocks;
    }

    static void
    append(Vector &value,
           std::size_t &offset,
           const Vector &block,
           const double factor = 1.0)
    {
      for (dealii::types::global_dof_index index = 0; index < block.size(); ++index)
        value[SerialBackend::checked_native_size(offset + index)] =
          factor * block[index];
      offset += static_cast<std::size_t>(block.size());
    }

    static void
    append_paired(Vector &value,
                  std::size_t &offset,
                  const contract::CovectorBlockT<SerialBackend> &blocks,
                  const contract::QuadraticKKTBlockPairing &pairing,
                  const double factor)
    {
      for (std::size_t index = 0; index < pairing.domain_blocks.size(); ++index)
        append(value,
               offset,
               blocks.block(pairing.range_blocks[index]),
               factor);
    }

    static void
    flatten(Vector &value,
            const Product::Residual &residual,
            const Product &         product,
            const double            factor)
    {
      std::size_t offset = 0;
      append_paired(value,
                    offset,
                    residual.stationarity,
                    product.layout().primal_stationarity_pairing,
                    factor);
      append_paired(value,
                    offset,
                    residual.equality,
                    product.layout().multiplier_equality_pairing,
                    factor);
      contract::require(offset == value.size(),
                        "Serial KKT action packing did not fill its output");
    }

    static void
    flatten(Vector &value,
            const Product::Residual &residual,
            const Product &         product)
    {
      flatten(value, residual, product, 1.0);
    }

    const Product &product_;
    const std::size_t domain_dimension_;
    const std::size_t range_dimension_;
  };

  inline double
  block_norm(const contract::CovectorBlockT<SerialBackend> &value)
  {
    double squared_norm = 0.0;
    for (std::size_t block = 0; block < value.n_blocks(); ++block)
      {
        const double norm = value.block(block).l2_norm();
        squared_norm += norm * norm;
      }
    return std::sqrt(squared_norm);
  }

  inline contract::QuadraticKKTSolveResultT<SerialBackend>
  solve_serial_quadratic_kkt(
    const SerialQuadraticKKTProduct &                 product,
    const contract::QuadraticKKTSolverPolicy &       policy)
  {
    contract::validate(product, policy);
    const SerialQuadraticKKTOperator operator_view(product);
    SerialKKTVector rhs(SerialBackend::checked_native_size(operator_view.n()));
    operator_view.right_hand_side(rhs);
    SerialKKTVector solution(rhs.size());
    solution = 0.0;

    contract::require(
      operator_view.n() <= std::numeric_limits<unsigned int>::max() / 10U,
      "Serial KKT solve dimension exceeds the iteration-policy range");
    contract::require(
      policy.maximum_iterations <= std::numeric_limits<unsigned int>::max(),
      "Serial KKT maximum iteration count exceeds the native range");
    const unsigned int maximum_iterations =
      policy.maximum_iterations == 0
        ? std::max(100U,
                   10U * static_cast<unsigned int>(operator_view.n()))
        : static_cast<unsigned int>(policy.maximum_iterations);
    const double requested_tolerance =
      std::max(policy.absolute_tolerance,
               policy.relative_tolerance * rhs.l2_norm());
    dealii::SolverControl control(maximum_iterations, requested_tolerance);
    bool converged = true;
    if (policy.method == contract::QuadraticKKTSolverMethod::minres)
      {
        dealii::SolverMinRes<SerialKKTVector> solver(control);
        try
          {
            solver.solve(operator_view,
                         solution,
                         rhs,
                         dealii::PreconditionIdentity());
          }
        catch (const dealii::SolverControl::NoConvergence &)
          {
            converged = false;
          }
      }
    else
      {
        dealii::SolverGMRES<SerialKKTVector>::AdditionalData additional_data(
          policy.gmres_maximum_basis);
        dealii::SolverGMRES<SerialKKTVector> solver(control, additional_data);
        try
          {
            solver.solve(operator_view,
                         solution,
                         rhs,
                         dealii::PreconditionIdentity());
          }
        catch (const dealii::SolverControl::NoConvergence &)
          {
            converged = false;
          }
      }

    const auto point = operator_view.point(solution);
    const auto residual = product.residual(point);
    const double stationarity_residual = block_norm(residual.stationarity);
    const double equality_residual = block_norm(residual.equality);
    const bool residuals_converged =
      stationarity_residual <= requested_tolerance &&
      equality_residual <= requested_tolerance;
    const auto method = policy.method == contract::QuadraticKKTSolverMethod::minres
                          ? "serial_minres"
                          : "serial_gmres";
    const contract::QuadraticKKTSolveReport report{
      {method,
       "identity",
       maximum_iterations,
       control.last_step(),
       policy.relative_tolerance,
       policy.absolute_tolerance,
       requested_tolerance,
       control.last_value(),
       converged ? contract::LinearSolveTermination::converged
                 : contract::LinearSolveTermination::failed},
      stationarity_residual,
      equality_residual,
      residuals_converged};
    return {point, report};
  }
} // namespace nmopt::dealii_backend
