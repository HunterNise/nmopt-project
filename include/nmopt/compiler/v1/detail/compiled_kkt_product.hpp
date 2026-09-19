#pragma once

// Private compiler implementation header.
// Do not include from application or reporting code.

#include "nmopt/contract/executable_model.hpp"
#include "nmopt/contract/quadratic_kkt.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::compiler::v1
{

  namespace compiled_kkt_detail
  {
    template <typename Backend>
    contract::PrimalBlockT<Backend>
    copy_primal(const contract::PrimalBlockT<Backend> &source,
                const contract::LayoutPtr              &layout)
    {
      std::vector<typename Backend::Vector> blocks;
      blocks.reserve(source.n_blocks());
      for (std::size_t block = 0; block < source.n_blocks(); ++block)
        blocks.push_back(source.block(block));
      return contract::PrimalBlockT<Backend>(layout, std::move(blocks));
    }

    template <typename Values>
    void
    negate(Values &value)
    {
      for (std::size_t block = 0; block < value.n_blocks(); ++block)
        value.scale_block(block, -1.0);
    }
  } // namespace compiled_kkt_detail

  template <typename Backend>
  std::shared_ptr<const contract::EqualityConstrainedQuadraticKKTProductT<Backend>>
  make_compiled_dto_kkt_product(
    std::shared_ptr<const contract::ExecutableModelT<Backend>> executable)
  {
    contract::require(static_cast<bool>(executable),
                      "A compiled DTO KKT product needs an executable model");
    using Product = contract::EqualityConstrainedQuadraticKKTProductT<Backend>;
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;

    const auto primal_layout = executable->variable_layout();
    const auto multiplier_layout = executable->test_layout();
    const auto adjoint_layout = executable->test_layout();
    const auto stationarity_layout = executable->variable_layout();
    const auto equality_layout = executable->test_layout();
    const typename Product::Layout layout(primal_layout,
                                          multiplier_layout,
                                          adjoint_layout,
                                          stationarity_layout,
                                          equality_layout,
                                          {"compiled_primal_stationarity",
                                           {0, 1},
                                           {0, 1},
                                           {"state_stationarity",
                                            "control_stationarity"}},
                                          {"compiled_multiplier_equality",
                                           {0},
                                           {0},
                                           {"state_equation"}});
    const Primal zero = Primal::zeros(primal_layout);

    const auto quadratic_action = [executable, zero](const Primal &primal) {
      return contract::subtract(executable->objective_derivative(primal),
                                executable->objective_derivative(zero));
    };
    const auto equality_action = [executable, zero](const Primal &primal) {
      return executable->residual_jvp(zero, primal);
    };
    const auto multiplier_action = [executable, zero](const Primal &multiplier) {
      return executable->residual_vjp(zero, multiplier);
    };
    const auto transpose_action = [executable,
                                   zero,
                                   quadratic_action,
                                   equality_action](const typename Product::Seed &seed) {
      const Primal stationarity_seed = seed.stationarity;
      Covector primal_action = quadratic_action(stationarity_seed);
      const Covector equality_pullback =
        executable->residual_vjp(zero, seed.equality);
      for (std::size_t block = 0; block < primal_action.n_blocks(); ++block)
        primal_action.add_scaled_block(block,
                                       1.0,
                                       equality_pullback.block(block));
      return typename Product::TransposeResult{
        primal_action,
        equality_action(stationarity_seed)};
    };

    Covector stationarity_rhs = executable->objective_derivative(zero);
    compiled_kkt_detail::negate(stationarity_rhs);
    Covector equality_rhs = executable->residual(zero);
    compiled_kkt_detail::negate(equality_rhs);
    const typename Product::MultiplierConversion conversion{
      "KKT multiplier equals negative framework adjoint",
      [adjoint_layout](const Primal &multiplier) {
        Primal result = multiplier;
        compiled_kkt_detail::negate(result);
        return compiled_kkt_detail::copy_primal(result, adjoint_layout);
      },
      [multiplier_layout](const Primal &adjoint) {
        Primal result = adjoint;
        compiled_kkt_detail::negate(result);
        return compiled_kkt_detail::copy_primal(result, multiplier_layout);
      }};
    const contract::QuadraticKKTAssumptions assumptions{
      true,
      true,
      "canonical scalar DTO equality Jacobian is full row rank",
      "canonical scalar DTO quadratic objective is positive on the equality-Jacobian kernel",
      true,
      true,
      "canonical scalar DTO D-transpose and KKT-transpose actions are declared exact under the listed block pairings"};

    return std::make_shared<const Product>(
      layout,
      quadratic_action,
      equality_action,
      multiplier_action,
      transpose_action,
      std::move(stationarity_rhs),
      std::move(equality_rhs),
      conversion,
      assumptions,
      contract::QuadraticKKTSymmetry::symmetric_indefinite);
  }
} // namespace nmopt::compiler::v1
