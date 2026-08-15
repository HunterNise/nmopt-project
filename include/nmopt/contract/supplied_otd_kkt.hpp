#pragma once

#include "nmopt/contract/quadratic_kkt.hpp"
#include "nmopt/contract/supplied_otd.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::contract
{
  namespace detail
  {
    template <typename Backend>
    PrimalBlockT<Backend>
    make_canonical_supplied_point(
      const SuppliedOTDSystemT<Backend> &system,
      const typename Backend::Vector &  state,
      const typename Backend::Vector &  adjoint,
      const typename Backend::Vector &  control)
    {
      const auto &selection = system.block_selection();
      std::vector<typename Backend::Vector> blocks;
      blocks.reserve(system.variable_layout()->n_blocks());
      for (std::size_t block = 0; block < system.variable_layout()->n_blocks();
           ++block)
        blocks.emplace_back(
          Backend::zeros(system.variable_layout()->dimension(block)));
      blocks.at(selection.state_variable) = state;
      blocks.at(selection.adjoint_variable) = adjoint;
      blocks.at(selection.control_variable) = control;
      return PrimalBlockT<Backend>(system.variable_layout(),
                                   std::move(blocks));
    }

    template <typename Backend>
    PrimalBlockT<Backend>
    make_canonical_supplied_seed(
      const SuppliedOTDSystemT<Backend> &system,
      const typename Backend::Vector &  state_seed,
      const typename Backend::Vector &  adjoint_seed,
      const typename Backend::Vector &  control_seed)
    {
      const auto &selection = system.block_selection();
      std::vector<typename Backend::Vector> blocks;
      blocks.reserve(system.residual_layout()->n_blocks());
      for (std::size_t block = 0; block < system.residual_layout()->n_blocks();
           ++block)
        blocks.emplace_back(
          Backend::zeros(system.residual_layout()->dimension(block)));
      blocks.at(selection.state_equation) = state_seed;
      blocks.at(selection.adjoint_equation) = adjoint_seed;
      blocks.at(selection.control_stationarity) = control_seed;
      return PrimalBlockT<Backend>(system.residual_layout(),
                                   std::move(blocks));
    }

    template <typename Backend>
    LayoutPtr
    make_canonical_supplied_primal_layout(
      const SuppliedOTDSystemT<Backend> &system,
      const std::string &               label)
    {
      const auto &selection = system.block_selection();
      const auto &source = system.variable_layout();
      return std::make_shared<const BlockLayout>(
        label,
        std::vector<SpaceId>{source->space(selection.state_variable),
                             source->space(selection.control_variable)},
        std::vector<std::size_t>{
          source->dimension(selection.state_variable),
          source->dimension(selection.control_variable)});
    }

    template <typename Backend>
    LayoutPtr
    make_canonical_supplied_stationarity_layout(
      const SuppliedOTDSystemT<Backend> &system,
      const std::string &               label)
    {
      const auto &selection = system.block_selection();
      const auto &source = system.residual_layout();
      return std::make_shared<const BlockLayout>(
        label,
        std::vector<SpaceId>{source->space(selection.adjoint_equation),
                             source->space(selection.control_stationarity)},
        std::vector<std::size_t>{
          source->dimension(selection.adjoint_equation),
          source->dimension(selection.control_stationarity)});
    }
  } // namespace detail

  /**
   * Lower the canonical supplied-OTD block system to the common quadratic
   * KKT boundary.
   *
   * The supplied system owns the state, adjoint, and control-stationarity
   * residual blocks.  This adapter only selects their canonical linearized
   * combinations: the negative adjoint equation supplies Q, the state
   * equation supplies D, and the supplied residual VJP supplies D^T.  It
   * therefore preserves supplied-OTD provenance without pretending that the
   * system was obtained by differentiating a DTO objective/residual pair.
   */
  template <typename Backend>
  EqualityConstrainedQuadraticKKTProductT<Backend>
  make_canonical_supplied_otd_kkt_product(
    const SuppliedOTDSystemT<Backend> &system)
  {
    using Product = EqualityConstrainedQuadraticKKTProductT<Backend>;
    using Vector = typename Backend::Vector;
    using Primal = PrimalBlockT<Backend>;
    using Covector = CovectorBlockT<Backend>;

    const auto supplied =
      std::make_shared<const SuppliedOTDSystemT<Backend>>(system);
    const auto &selection = supplied->block_selection();
    const auto &variables = supplied->variable_layout();
    const auto &residuals = supplied->residual_layout();
    const LayoutPtr primal_layout =
      detail::make_canonical_supplied_primal_layout(
        *supplied, "supplied_otd_kkt_primal");
    const LayoutPtr multiplier_layout = variables->single_block(
      selection.adjoint_variable, "supplied_otd_kkt_multiplier");
    const LayoutPtr adjoint_layout = variables->single_block(
      selection.adjoint_variable, "supplied_otd_kkt_adjoint");
    const LayoutPtr stationarity_layout =
      detail::make_canonical_supplied_stationarity_layout(
        *supplied, "supplied_otd_kkt_stationarity");
    const LayoutPtr equality_layout = residuals->single_block(
      selection.state_equation, "supplied_otd_kkt_equality");
    const typename Product::Layout layout(primal_layout,
                                          multiplier_layout,
                                          adjoint_layout,
                                          stationarity_layout,
                                          equality_layout);

    const auto zero_point = [supplied] {
      return Primal::zeros(supplied->variable_layout());
    };
    const auto quadratic_action = [supplied, stationarity_layout, zero_point](
                                    const typename Product::Primal &primal) {
      const auto &block_selection = supplied->block_selection();
      const Primal tangent = detail::make_canonical_supplied_point(
        *supplied,
        primal.block(0),
        Backend::zeros(supplied->variable_layout()->dimension(
          block_selection.adjoint_variable)),
        primal.block(1));
      const Covector linearised = supplied->residual_jvp(zero_point(), tangent);
      Vector state = linearised.block(block_selection.adjoint_equation);
      Backend::scale(state, -1.0);
      Vector control = linearised.block(block_selection.control_stationarity);
      return typename Product::Covector(stationarity_layout,
                                        {std::move(state), std::move(control)});
    };
    const auto equality_action = [supplied, equality_layout, zero_point](
                                  const typename Product::Primal &primal) {
      const auto &block_selection = supplied->block_selection();
      const Primal tangent = detail::make_canonical_supplied_point(
        *supplied,
        primal.block(0),
        Backend::zeros(supplied->variable_layout()->dimension(
          block_selection.adjoint_variable)),
        primal.block(1));
      const Covector linearised = supplied->residual_jvp(zero_point(), tangent);
      return typename Product::Covector(
        equality_layout,
        {linearised.block(block_selection.state_equation)});
    };
    const auto multiplier_action = [supplied, stationarity_layout, zero_point](
                                    const typename Product::Primal &multiplier) {
      const auto &block_selection = supplied->block_selection();
      const Primal tangent = detail::make_canonical_supplied_point(
        *supplied,
        Backend::zeros(supplied->variable_layout()->dimension(
          block_selection.state_variable)),
        multiplier.block(0),
        Backend::zeros(supplied->variable_layout()->dimension(
          block_selection.control_variable)));
      const Covector linearised = supplied->residual_jvp(zero_point(), tangent);
      Vector state = linearised.block(block_selection.adjoint_equation);
      Vector control = linearised.block(block_selection.control_stationarity);
      Backend::scale(control, -1.0);
      return typename Product::Covector(stationarity_layout,
                                        {std::move(state), std::move(control)});
    };
    const auto transpose_action = [supplied,
                                   primal_layout,
                                   multiplier_layout](const typename Product::Seed &seed) {
      const auto &block_selection = supplied->block_selection();
      Vector adjoint_seed = seed.stationarity.block(0);
      Backend::scale(adjoint_seed, -1.0);
      const Primal residual_seed = detail::make_canonical_supplied_seed(
        *supplied,
        seed.equality.block(0),
        adjoint_seed,
        seed.stationarity.block(1));
      const Covector transposed = supplied->residual_vjp(
        Primal::zeros(supplied->variable_layout()), residual_seed);
      Vector state = transposed.block(block_selection.state_variable);
      Vector control = transposed.block(block_selection.control_variable);
      Vector multiplier = transposed.block(block_selection.adjoint_variable);
      Backend::scale(multiplier, -1.0);
      return typename Product::TransposeResult{
        typename Product::Covector(primal_layout,
                                   {std::move(state), std::move(control)}),
        typename Product::Covector(multiplier_layout,
                                   {std::move(multiplier)})};
    };

    const Covector zero_residual = supplied->residual(zero_point());
    Vector state_rhs = zero_residual.block(selection.adjoint_equation);
    Vector control_rhs = zero_residual.block(selection.control_stationarity);
    Vector equality_rhs = zero_residual.block(selection.state_equation);
    Backend::scale(equality_rhs, -1.0);
    const typename Product::MultiplierConversion conversion{
      "canonical supplied OTD multiplier lambda equals negative framework adjoint",
      [adjoint_layout](const Primal &multiplier) {
        Vector value = multiplier.block(0);
        Backend::scale(value, -1.0);
        return Primal(adjoint_layout, {std::move(value)});
      },
      [multiplier_layout](const Primal &adjoint) {
        Vector value = adjoint.block(0);
        Backend::scale(value, -1.0);
        return Primal(multiplier_layout, {std::move(value)});
      }};
    const QuadraticKKTAssumptions assumptions{
      true,
      true,
      "canonical supplied OTD equality block declares full row rank",
      "canonical supplied OTD quadratic objective is positive on ker(D)"};

    return Product(
      layout,
      quadratic_action,
      equality_action,
      multiplier_action,
      transpose_action,
      typename Product::Covector(stationarity_layout,
                                 {std::move(state_rhs), std::move(control_rhs)}),
      typename Product::Covector(equality_layout, {std::move(equality_rhs)}),
      conversion,
      assumptions,
      QuadraticKKTSymmetry::symmetric_indefinite);
  }
} // namespace nmopt::contract
