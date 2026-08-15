#pragma once

#include "nmopt/contract/quadratic_kkt.hpp"
#include "nmopt/contract/supplied_otd.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::reference
{
  using namespace nmopt::contract;

  struct LinearQuadraticKKTData
  {
    DenseMatrix A;
    DenseMatrix B;
    DenseVector f;
    DenseMatrix C;
    DenseVector desired_observation;
    DenseVector observation_weights;
    DenseVector regularisation_weights;
    double      alpha = 0.0;

    LinearQuadraticKKTData(DenseMatrix matrix_a,
                           DenseMatrix matrix_b,
                           DenseVector load,
                           DenseMatrix observation,
                           DenseVector desired,
                           DenseVector observation_weights_input,
                           DenseVector regularisation_weights_input,
                           const double alpha_input)
      : A(std::move(matrix_a))
      , B(std::move(matrix_b))
      , f(std::move(load))
      , C(std::move(observation))
      , desired_observation(std::move(desired))
      , observation_weights(std::move(observation_weights_input))
      , regularisation_weights(std::move(regularisation_weights_input))
      , alpha(alpha_input)
    {
      require(A.rows() == A.columns(),
              "Reference KKT state matrix must be square");
      require(B.rows() == A.rows(),
              "Reference KKT control coupling has the wrong row count");
      require(f.size() == A.rows(),
              "Reference KKT load has the wrong size");
      require(C.columns() == A.columns(),
              "Reference KKT observation has the wrong column count");
      require(desired_observation.size() == C.rows() &&
                observation_weights.size() == C.rows(),
              "Reference KKT observation data has the wrong size");
      require(regularisation_weights.size() == B.columns(),
              "Reference KKT regularisation data has the wrong size");
      require(alpha > 0.0, "Reference KKT alpha must be positive");
      for (std::size_t index = 0; index < observation_weights.size(); ++index)
        require(observation_weights[index] > 0.0,
                "Reference KKT observation weights must be positive");
      for (std::size_t index = 0; index < regularisation_weights.size();
           ++index)
        require(regularisation_weights[index] > 0.0,
                "Reference KKT regularisation weights must be positive");
    }
  };

  namespace detail
  {
    inline DenseVector
    observation_normal_action(const LinearQuadraticKKTData &data,
                              const DenseVector &             state)
    {
      DenseVector observation = data.C.vmult(state);
      for (std::size_t index = 0; index < observation.size(); ++index)
        observation[index] *= data.observation_weights[index];
      return data.C.transpose_vmult(observation);
    }

    inline DenseVector
    weighted_observation_target(const LinearQuadraticKKTData &data)
    {
      DenseVector target = data.desired_observation;
      for (std::size_t index = 0; index < target.size(); ++index)
        target[index] *= data.observation_weights[index];
      return data.C.transpose_vmult(target);
    }

    inline DenseVector
    regularisation_action(const LinearQuadraticKKTData &data,
                          const DenseVector &             control)
    {
      DenseVector value = control;
      for (std::size_t index = 0; index < value.size(); ++index)
        value[index] *= data.alpha * data.regularisation_weights[index];
      return value;
    }

    inline LayoutPtr
    make_two_block_layout(const std::string &label,
                          const SpaceId &   first_space,
                          const std::size_t first_dimension,
                          const SpaceId &   second_space,
                          const std::size_t second_dimension)
    {
      return std::make_shared<const BlockLayout>(
        label,
        std::vector<SpaceId>{first_space, second_space},
        std::vector<std::size_t>{first_dimension, second_dimension});
    }

    inline PrimalBlock
    make_supplied_point(const SuppliedOTDSystem &system,
                        const DenseVector &       state,
                        const DenseVector &       adjoint,
                        const DenseVector &       control)
    {
      const auto &selection = system.block_selection();
      std::vector<DenseVector> blocks;
      blocks.reserve(system.variable_layout()->n_blocks());
      for (std::size_t block = 0; block < system.variable_layout()->n_blocks();
           ++block)
        blocks.emplace_back(
          DenseVector(system.variable_layout()->dimension(block)));
      blocks.at(selection.state_variable) = state;
      blocks.at(selection.adjoint_variable) = adjoint;
      blocks.at(selection.control_variable) = control;
      return PrimalBlock(system.variable_layout(), std::move(blocks));
    }

    inline PrimalBlock
    make_supplied_seed(const SuppliedOTDSystem &system,
                       const DenseVector &       state_seed,
                       const DenseVector &       adjoint_seed,
                       const DenseVector &       control_seed)
    {
      const auto &selection = system.block_selection();
      std::vector<DenseVector> blocks;
      blocks.reserve(system.residual_layout()->n_blocks());
      for (std::size_t block = 0; block < system.residual_layout()->n_blocks();
           ++block)
        blocks.emplace_back(
          DenseVector(system.residual_layout()->dimension(block)));
      blocks.at(selection.state_equation) = state_seed;
      blocks.at(selection.adjoint_equation) = adjoint_seed;
      blocks.at(selection.control_stationarity) = control_seed;
      return PrimalBlock(system.residual_layout(), std::move(blocks));
    }

    inline LayoutPtr
    make_selected_variable_layout(const SuppliedOTDSystem &system,
                                  const std::string &       label)
    {
      const auto &selection = system.block_selection();
      const auto &source = system.variable_layout();
      return make_two_block_layout(
        label,
        source->space(selection.state_variable),
        source->dimension(selection.state_variable),
        source->space(selection.control_variable),
        source->dimension(selection.control_variable));
    }

    inline LayoutPtr
    make_selected_residual_layout(const SuppliedOTDSystem &system,
                                  const std::string &       label)
    {
      const auto &selection = system.block_selection();
      const auto &source = system.residual_layout();
      return make_two_block_layout(
        label,
        source->space(selection.adjoint_equation),
        source->dimension(selection.adjoint_equation),
        source->space(selection.control_stationarity),
        source->dimension(selection.control_stationarity));
    }
  } // namespace detail

  inline EqualityConstrainedQuadraticKKTProduct
  make_dto_kkt_product(const LinearQuadraticKKTData &input)
  {
    const auto data = std::make_shared<const LinearQuadraticKKTData>(input);
    const std::size_t state_dimension = data->A.columns();
    const std::size_t control_dimension = data->B.columns();

    const LayoutPtr primal_layout = detail::make_two_block_layout(
      "dto_kkt_primal",
      SpaceId{"state"},
      state_dimension,
      SpaceId{"control"},
      control_dimension);
    const LayoutPtr multiplier_layout = std::make_shared<const BlockLayout>(
      "dto_kkt_multiplier",
      std::vector<SpaceId>{{"multiplier"}},
      std::vector<std::size_t>{state_dimension});
    const LayoutPtr adjoint_layout = std::make_shared<const BlockLayout>(
      "dto_kkt_adjoint",
      std::vector<SpaceId>{{"adjoint"}},
      std::vector<std::size_t>{state_dimension});
    const LayoutPtr stationarity_layout = detail::make_two_block_layout(
      "dto_kkt_stationarity",
      SpaceId{"state_stationarity"},
      state_dimension,
      SpaceId{"control_stationarity"},
      control_dimension);
    const LayoutPtr equality_layout = std::make_shared<const BlockLayout>(
      "dto_kkt_equality",
      std::vector<SpaceId>{{"state_equation"}},
      std::vector<std::size_t>{state_dimension});
    const EqualityConstrainedQuadraticKKTProduct::Layout layout(
      primal_layout,
      multiplier_layout,
      adjoint_layout,
      stationarity_layout,
      equality_layout);

    const auto quadratic_action = [data, stationarity_layout](
                                    const EqualityConstrainedQuadraticKKTProduct::Primal &
                                      primal) {
      DenseVector state = detail::observation_normal_action(*data,
                                                            primal.block(0));
      DenseVector control = detail::regularisation_action(*data,
                                                          primal.block(1));
      return EqualityConstrainedQuadraticKKTProduct::Covector(
        stationarity_layout, {std::move(state), std::move(control)});
    };
    const auto equality_action = [data, equality_layout](
                                  const EqualityConstrainedQuadraticKKTProduct::Primal &
                                    primal) {
      DenseVector value = data->A.vmult(primal.block(0));
      value.add_scaled(-1.0, data->B.vmult(primal.block(1)));
      return EqualityConstrainedQuadraticKKTProduct::Covector(
        equality_layout, {std::move(value)});
    };
    const auto multiplier_action = [data, stationarity_layout](
                                    const EqualityConstrainedQuadraticKKTProduct::Primal &
                                      multiplier) {
      DenseVector state = data->A.transpose_vmult(multiplier.block(0));
      DenseVector control = data->B.transpose_vmult(multiplier.block(0));
      control.scale(-1.0);
      return EqualityConstrainedQuadraticKKTProduct::Covector(
        stationarity_layout, {std::move(state), std::move(control)});
    };
    const auto transpose_action = [data, primal_layout, multiplier_layout](
                                   const EqualityConstrainedQuadraticKKTProduct::Seed &
                                     seed) {
      DenseVector state = detail::observation_normal_action(
        *data, seed.stationarity.block(0));
      state.add_scaled(1.0,
                       data->A.transpose_vmult(seed.equality.block(0)));
      DenseVector control = detail::regularisation_action(
        *data, seed.stationarity.block(1));
      control.add_scaled(-1.0,
                         data->B.transpose_vmult(seed.equality.block(0)));
      DenseVector multiplier = data->A.vmult(seed.stationarity.block(0));
      multiplier.add_scaled(-1.0,
                            data->B.vmult(seed.stationarity.block(1)));
      return EqualityConstrainedQuadraticKKTProduct::TransposeResult{
        EqualityConstrainedQuadraticKKTProduct::Covector(
          primal_layout, {std::move(state), std::move(control)}),
        EqualityConstrainedQuadraticKKTProduct::Covector(
          multiplier_layout, {std::move(multiplier)})};
    };

    DenseVector state_rhs = detail::weighted_observation_target(*data);
    DenseVector control_rhs(control_dimension);
    const QuadraticKKTMultiplierConversion conversion{
      "reference DTO multiplier lambda equals negative framework adjoint",
      [adjoint_layout](const PrimalBlock &multiplier) {
        DenseVector value = multiplier.block(0);
        value.scale(-1.0);
        return PrimalBlock(adjoint_layout, {std::move(value)});
      },
      [multiplier_layout](const PrimalBlock &adjoint) {
        DenseVector value = adjoint.block(0);
        value.scale(-1.0);
        return PrimalBlock(multiplier_layout, {std::move(value)});
      }};
    const QuadraticKKTAssumptions assumptions{
      true,
      true,
      "reference A is declared full row rank",
      "reference quadratic objective is positive on ker(D)"};

    return EqualityConstrainedQuadraticKKTProduct(
      layout,
      quadratic_action,
      equality_action,
      multiplier_action,
      transpose_action,
      EqualityConstrainedQuadraticKKTProduct::Covector(
        stationarity_layout, {std::move(state_rhs), std::move(control_rhs)}),
      EqualityConstrainedQuadraticKKTProduct::Covector(
        equality_layout, {data->f}),
      conversion,
      assumptions,
      QuadraticKKTSymmetry::symmetric_indefinite);
  }

  inline EqualityConstrainedQuadraticKKTProduct
  make_canonical_supplied_otd_kkt_product(const SuppliedOTDSystem &system)
  {
    const auto supplied = std::make_shared<const SuppliedOTDSystem>(system);
    const auto &layout_selection = supplied->block_selection();
    const auto &variables = supplied->variable_layout();
    const auto &residuals = supplied->residual_layout();
    const LayoutPtr primal_layout = detail::make_selected_variable_layout(
      *supplied, "supplied_otd_kkt_primal");
    const LayoutPtr multiplier_layout = variables->single_block(
      layout_selection.adjoint_variable, "supplied_otd_kkt_multiplier");
    const LayoutPtr adjoint_layout = variables->single_block(
      layout_selection.adjoint_variable, "supplied_otd_kkt_adjoint");
    const LayoutPtr stationarity_layout = detail::make_selected_residual_layout(
      *supplied, "supplied_otd_kkt_stationarity");
    const LayoutPtr equality_layout = residuals->single_block(
      layout_selection.state_equation, "supplied_otd_kkt_equality");
    const EqualityConstrainedQuadraticKKTProduct::Layout layout(
      primal_layout,
      multiplier_layout,
      adjoint_layout,
      stationarity_layout,
      equality_layout);

    const auto zero_point = [supplied] {
      return PrimalBlock::zeros(supplied->variable_layout());
    };
    const auto quadratic_action = [supplied, stationarity_layout, zero_point](
                                    const EqualityConstrainedQuadraticKKTProduct::Primal &
                                      primal) {
      const auto &selection = supplied->block_selection();
      const PrimalBlock tangent = detail::make_supplied_point(
        *supplied,
        primal.block(0),
        DenseVector(supplied->variable_layout()->dimension(
          selection.adjoint_variable)),
        primal.block(1));
      const CovectorBlock linearised = supplied->residual_jvp(
        zero_point(), tangent);
      DenseVector state = linearised.block(selection.adjoint_equation);
      state.scale(-1.0);
      DenseVector control = linearised.block(selection.control_stationarity);
      return EqualityConstrainedQuadraticKKTProduct::Covector(
        stationarity_layout, {std::move(state), std::move(control)});
    };
    const auto equality_action = [supplied, equality_layout, zero_point](
                                  const EqualityConstrainedQuadraticKKTProduct::Primal &
                                    primal) {
      const auto &selection = supplied->block_selection();
      const PrimalBlock tangent = detail::make_supplied_point(
        *supplied,
        primal.block(0),
        DenseVector(supplied->variable_layout()->dimension(
          selection.adjoint_variable)),
        primal.block(1));
      const CovectorBlock linearised = supplied->residual_jvp(
        zero_point(), tangent);
      return EqualityConstrainedQuadraticKKTProduct::Covector(
        equality_layout,
        {linearised.block(selection.state_equation)});
    };
    const auto multiplier_action = [supplied, stationarity_layout, zero_point](
                                    const EqualityConstrainedQuadraticKKTProduct::Primal &
                                      multiplier) {
      const auto &selection = supplied->block_selection();
      const PrimalBlock tangent = detail::make_supplied_point(
        *supplied,
        DenseVector(supplied->variable_layout()->dimension(
          selection.state_variable)),
        multiplier.block(0),
        DenseVector(supplied->variable_layout()->dimension(
          selection.control_variable)));
      const CovectorBlock linearised = supplied->residual_jvp(
        zero_point(), tangent);
      DenseVector state = linearised.block(selection.adjoint_equation);
      DenseVector control = linearised.block(selection.control_stationarity);
      control.scale(-1.0);
      return EqualityConstrainedQuadraticKKTProduct::Covector(
        stationarity_layout, {std::move(state), std::move(control)});
    };
    const auto transpose_action = [supplied, primal_layout, multiplier_layout](
                                   const EqualityConstrainedQuadraticKKTProduct::Seed &
                                     seed) {
      const auto &selection = supplied->block_selection();
      DenseVector adjoint_seed = seed.stationarity.block(0);
      adjoint_seed.scale(-1.0);
      const PrimalBlock residual_seed = detail::make_supplied_seed(
        *supplied,
        seed.equality.block(0),
        adjoint_seed,
        seed.stationarity.block(1));
      const CovectorBlock transposed = supplied->residual_vjp(
        PrimalBlock::zeros(supplied->variable_layout()), residual_seed);
      DenseVector state = transposed.block(selection.state_variable);
      DenseVector control = transposed.block(selection.control_variable);
      DenseVector multiplier = transposed.block(selection.adjoint_variable);
      multiplier.scale(-1.0);
      return EqualityConstrainedQuadraticKKTProduct::TransposeResult{
        EqualityConstrainedQuadraticKKTProduct::Covector(
          primal_layout, {std::move(state), std::move(control)}),
        EqualityConstrainedQuadraticKKTProduct::Covector(
          multiplier_layout, {std::move(multiplier)})};
    };

    const CovectorBlock zero_residual = supplied->residual(zero_point());
    DenseVector state_rhs = zero_residual.block(
      layout_selection.adjoint_equation);
    DenseVector control_rhs = zero_residual.block(
      layout_selection.control_stationarity);
    DenseVector equality_rhs = zero_residual.block(
      layout_selection.state_equation);
    equality_rhs.scale(-1.0);
    const QuadraticKKTMultiplierConversion conversion{
      "canonical supplied OTD multiplier lambda equals negative framework adjoint",
      [adjoint_layout](const PrimalBlock &multiplier) {
        DenseVector value = multiplier.block(0);
        value.scale(-1.0);
        return PrimalBlock(adjoint_layout, {std::move(value)});
      },
      [multiplier_layout](const PrimalBlock &adjoint) {
        DenseVector value = adjoint.block(0);
        value.scale(-1.0);
        return PrimalBlock(multiplier_layout, {std::move(value)});
      }};
    const QuadraticKKTAssumptions assumptions{
      true,
      true,
      "canonical supplied OTD reference has full row rank",
      "canonical supplied OTD quadratic objective is positive on ker(D)"};

    return EqualityConstrainedQuadraticKKTProduct(
      layout,
      quadratic_action,
      equality_action,
      multiplier_action,
      transpose_action,
      EqualityConstrainedQuadraticKKTProduct::Covector(
        stationarity_layout, {std::move(state_rhs), std::move(control_rhs)}),
      EqualityConstrainedQuadraticKKTProduct::Covector(
        equality_layout, {std::move(equality_rhs)}),
      conversion,
      assumptions,
      QuadraticKKTSymmetry::symmetric_indefinite);
  }
} // namespace nmopt::reference
