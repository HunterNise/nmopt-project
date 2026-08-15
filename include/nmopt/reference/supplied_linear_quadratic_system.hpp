#pragma once

#include "nmopt/contract/supplied_otd.hpp"
#include "nmopt/reference/quadratic_kkt.hpp"

#include <memory>
#include <utility>
#include <vector>

namespace nmopt::reference
{
  using namespace nmopt::contract;

  // Reference supplied-OTD target:
  //
  //   A y - f - B u = 0,
  //   A^T p - C^T W(C y-d) = 0,
  //   B^T p + alpha R u = 0.
  //
  // The three equations are deliberately registered as independent supplied
  // actions.  This target shares data with LinearQuadraticModel for comparison
  // tests, but it does not obtain its blocks by differentiating that DTO
  // model.
  class SuppliedLinearQuadraticSystem final
  {
  public:
    SuppliedLinearQuadraticSystem(DenseMatrix A,
                                  DenseMatrix B,
                                  DenseVector f,
                                  DenseMatrix C,
                                  DenseVector desired_observation,
                                  DenseVector observation_weights,
                                  DenseVector regularisation_weights,
                                  const double alpha)
      : data_(std::make_shared<const Data>(std::move(A),
                                           std::move(B),
                                           std::move(f),
                                           std::move(C),
                                           std::move(desired_observation),
                                           std::move(observation_weights),
                                           std::move(regularisation_weights),
                                           alpha))
      , system_(make_system(data_))
    {}

    const SuppliedOTDSystem &
    system() const
    {
      return system_;
    }

    LinearQuadraticKKTData
    kkt_data() const
    {
      return LinearQuadraticKKTData(data_->A,
                                    data_->B,
                                    data_->f,
                                    data_->C,
                                    data_->desired_observation,
                                    data_->observation_weights,
                                    data_->regularisation_weights,
                                    data_->alpha);
    }

  private:
    struct Data
    {
      Data(DenseMatrix matrix_a,
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
        , variable_layout(std::make_shared<const BlockLayout>(
            "supplied_otd_variables",
            std::vector<SpaceId>{{"state"}, {"state_test"}, {"control"}},
            std::vector<std::size_t>{this->A.rows(),
                                     this->A.columns(),
                                     this->B.columns()}))
        , residual_layout(std::make_shared<const BlockLayout>(
            "supplied_otd_residuals",
            std::vector<SpaceId>{{"state_equation"},
                                 {"adjoint_equation"},
                                 {"control_stationarity"}},
            std::vector<std::size_t>{this->A.rows(),
                                     this->A.columns(),
                                     this->B.columns()}))
      {
        require(A.rows() == A.columns(),
                "Supplied reference state matrix must be square");
        require(B.rows() == A.rows(),
                "Supplied reference control coupling has the wrong row count");
        require(f.size() == A.rows(),
                "Supplied reference load has the wrong size");
        require(C.columns() == A.columns(),
                "Supplied reference observation has the wrong column count");
        require(desired_observation.size() == C.rows() &&
                  observation_weights.size() == C.rows(),
                "Supplied reference observation data has the wrong size");
        require(regularisation_weights.size() == B.columns(),
                "Supplied reference regularisation has the wrong size");
        require(this->alpha > 0.0,
                "Supplied reference alpha must be positive");
        for (std::size_t index = 0; index < observation_weights.size(); ++index)
          require(observation_weights[index] > 0.0,
                  "Supplied reference observation weights must be positive");
        for (std::size_t index = 0; index < regularisation_weights.size();
             ++index)
          require(regularisation_weights[index] > 0.0,
                  "Supplied reference regularisation weights must be positive");
      }

      DenseMatrix A;
      DenseMatrix B;
      DenseVector f;
      DenseMatrix C;
      DenseVector desired_observation;
      DenseVector observation_weights;
      DenseVector regularisation_weights;
      double      alpha;
      LayoutPtr   variable_layout;
      LayoutPtr   residual_layout;
    };

    static DenseVector
    weighted_observation(const Data &data, const DenseVector &state)
    {
      DenseVector observation = data.C.vmult(state);
      observation.add_scaled(-1.0, data.desired_observation);
      for (std::size_t index = 0; index < observation.size(); ++index)
        observation[index] *= data.observation_weights[index];
      return observation;
    }

    static DenseVector
    observation_normal_action(const Data &data, const DenseVector &state)
    {
      DenseVector observation = data.C.vmult(state);
      for (std::size_t index = 0; index < observation.size(); ++index)
        observation[index] *= data.observation_weights[index];
      return data.C.transpose_vmult(observation);
    }

    static DenseVector
    weighted_observation_target(const Data &data)
    {
      DenseVector target = data.desired_observation;
      for (std::size_t index = 0; index < target.size(); ++index)
        target[index] *= data.observation_weights[index];
      return data.C.transpose_vmult(target);
    }

    static DenseVector
    regularisation_action(const Data &data, const DenseVector &control)
    {
      DenseVector value = control;
      for (std::size_t index = 0; index < value.size(); ++index)
        value[index] *= data.alpha * data.regularisation_weights[index];
      return value;
    }

    static DenseMatrix
    observation_normal_matrix(const Data &data)
    {
      const std::size_t dimension = data.A.columns();
      std::vector<double> entries(dimension * dimension, 0.0);
      for (std::size_t row = 0; row < dimension; ++row)
        for (std::size_t column = 0; column < dimension; ++column)
          for (std::size_t observation = 0; observation < data.C.rows();
               ++observation)
            entries[row * dimension + column] +=
              data.C(observation, row) * data.observation_weights[observation] *
              data.C(observation, column);
      return DenseMatrix(dimension, dimension, std::move(entries));
    }

    static SuppliedOTDSystem
    make_system(const std::shared_ptr<const Data> &data)
    {
      const SuppliedOTDLayout layout(data->variable_layout,
                                     data->residual_layout);

      const auto residual = [data](const PrimalBlock &point) {
        DenseVector state = data->A.vmult(point.block(0));
        state.add_scaled(-1.0, data->f);
        state.add_scaled(-1.0, data->B.vmult(point.block(2)));

        DenseVector adjoint = data->A.transpose_vmult(point.block(1));
        adjoint.add_scaled(
          -1.0, data->C.transpose_vmult(weighted_observation(*data,
                                                              point.block(0))));

        DenseVector stationarity = data->B.transpose_vmult(point.block(1));
        stationarity.add_scaled(
          1.0, regularisation_action(*data, point.block(2)));

        return CovectorBlock(data->residual_layout,
                             {std::move(state),
                              std::move(adjoint),
                              std::move(stationarity)});
      };

      const auto residual_jvp = [data](const PrimalBlock &,
                                       const PrimalBlock &tangent) {
        DenseVector state = data->A.vmult(tangent.block(0));
        state.add_scaled(-1.0, data->B.vmult(tangent.block(2)));

        DenseVector adjoint = data->A.transpose_vmult(tangent.block(1));
        adjoint.add_scaled(
          -1.0, observation_normal_action(*data, tangent.block(0)));

        DenseVector stationarity = data->B.transpose_vmult(tangent.block(1));
        stationarity.add_scaled(
          1.0, regularisation_action(*data, tangent.block(2)));

        return CovectorBlock(data->residual_layout,
                             {std::move(state),
                              std::move(adjoint),
                              std::move(stationarity)});
      };

      const auto residual_vjp = [data](const PrimalBlock &,
                                       const PrimalBlock &seed) {
        DenseVector state = data->A.transpose_vmult(seed.block(0));
        state.add_scaled(
          -1.0, observation_normal_action(*data, seed.block(1)));

        DenseVector adjoint = data->A.vmult(seed.block(1));
        adjoint.add_scaled(1.0, data->B.vmult(seed.block(2)));

        DenseVector control = data->B.transpose_vmult(seed.block(0));
        control.scale(-1.0);
        control.add_scaled(1.0, regularisation_action(*data, seed.block(2)));

        return CovectorBlock(data->variable_layout,
                             {std::move(state),
                              std::move(adjoint),
                              std::move(control)});
      };

      const auto solve = [data](const PrimalBlock &) {
        const std::size_t state_dimension = data->A.rows();
        const std::size_t control_dimension = data->B.columns();
        const std::size_t system_dimension =
          2 * state_dimension + control_dimension;
        std::vector<double> entries(system_dimension * system_dimension, 0.0);
        const DenseMatrix normal = observation_normal_matrix(*data);

        for (std::size_t row = 0; row < state_dimension; ++row)
          {
            for (std::size_t column = 0; column < state_dimension; ++column)
              entries[row * system_dimension + column] = data->A(row, column);
            for (std::size_t column = 0; column < control_dimension; ++column)
              entries[row * system_dimension + 2 * state_dimension + column] =
                -data->B(row, column);
          }

        for (std::size_t row = 0; row < state_dimension; ++row)
          {
            for (std::size_t column = 0; column < state_dimension; ++column)
              {
                entries[(state_dimension + row) * system_dimension + column] =
                  -normal(row, column);
                entries[(state_dimension + row) * system_dimension +
                        state_dimension + column] = data->A(column, row);
              }
          }

        for (std::size_t row = 0; row < control_dimension; ++row)
          {
            for (std::size_t column = 0; column < state_dimension; ++column)
              entries[(2 * state_dimension + row) * system_dimension +
                      state_dimension + column] = data->B(column, row);
            entries[(2 * state_dimension + row) * system_dimension +
                    2 * state_dimension + row] =
              data->alpha * data->regularisation_weights[row];
          }

        DenseVector right_hand_side(system_dimension);
        for (std::size_t row = 0; row < state_dimension; ++row)
          right_hand_side[row] = data->f[row];
        const DenseVector target = weighted_observation_target(*data);
        for (std::size_t row = 0; row < state_dimension; ++row)
          right_hand_side[state_dimension + row] = -target[row];

        const DenseVector solution =
          DenseMatrix(system_dimension,
                      system_dimension,
                      std::move(entries))
            .solve(std::move(right_hand_side));
        return SuppliedOTDSystem::SolveResult(
          PrimalBlock(data->variable_layout,
                      {DenseVector(std::vector<double>(
                         solution.values().begin(),
                         solution.values().begin() + state_dimension)),
                       DenseVector(std::vector<double>(
                         solution.values().begin() + state_dimension,
                         solution.values().begin() + 2 * state_dimension)),
                       DenseVector(std::vector<double>(
                         solution.values().begin() + 2 * state_dimension,
                         solution.values().end()))}),
          LinearSolveReport{"dense_gaussian",
                            "not applicable",
                            1,
                            1,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            LinearSolveTermination::converged});
      };

      return SuppliedOTDSystem(
        layout, residual, residual_jvp, residual_vjp, solve);
    }

    std::shared_ptr<const Data> data_;
    SuppliedOTDSystem           system_;
  };
} // namespace nmopt::reference
