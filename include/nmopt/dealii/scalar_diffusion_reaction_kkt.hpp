#pragma once

#include "nmopt/contract/quadratic_kkt.hpp"
#include "nmopt/dealii/scalar_diffusion_reaction.hpp"
#include "nmopt/dealii/serial_kkt_solver.hpp"

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::dealii_backend
{
  /**
   * The named block provenance for the serial scalar DTO KKT realization.
   *
   * The assembled matrix has row and column ordering
   * [state, control, multiplier], with the row roles
   * [state_stationarity, control_stationarity, equality].  Keeping this
   * record beside the deal.II block matrix makes the lowerer's assembly
   * explainable without making the generic KKT contract depend on deal.II.
   */
  struct ScalarDiffusionReactionKKTBlockProvenance
  {
    unsigned int row_block = 0;
    unsigned int column_block = 0;
    std::string row_role;
    std::string column_role;
    std::string source;
    double      factor = 1.0;
    bool        transposed = false;
  };

  /**
   * Serial deal.II realization of the scalar DTO equality-constrained
   * quadratic KKT product.
   *
   * Its matrix-free product is the backend-neutral contract oracle.  The
   * BlockSparseMatrix is an explicit lowerer representation of the same
   * operator and is intentionally exposed for assembled/matrix-free checks
   * and future serial Krylov policies.
   */
  template <int dim>
  class ScalarDiffusionReactionKKT final
  {
  public:
    using Model = ScalarDiffusionReactionModel<dim>;
    using Product = contract::EqualityConstrainedQuadraticKKTProductT<
      SerialBackend>;
    using SolveResult = contract::QuadraticKKTSolveResultT<SerialBackend>;
    using Vector = dealii::Vector<double>;
    using BlockVector = dealii::BlockVector<double>;

    explicit ScalarDiffusionReactionKKT(
      std::shared_ptr<const Model> model)
      : model_(std::move(model))
      , product_(make_product(model_))
    {
      assemble_matrix();
      assemble_rhs();
    }

    const Product &
    product() const
    {
      return product_;
    }

    SolveResult
    solve(const contract::QuadraticKKTSolverPolicy &policy) const
    {
      return solve_serial_quadratic_kkt(product_, policy);
    }

    const dealii::BlockSparseMatrix<double> &
    matrix() const
    {
      return matrix_;
    }

    const BlockVector &
    right_hand_side() const
    {
      return right_hand_side_;
    }

    const std::vector<std::string> &
    block_names() const
    {
      return block_names_;
    }

    const std::vector<ScalarDiffusionReactionKKTBlockProvenance> &
    block_provenance() const
    {
      return block_provenance_;
    }

  private:
    static constexpr unsigned int n_blocks = 3;

    static contract::LayoutPtr
    make_two_block_layout(const std::string &label,
                          const std::size_t state_dimension,
                          const std::size_t control_dimension)
    {
      return std::make_shared<const contract::BlockLayout>(
        label,
        std::vector<contract::SpaceId>{{"state"}, {"control"}},
        std::vector<std::size_t>{state_dimension, control_dimension});
    }

    static contract::LayoutPtr
    make_single_block_layout(const std::string &label,
                             const std::string &space,
                             const std::size_t dimension)
    {
      return std::make_shared<const contract::BlockLayout>(
        label,
        std::vector<contract::SpaceId>{{space}},
        std::vector<std::size_t>{dimension});
    }

    static Product
    make_product(const std::shared_ptr<const Model> &model)
    {
      contract::require(static_cast<bool>(model),
                        "Scalar DTO KKT realization needs a model");

      const std::size_t state_dimension =
        static_cast<std::size_t>(model->system_matrix_.m());
      const std::size_t control_dimension =
        static_cast<std::size_t>(model->control_mass_->m());
      const contract::LayoutPtr primal_layout = make_two_block_layout(
        "scalar_dto_kkt_primal", state_dimension, control_dimension);
      const contract::LayoutPtr multiplier_layout = make_single_block_layout(
        "scalar_dto_kkt_multiplier", "multiplier", state_dimension);
      const contract::LayoutPtr adjoint_layout = make_single_block_layout(
        "scalar_dto_kkt_adjoint", "adjoint", state_dimension);
      const contract::LayoutPtr stationarity_layout = make_two_block_layout(
        "scalar_dto_kkt_stationarity", state_dimension, control_dimension);
      const contract::LayoutPtr equality_layout = make_single_block_layout(
        "scalar_dto_kkt_equality", "state_equation", state_dimension);
      const Product::Layout layout(primal_layout,
                                   multiplier_layout,
                                   adjoint_layout,
                                   stationarity_layout,
                                   equality_layout);

      const auto quadratic_action = [model, stationarity_layout](
                                      const Product::Primal &primal) {
        Vector state(model->system_matrix_.m());
        model->state_mass_.vmult(state, primal.block(0));
        Vector control(model->control_mass_->m());
        model->control_mass_->vmult(control, primal.block(1));
        control *= model->regularisation_weight_;
        return Product::Covector(stationarity_layout,
                                 {std::move(state), std::move(control)});
      };

      const auto equality_action = [model, equality_layout](
                                      const Product::Primal &primal) {
        Vector value(model->system_matrix_.m());
        model->system_matrix_.vmult(value, primal.block(0));
        Vector control_contribution(model->system_matrix_.m());
        model->control_coupling_.vmult(control_contribution, primal.block(1));
        value.add(-1.0, control_contribution);
        return Product::Covector(equality_layout, {std::move(value)});
      };

      const auto multiplier_action = [model, stationarity_layout](
                                       const Product::Primal &multiplier) {
        Vector state(model->system_matrix_.m());
        model->system_matrix_.Tvmult(state, multiplier.block(0));
        Vector control(model->control_mass_->m());
        model->control_coupling_.Tvmult(control, multiplier.block(0));
        control *= -1.0;
        return Product::Covector(stationarity_layout,
                                 {std::move(state), std::move(control)});
      };

      const auto transpose_action = [model,
                                     primal_layout,
                                     multiplier_layout](const Product::Seed &seed) {
        Vector state(model->system_matrix_.m());
        model->state_mass_.Tvmult(state, seed.stationarity.block(0));
        Vector equality_state(model->system_matrix_.m());
        model->system_matrix_.Tvmult(equality_state, seed.equality.block(0));
        state.add(1.0, equality_state);

        Vector control(model->control_mass_->m());
        model->control_mass_->Tvmult(control, seed.stationarity.block(1));
        control *= model->regularisation_weight_;
        Vector equality_control(model->control_mass_->m());
        model->control_coupling_.Tvmult(equality_control,
                                        seed.equality.block(0));
        control.add(-1.0, equality_control);

        Vector multiplier(model->system_matrix_.m());
        model->system_matrix_.vmult(multiplier, seed.stationarity.block(0));
        Vector control_contribution(model->system_matrix_.m());
        model->control_coupling_.vmult(control_contribution,
                                       seed.stationarity.block(1));
        multiplier.add(-1.0, control_contribution);

        return Product::TransposeResult{
          Product::Covector(primal_layout,
                            {std::move(state), std::move(control)}),
          Product::Covector(multiplier_layout, {std::move(multiplier)})};
      };

      Vector state_rhs(model->desired_state_load_);
      Vector control_rhs(control_dimension);
      Vector equality_rhs(model->forcing_load_);
      const Product::MultiplierConversion conversion{
        "KKT multiplier lambda equals negative framework adjoint",
        [adjoint_layout](const Product::Primal &multiplier) {
          Vector value = multiplier.block(0);
          value *= -1.0;
          return Product::Primal(adjoint_layout, {std::move(value)});
        },
        [multiplier_layout](const Product::Primal &adjoint) {
          Vector value = adjoint.block(0);
          value *= -1.0;
          return Product::Primal(multiplier_layout, {std::move(value)});
        }};
      const contract::QuadraticKKTAssumptions assumptions{
        true,
        true,
        "serial scalar diffusion-reaction equality Jacobian is full row rank",
        "quadratic objective is positive on the equality-Jacobian kernel"};

      return Product(
        layout,
        quadratic_action,
        equality_action,
        multiplier_action,
        transpose_action,
        Product::Covector(stationarity_layout,
                          {std::move(state_rhs), std::move(control_rhs)}),
        Product::Covector(equality_layout, {std::move(equality_rhs)}),
        conversion,
        assumptions,
        contract::QuadraticKKTSymmetry::symmetric_indefinite);
    }

    static void
    add_sparsity_entries(dealii::DynamicSparsityPattern &target,
                         const dealii::SparseMatrix<double> &source,
                         const bool transpose)
    {
      for (dealii::types::global_dof_index row = 0; row < source.m(); ++row)
        for (auto entry = source.begin(row); entry != source.end(row); ++entry)
          target.add(transpose ? entry->column() : row,
                    transpose ? row : entry->column());
    }

    static void
    add_matrix_entries(dealii::SparseMatrix<double> &target,
                       const dealii::SparseMatrix<double> &source,
                       const double factor,
                       const bool transpose)
    {
      for (dealii::types::global_dof_index row = 0; row < source.m(); ++row)
        for (auto entry = source.begin(row); entry != source.end(row); ++entry)
          target.add(transpose ? entry->column() : row,
                     transpose ? row : entry->column(),
                     factor * entry->value());
    }

    void
    assemble_matrix()
    {
      const auto state_dimension =
        static_cast<dealii::types::global_dof_index>(
          model_->system_matrix_.m());
      const auto control_dimension =
        static_cast<dealii::types::global_dof_index>(
          model_->control_mass_->m());
      const std::vector<dealii::types::global_dof_index> sizes{
        state_dimension, control_dimension, state_dimension};

      dealii::BlockDynamicSparsityPattern dynamic(n_blocks, n_blocks);
      for (unsigned int row = 0; row < n_blocks; ++row)
        for (unsigned int column = 0; column < n_blocks; ++column)
          dynamic.block(row, column).reinit(sizes[row], sizes[column]);

      add_sparsity_entries(dynamic.block(0, 0), model_->state_mass_, false);
      add_sparsity_entries(dynamic.block(0, 2), model_->system_matrix_, true);
      add_sparsity_entries(dynamic.block(1, 2),
                           model_->control_coupling_,
                           true);
      add_sparsity_entries(dynamic.block(1, 1), *model_->control_mass_, false);
      add_sparsity_entries(dynamic.block(2, 0), model_->system_matrix_, false);
      add_sparsity_entries(dynamic.block(2, 1),
                           model_->control_coupling_,
                           false);
      dynamic.collect_sizes();
      matrix_sparsity_.copy_from(dynamic);
      matrix_.reinit(matrix_sparsity_);

      add_matrix_entries(matrix_.block(0, 0), model_->state_mass_, 1.0, false);
      add_matrix_entries(matrix_.block(0, 2),
                         model_->system_matrix_,
                         1.0,
                         true);
      add_matrix_entries(matrix_.block(1, 1),
                         *model_->control_mass_,
                         model_->regularisation_weight_,
                         false);
      add_matrix_entries(matrix_.block(1, 2),
                         model_->control_coupling_,
                         -1.0,
                         true);
      add_matrix_entries(matrix_.block(2, 0),
                         model_->system_matrix_,
                         1.0,
                         false);
      add_matrix_entries(matrix_.block(2, 1),
                         model_->control_coupling_,
                         -1.0,
                         false);

      block_names_ = {"state", "control", "multiplier"};
      block_provenance_ = {
        {0, 0, "state_stationarity", "state",
         "objective state tracking Hessian: state_mass", 1.0, false},
        {0, 1, "state_stationarity", "control",
         "zero objective state-control Hessian block", 0.0, false},
        {0, 2, "state_stationarity", "multiplier",
         "constraint Jacobian transpose: system_matrix", 1.0, true},
        {1, 0, "control_stationarity", "state",
         "zero objective control-state Hessian block", 0.0, false},
        {1, 1, "control_stationarity", "control",
         "objective control regularisation: regularisation_weight * control_mass",
         model_->regularisation_weight_, false},
        {1, 2, "control_stationarity", "multiplier",
         "constraint Jacobian transpose: -control_coupling", -1.0, true},
        {2, 0, "equality", "state",
         "constraint Jacobian: system_matrix", 1.0, false},
        {2, 1, "equality", "control",
         "constraint Jacobian: -control_coupling", -1.0, false},
        {2, 2, "equality", "multiplier", "zero multiplier block", 0.0, false}};
    }

    void
    assemble_rhs()
    {
      right_hand_side_.reinit({
        static_cast<dealii::types::global_dof_index>(
          model_->system_matrix_.m()),
        static_cast<dealii::types::global_dof_index>(
          model_->control_mass_->m()),
        static_cast<dealii::types::global_dof_index>(
          model_->system_matrix_.m())});
      right_hand_side_.block(0) = model_->desired_state_load_;
      right_hand_side_.block(1) = 0.0;
      right_hand_side_.block(2) = model_->forcing_load_;
    }

    std::shared_ptr<const Model> model_;
    Product product_;
    dealii::BlockSparsityPattern matrix_sparsity_;
    dealii::BlockSparseMatrix<double> matrix_;
    BlockVector right_hand_side_;
    std::vector<std::string> block_names_;
    std::vector<ScalarDiffusionReactionKKTBlockProvenance> block_provenance_;
  };
} // namespace nmopt::dealii_backend
