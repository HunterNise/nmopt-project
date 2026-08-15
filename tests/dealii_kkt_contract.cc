#include "nmopt/dealii/scalar_diffusion_reaction_kkt.hpp"
#include "test_support/scenario_dispatch.hpp"

#include <deal.II/base/function_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <cmath>
#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
  using Backend = nmopt::dealii_backend::SerialBackend;
  using KKT = nmopt::dealii_backend::ScalarDiffusionReactionKKT<2>;
  using Product = KKT::Product;
  using Vector = dealii::Vector<double>;
  using BlockVector = dealii::BlockVector<double>;

  void
  require_close(const double actual,
                const double expected,
                const double tolerance,
                const std::string &description)
  {
    nmopt::contract::require(std::abs(actual - expected) <= tolerance,
                             description);
  }

  void
  require_vector_close(const Vector &actual,
                       const Vector &expected,
                       const double tolerance,
                       const std::string &description)
  {
    nmopt::contract::require(actual.size() == expected.size(),
                             description + " has incompatible sizes");
    Vector difference = actual;
    difference.add(-1.0, expected);
    require_close(difference.l2_norm(), 0.0, tolerance, description);
  }

  std::shared_ptr<const KKT::Model>
  make_model(dealii::Triangulation<2> &triangulation)
  {
    dealii::Functions::ConstantFunction<2> forcing(1.0);
    dealii::Functions::ConstantFunction<2> desired_state(0.25);
    return std::make_shared<const KKT::Model>(triangulation,
                                              forcing,
                                              desired_state,
                                              1.3,
                                              0.2,
                                              0.7,
                                              1);
  }

  Product::Point
  make_point(const Product &product,
             const Vector &state,
             const Vector &control,
             const Vector &multiplier)
  {
    return Product::Point{
      Product::Primal(product.layout().primal, {state, control}),
      Product::Primal(product.layout().multiplier, {multiplier})};
  }

  BlockVector
  make_block_vector(const std::vector<dealii::types::global_dof_index> &sizes,
                    const Vector &state,
                    const Vector &control,
                    const Vector &multiplier)
  {
    BlockVector result(sizes);
    result.block(0) = state;
    result.block(1) = control;
    result.block(2) = multiplier;
    return result;
  }

  void
  test_named_block_assembly_matches_product()
  {
    dealii::Triangulation<2> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
    triangulation.refine_global(1);
    const auto model = make_model(triangulation);
    const KKT kkt(model);
    const Product &product = kkt.product();

    nmopt::contract::require(kkt.block_names() ==
                               std::vector<std::string>{"state",
                                                        "control",
                                                        "multiplier"},
                             "KKT block names do not expose the ordering");
    nmopt::contract::require(kkt.block_provenance().size() == 9,
                             "KKT provenance does not describe all blocks");
    const auto &provenance = kkt.block_provenance();
    nmopt::contract::require(
      provenance[0].source.find("state_mass") != std::string::npos,
      "KKT state Hessian provenance is missing its source matrix");
    nmopt::contract::require(provenance[2].transposed &&
                               provenance[2].source.find("system_matrix") !=
                                 std::string::npos,
                             "KKT state multiplier block is not identified as A^T");
    nmopt::contract::require(provenance[5].factor == -1.0 &&
                               provenance[5].transposed,
                             "KKT control multiplier block lost its sign or transpose");
    nmopt::contract::require(provenance[8].factor == 0.0,
                             "KKT multiplier block is not explicitly zero");

    const auto state_dimension = product.layout().primal->dimension(0);
    const auto control_dimension = product.layout().primal->dimension(1);
    const std::vector<dealii::types::global_dof_index> sizes{
      static_cast<dealii::types::global_dof_index>(state_dimension),
      static_cast<dealii::types::global_dof_index>(control_dimension),
      static_cast<dealii::types::global_dof_index>(state_dimension)};
    nmopt::contract::require(kkt.matrix().n_block_rows() == 3 &&
                               kkt.matrix().n_block_cols() == 3,
                             "KKT matrix is not a three-by-three block matrix");
    for (unsigned int block = 0; block < 3; ++block)
      nmopt::contract::require(kkt.matrix().block(block, block).m() ==
                                 sizes[block],
                               "KKT diagonal block has the wrong size");

    Vector state(state_dimension);
    Vector control(control_dimension);
    Vector multiplier(state_dimension);
    for (std::size_t index = 0; index < state_dimension; ++index)
      {
        state[index] = 0.1 + 0.03 * static_cast<double>(index);
        multiplier[index] = -0.2 + 0.02 * static_cast<double>(index);
      }
    for (std::size_t index = 0; index < control_dimension; ++index)
      control[index] = -0.15 + 0.04 * static_cast<double>(index);

    const Product::Point point =
      make_point(product, state, control, multiplier);
    const auto product_action = product.apply_kkt(point);
    const BlockVector source =
      make_block_vector(sizes, state, control, multiplier);
    BlockVector assembled_action(sizes);
    kkt.matrix().vmult(assembled_action, source);
    require_vector_close(assembled_action.block(0),
                         product_action.stationarity.block(0),
                         1e-12,
                         "assembled state stationarity action");
    require_vector_close(assembled_action.block(1),
                         product_action.stationarity.block(1),
                         1e-12,
                         "assembled control stationarity action");
    require_vector_close(assembled_action.block(2),
                         product_action.equality.block(0),
                         1e-12,
                         "assembled equality action");

    Vector stationarity_state(state_dimension);
    Vector stationarity_control(control_dimension);
    Vector equality(state_dimension);
    for (std::size_t index = 0; index < state_dimension; ++index)
      {
        stationarity_state[index] = 0.07 - 0.01 * static_cast<double>(index);
        equality[index] = 0.11 + 0.015 * static_cast<double>(index);
      }
    for (std::size_t index = 0; index < control_dimension; ++index)
      stationarity_control[index] = 0.09 + 0.02 * static_cast<double>(index);
    const Product::Seed seed{
      Product::Primal(product.layout().stationarity,
                      {stationarity_state, stationarity_control}),
      Product::Primal(product.layout().equality, {equality})};
    const auto product_transpose = product.apply_kkt_transpose(seed);
    const BlockVector transpose_source = make_block_vector(
      sizes, stationarity_state, stationarity_control, equality);
    BlockVector assembled_transpose(sizes);
    kkt.matrix().Tvmult(assembled_transpose, transpose_source);
    require_vector_close(assembled_transpose.block(0),
                         product_transpose.primal.block(0),
                         1e-12,
                         "assembled transposed state action");
    require_vector_close(assembled_transpose.block(1),
                         product_transpose.primal.block(1),
                         1e-12,
                         "assembled transposed control action");
    require_vector_close(assembled_transpose.block(2),
                         product_transpose.multiplier.block(0),
                         1e-12,
                         "assembled transposed multiplier action");

    const Product::Point zero_point = Product::Point{
      Product::Primal::zeros(product.layout().primal),
      Product::Primal::zeros(product.layout().multiplier)};
    const auto zero_residual = product.residual(zero_point);
    for (std::size_t index = 0; index < state_dimension; ++index)
      {
        require_close(kkt.right_hand_side().block(0)[index],
                      -zero_residual.stationarity.block(0)[index],
                      1e-12,
                      "assembled state right-hand side");
        require_close(kkt.right_hand_side().block(2)[index],
                      -zero_residual.equality.block(0)[index],
                      1e-12,
                      "assembled equality right-hand side");
      }
    for (std::size_t index = 0; index < control_dimension; ++index)
      require_close(kkt.right_hand_side().block(1)[index],
                    -zero_residual.stationarity.block(1)[index],
                    1e-12,
                    "assembled control right-hand side");
  }
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<nmopt::test_support::Scenario> scenarios{
        {"scalar_dto_kkt_blocks",
         "nmopt.dealii.scalar_dto_kkt_blocks",
         {"dealii", "contract", "kkt", "assembled"},
         60,
         test_named_block_assembly_matches_product}};
      const auto result = nmopt::test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "deal.II scalar DTO KKT block contract scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "deal.II scalar DTO KKT block contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
