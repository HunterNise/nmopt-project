#include "nmopt/dealii/scalar_diffusion_reaction_kkt.hpp"
#include "nmopt/contract/supplied_otd_kkt.hpp"
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
  using SuppliedSystem = KKT::Model::SuppliedSystem;
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

  void
  test_supplied_otd_bridge_matches_dto_product()
  {
    dealii::Triangulation<2> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
    triangulation.refine_global(1);
    const auto model = make_model(triangulation);
    const KKT dto_kkt(model);
    const SuppliedSystem supplied =
      KKT::Model::make_supplied_otd_system(model);
    const Product supplied_product =
      nmopt::contract::make_canonical_supplied_otd_kkt_product(supplied);

    const std::size_t state_dimension = dto_kkt.product().layout().primal->dimension(0);
    const std::size_t control_dimension = dto_kkt.product().layout().primal->dimension(1);
    Vector state(state_dimension);
    Vector control(control_dimension);
    Vector multiplier(state_dimension);
    for (std::size_t index = 0; index < state_dimension; ++index)
      {
        state[index] = 0.12 + 0.01 * static_cast<double>(index);
        multiplier[index] = -0.08 + 0.015 * static_cast<double>(index);
      }
    for (std::size_t index = 0; index < control_dimension; ++index)
      control[index] = -0.06 + 0.025 * static_cast<double>(index);

    const Product::Point dto_point =
      make_point(dto_kkt.product(), state, control, multiplier);
    const Product::Point supplied_point =
      make_point(supplied_product, state, control, multiplier);
    const auto dto_action = dto_kkt.product().apply_kkt(dto_point);
    const auto supplied_action = supplied_product.apply_kkt(supplied_point);
    require_vector_close(supplied_action.stationarity.block(0),
                         dto_action.stationarity.block(0),
                         1e-12,
                         "supplied-OTD KKT state action differs from DTO");
    require_vector_close(supplied_action.stationarity.block(1),
                         dto_action.stationarity.block(1),
                         1e-12,
                         "supplied-OTD KKT control action differs from DTO");
    require_vector_close(supplied_action.equality.block(0),
                         dto_action.equality.block(0),
                         1e-12,
                         "supplied-OTD KKT equality action differs from DTO");

    const auto dto_residual = dto_kkt.product().residual(dto_point);
    const auto supplied_residual = supplied_product.residual(supplied_point);
    require_vector_close(supplied_residual.stationarity.block(0),
                         dto_residual.stationarity.block(0),
                         1e-12,
                         "supplied-OTD KKT state residual differs from DTO");
    require_vector_close(supplied_residual.stationarity.block(1),
                         dto_residual.stationarity.block(1),
                         1e-12,
                         "supplied-OTD KKT control residual differs from DTO");
    require_vector_close(supplied_residual.equality.block(0),
                         dto_residual.equality.block(0),
                         1e-12,
                         "supplied-OTD KKT equality residual differs from DTO");

    Product::Seed dto_seed{
      Product::Primal(dto_kkt.product().layout().stationarity,
                      {Vector(state_dimension), Vector(control_dimension)}),
      Product::Primal(dto_kkt.product().layout().equality,
                      {Vector(state_dimension)})};
    for (std::size_t index = 0; index < state_dimension; ++index)
      {
        dto_seed.stationarity.block(0)[index] =
          0.04 - 0.006 * static_cast<double>(index);
        dto_seed.equality.block(0)[index] =
          -0.03 + 0.009 * static_cast<double>(index);
      }
    for (std::size_t index = 0; index < control_dimension; ++index)
      dto_seed.stationarity.block(1)[index] =
        0.05 + 0.008 * static_cast<double>(index);
    const Product::Seed supplied_seed{
      Product::Primal(supplied_product.layout().stationarity,
                      {dto_seed.stationarity.block(0),
                       dto_seed.stationarity.block(1)}),
      Product::Primal(supplied_product.layout().equality,
                      {dto_seed.equality.block(0)})};
    const auto dto_transpose =
      dto_kkt.product().apply_kkt_transpose(dto_seed);
    const auto supplied_transpose =
      supplied_product.apply_kkt_transpose(supplied_seed);
    require_vector_close(supplied_transpose.primal.block(0),
                         dto_transpose.primal.block(0),
                         1e-12,
                         "supplied-OTD transposed state action differs from DTO");
    require_vector_close(supplied_transpose.primal.block(1),
                         dto_transpose.primal.block(1),
                         1e-12,
                         "supplied-OTD transposed control action differs from DTO");
    require_vector_close(supplied_transpose.multiplier.block(0),
                         dto_transpose.multiplier.block(0),
                         1e-12,
                         "supplied-OTD transposed multiplier differs from DTO");

    const auto supplied_adjoint =
      supplied_product.multiplier_to_adjoint(supplied_point.multiplier);
    const auto dto_adjoint =
      dto_kkt.product().multiplier_to_adjoint(dto_point.multiplier);
    require_vector_close(supplied_adjoint.block(0),
                         dto_adjoint.block(0),
                         1e-12,
                         "supplied-OTD multiplier conversion differs from DTO");

    const auto solution = supplied.solve(
      SuppliedSystem::Primal::zeros(supplied.variable_layout()));
    nmopt::contract::require(solution.report.converged(),
                             "supplied-OTD bridge solve did not converge");
    const auto &selection = supplied.block_selection();
    Vector solved_state = solution.solution.block(selection.state_variable);
    Vector solved_adjoint = solution.solution.block(selection.adjoint_variable);
    Vector solved_control = solution.solution.block(selection.control_variable);
    Vector solved_multiplier = solved_adjoint;
    solved_multiplier *= -1.0;
    const Product::Point solved_point = make_point(
      supplied_product, solved_state, solved_control, solved_multiplier);
    const auto solved_residual = supplied_product.residual(solved_point);
    require_close(solved_residual.stationarity.block(0).l2_norm(),
                  0.0,
                  1e-10,
                  "supplied-OTD bridge state solution residual");
    require_close(solved_residual.stationarity.block(1).l2_norm(),
                  0.0,
                  1e-10,
                  "supplied-OTD bridge control solution residual");
    require_close(solved_residual.equality.block(0).l2_norm(),
                  0.0,
                  1e-10,
                  "supplied-OTD bridge equality solution residual");
  }

  void
  test_krylov_solvers_cross_dto_and_supplied_otd()
  {
    dealii::Triangulation<2> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
    triangulation.refine_global(1);
    const auto model = make_model(triangulation);
    const KKT dto_kkt(model);
    const SuppliedSystem supplied =
      KKT::Model::make_supplied_otd_system(model);
    const Product supplied_product =
      nmopt::contract::make_canonical_supplied_otd_kkt_product(supplied);

    nmopt::contract::QuadraticKKTSolverPolicy minres_policy;
    minres_policy.maximum_iterations = 500;
    minres_policy.relative_tolerance = 1e-11;
    minres_policy.absolute_tolerance = 1e-13;
    const auto dto_minres = dto_kkt.solve(minres_policy);
    const auto supplied_minres =
      nmopt::dealii_backend::solve_serial_quadratic_kkt(supplied_product,
                                                         minres_policy);
    nmopt::contract::require(dto_minres.report.converged(),
                             "serial DTO MINRES did not converge");
    nmopt::contract::require(supplied_minres.report.converged(),
                             "serial supplied-OTD MINRES did not converge");
    nmopt::contract::require(
      dto_minres.report.linear_solve.algorithm == "serial_minres" &&
        supplied_minres.report.linear_solve.algorithm == "serial_minres",
      "MINRES solve report did not identify its algorithm");

    nmopt::contract::QuadraticKKTSolverPolicy gmres_policy = minres_policy;
    gmres_policy.method =
      nmopt::contract::QuadraticKKTSolverMethod::gmres;
    gmres_policy.gmres_maximum_basis = 12;
    const auto dto_gmres =
      nmopt::dealii_backend::solve_serial_quadratic_kkt(dto_kkt.product(),
                                                         gmres_policy);
    const auto supplied_gmres =
      nmopt::dealii_backend::solve_serial_quadratic_kkt(supplied_product,
                                                         gmres_policy);
    const auto check_report =
      [](const Product &product, const auto &result, const char *description) {
        const auto residual = product.residual(result.solution);
        const double stationarity_residual =
          nmopt::dealii_backend::block_norm(residual.stationarity);
        const double equality_residual =
          nmopt::dealii_backend::block_norm(residual.equality);
        require_close(result.report.stationarity_residual,
                      stationarity_residual,
                      1e-12,
                      std::string(description) + " stationarity report");
        require_close(result.report.equality_residual,
                      equality_residual,
                      1e-12,
                      std::string(description) + " equality report");
        const bool residuals_converged =
          stationarity_residual <=
            result.report.linear_solve.requested_tolerance &&
          equality_residual <=
            result.report.linear_solve.requested_tolerance;
        nmopt::contract::require(
          result.report.residuals_converged == residuals_converged,
          std::string(description) + " residual convergence flag");
      };
    check_report(dto_kkt.product(), dto_minres, "DTO MINRES");
    check_report(supplied_product, supplied_minres, "supplied MINRES");
    check_report(dto_kkt.product(), dto_gmres, "DTO GMRES");
    check_report(supplied_product, supplied_gmres, "supplied GMRES");
    nmopt::contract::require(dto_gmres.report.converged(),
                             "serial DTO GMRES did not converge");
    nmopt::contract::require(supplied_gmres.report.converged(),
                             "serial supplied-OTD GMRES did not converge");
    nmopt::contract::require(
      dto_gmres.report.linear_solve.algorithm == "serial_gmres" &&
        supplied_gmres.report.linear_solve.algorithm == "serial_gmres",
      "GMRES solve report did not identify its algorithm");

    nmopt::contract::QuadraticKKTSolverPolicy limited_policy = minres_policy;
    limited_policy.maximum_iterations = 1;
    limited_policy.relative_tolerance = 1e-14;
    limited_policy.absolute_tolerance = 1e-16;
    const auto limited =
      nmopt::dealii_backend::solve_serial_quadratic_kkt(
        dto_kkt.product(), limited_policy);
    check_report(dto_kkt.product(), limited, "limited MINRES");
    nmopt::contract::require(
      !limited.report.linear_solve.converged() &&
        limited.report.linear_solve.termination ==
          nmopt::contract::LinearSolveTermination::failed &&
        limited.report.linear_solve.maximum_iterations == 1 &&
        limited.report.linear_solve.iterations <=
          limited.report.linear_solve.maximum_iterations &&
        !limited.report.converged(),
      "iteration-limited KKT solve did not report linear failure separately");

    require_vector_close(dto_minres.solution.primal.block(0),
                         supplied_minres.solution.primal.block(0),
                         1e-8,
                         "DTO and supplied-OTD MINRES state solutions differ");
    require_vector_close(dto_minres.solution.primal.block(1),
                         supplied_minres.solution.primal.block(1),
                         1e-8,
                         "DTO and supplied-OTD MINRES control solutions differ");
    require_vector_close(dto_minres.solution.multiplier.block(0),
                         supplied_minres.solution.multiplier.block(0),
                         1e-8,
                         "DTO and supplied-OTD MINRES multipliers differ");
    require_vector_close(dto_gmres.solution.primal.block(0),
                         supplied_gmres.solution.primal.block(0),
                         1e-8,
                         "DTO and supplied-OTD GMRES state solutions differ");
    require_vector_close(dto_gmres.solution.primal.block(1),
                         supplied_gmres.solution.primal.block(1),
                         1e-8,
                         "DTO and supplied-OTD GMRES control solutions differ");
    require_vector_close(dto_gmres.solution.multiplier.block(0),
                         supplied_gmres.solution.multiplier.block(0),
                         1e-8,
                         "DTO and supplied-OTD GMRES multipliers differ");
  }

  void
  test_serial_operator_uses_declared_pairings()
  {
    using ContractProduct =
      nmopt::contract::EqualityConstrainedQuadraticKKTProductT<Backend>;
    const auto primal = std::make_shared<const nmopt::contract::BlockLayout>(
      "pairing_primal",
      std::vector<nmopt::contract::SpaceId>{{"primal_0"}, {"primal_1"}},
      std::vector<std::size_t>{1, 1});
    const auto multiplier = std::make_shared<const nmopt::contract::BlockLayout>(
      "pairing_multiplier",
      std::vector<nmopt::contract::SpaceId>{{"multiplier"}},
      std::vector<std::size_t>{1});
    const auto adjoint = std::make_shared<const nmopt::contract::BlockLayout>(
      "pairing_adjoint",
      std::vector<nmopt::contract::SpaceId>{{"adjoint"}},
      std::vector<std::size_t>{1});
    const auto stationarity =
      std::make_shared<const nmopt::contract::BlockLayout>(
        "pairing_stationarity",
        std::vector<nmopt::contract::SpaceId>{{"stationarity_0"},
                                              {"stationarity_1"}},
        std::vector<std::size_t>{1, 1});
    const auto equality = std::make_shared<const nmopt::contract::BlockLayout>(
      "pairing_equality",
      std::vector<nmopt::contract::SpaceId>{{"equality"}},
      std::vector<std::size_t>{1});
    const ContractProduct::Layout layout(
      primal,
      multiplier,
      adjoint,
      stationarity,
      equality,
      {"permuted_primal_stationarity", {0, 1}, {1, 0}, {"p0_s1", "p1_s0"}},
      {"multiplier_equality", {0}, {0}, {"m_e"}});
    const auto quadratic_action = [stationarity](const ContractProduct::Primal &value) {
      return ContractProduct::Covector(
        stationarity,
        {value.block(0), value.block(1)});
    };
    const auto equality_action = [equality](const ContractProduct::Primal &value) {
      return ContractProduct::Covector(equality, {value.block(0)});
    };
    const auto multiplier_action =
      [stationarity](const ContractProduct::Primal &value) {
        return ContractProduct::Covector(
          stationarity,
          {value.block(0), value.block(0)});
      };
    const auto transpose_action = [primal, multiplier](
                                    const ContractProduct::Seed &seed) {
      return ContractProduct::TransposeResult{
        ContractProduct::Covector(
          primal,
          {seed.stationarity.block(0), seed.stationarity.block(1)}),
        ContractProduct::Covector(multiplier, {seed.equality.block(0)})};
    };
    const ContractProduct::MultiplierConversion conversion{
      "pairing test multiplier conversion",
      [adjoint](const ContractProduct::Primal &value) {
        return ContractProduct::Primal(adjoint, {value.block(0)});
      },
      [multiplier](const ContractProduct::Primal &value) {
        return ContractProduct::Primal(multiplier, {value.block(0)});
      }};
    const nmopt::contract::QuadraticKKTAssumptions assumptions{
      true,
      true,
      "pairing test rank",
      "pairing test kernel positivity",
      true,
      true,
      "pairing test transpose declarations"};
    const ContractProduct product(
      layout,
      quadratic_action,
      equality_action,
      multiplier_action,
      transpose_action,
      ContractProduct::Covector(stationarity, {Vector(1), Vector(1)}),
      ContractProduct::Covector(equality, {Vector(1)}),
      conversion,
      assumptions,
      nmopt::contract::QuadraticKKTSymmetry::symmetric_indefinite);

    nmopt::dealii_backend::SerialQuadraticKKTOperator operator_view(product);
    Vector source(3);
    source[0] = 1.0;
    source[1] = 2.0;
    source[2] = 3.0;
    Vector destination(3);
    operator_view.vmult(destination, source);
    require_close(destination[0], 2.0 + 3.0, 1e-12,
                  "serial KKT operator ignored the declared stationarity pairing");
    require_close(destination[1], 1.0 + 3.0, 1e-12,
                  "serial KKT operator used the physical stationarity order");
    require_close(destination[2], 1.0, 1e-12,
                  "serial KKT operator ignored the declared equality pairing");
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
         test_named_block_assembly_matches_product},
        {"scalar_supplied_otd_kkt_bridge",
         "nmopt.dealii.scalar_supplied_otd_kkt_bridge",
         {"dealii", "contract", "kkt", "supplied-otd"},
         60,
         test_supplied_otd_bridge_matches_dto_product},
        {"scalar_krylov_cross_path",
         "nmopt.dealii.scalar_krylov_cross_path",
         {"dealii", "contract", "kkt", "solver", "supplied-otd"},
         60,
         test_krylov_solvers_cross_dto_and_supplied_otd},
        {"serial_operator_pairings",
         "nmopt.dealii.serial_operator_pairings",
         {"dealii", "contract", "kkt", "solver", "pairing"},
         30,
         test_serial_operator_uses_declared_pairings}};
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
