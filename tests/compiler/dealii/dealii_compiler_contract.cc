#include "dealii_compiler_support.hpp"

namespace
{
#include "scenarios/compiled_formulations.hpp"
#include "scenarios/dirichlet_lowering.hpp"
#include "scenarios/natural_boundary_lowering.hpp"
#include "scenarios/observations.hpp"
#include "scenarios/control_metrics.hpp"
#include "scenarios/diagnostics.hpp"
#include "scenarios/session_ownership.hpp"
} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<test_support::Scenario> scenarios{
        {"canonical_volume_control",
         "nmopt.dealii.canonical_volume_control",
         {"dealii", "compiler"},
         180,
         []() { run_canonical_volume_control_contract_test<2>(); }},
        {"compiled_dto_kkt",
         "nmopt.dealii.compiled_dto_kkt",
         {"dealii", "compiler", "kkt"},
         60,
         []() { run_compiled_dto_kkt_contract_test<2>(); }},
        {"compiled_pdas",
         "nmopt.dealii.compiled_pdas",
         {"dealii", "compiler", "pdas", "active-set"},
         120,
         []() { run_compiled_pdas_contract_test<2>(); }},
        {"fixed_dirichlet",
         "nmopt.dealii.fixed_dirichlet",
         {"dealii", "compiler"},
         60,
         []() { run_fixed_dirichlet_contract_test<2>(); }},
        {"dirichlet_control",
         "nmopt.dealii.dirichlet_control",
         {"dealii", "compiler"},
         60,
         []() { run_dirichlet_control_contract_test<2>(); }},
        {"l2_dirichlet_transposition",
         "nmopt.dealii.l2_dirichlet_transposition",
         {"dealii", "compiler"},
         60,
         []() { run_l2_dirichlet_transposition_lowering_test<2>(); }},
        {"partial_dirichlet_control",
         "nmopt.dealii.partial_dirichlet_control",
         {"dealii", "compiler"},
         60,
         []() { run_partial_dirichlet_control_contract_test<2>(); }},
        {"subdomain_observation",
         "nmopt.dealii.subdomain_observation",
         {"dealii", "compiler"},
         60,
         []() { run_subdomain_observation_contract_test<2>(); }},
        {"point_sensor",
         "nmopt.dealii.point_sensor",
         {"dealii", "compiler", "observation"},
         60,
         []() { run_point_sensor_contract_test<2>(); }},
        {"normal_flux",
         "nmopt.dealii.normal_flux",
         {"dealii", "compiler", "observation"},
         60,
         []() { run_normal_flux_contract_test<2>(); }},
        {"h1_state_observation",
         "nmopt.dealii.h1_state_observation",
         {"dealii", "compiler"},
         60,
         []() { run_h1_state_observation_contract_test<2>(); }},
        {"neumann_boundary",
         "nmopt.dealii.neumann_boundary",
         {"dealii", "compiler"},
         60,
         []() { run_neumann_boundary_contract_test<2>(); }},
        {"neumann_convection_subdomain",
         "nmopt.dealii.neumann_convection_subdomain",
         {"dealii", "compiler"},
         60,
         []() { run_neumann_convection_subdomain_contract_test<2>(); }},
        {"weighted_boundary_trace",
         "nmopt.dealii.weighted_boundary_trace",
         {"dealii", "compiler"},
         60,
         []() { run_weighted_boundary_trace_contract_test<2>(); }},
        {"h1_control",
         "nmopt.dealii.h1_control",
         {"dealii", "compiler"},
         60,
         []() { run_h1_control_regularisation_contract_test<2>(); }},
        {"continuous_control_components",
         "nmopt.dealii.continuous_control_components",
         {"dealii", "compiler", "metric"},
         30,
         []() { run_continuous_control_component_contract_test<2>(); }},
        {"volume_observation_assembly",
         "nmopt.dealii.volume_observation_assembly",
         {"dealii", "compiler", "observation"},
         30,
         []() { run_volume_observation_assembly_contract_test<2>(); }},
        {"continuous_neumann_control_lowering",
         "nmopt.dealii.continuous_neumann_control_lowering",
         {"dealii", "compiler", "metric"},
         60,
         []() { run_continuous_neumann_control_lowering_contract_test<2>(); }},
        {"l2_tracking_continuous_control",
         "nmopt.dealii.l2_tracking_continuous_control",
         {"dealii", "compiler", "metric"},
         60,
         []() { run_l2_tracking_continuous_control_contract_test<2>(); }},
        {"simplex_continuous_control",
         "nmopt.dealii.simplex_continuous_control",
         {"dealii", "compiler", "metric"},
         60,
         []() { run_simplex_continuous_control_contract_test<2>(); }},
        {"hminus1_compilation",
         "nmopt.dealii.hminus1_compilation",
         {"dealii", "compiler", "metric"},
         60,
         []() { run_hminus1_compilation_contract_test<2>(); }},
        {"coefficient_identification",
         "nmopt.dealii.coefficient_identification",
         {"dealii", "compiler"},
         60,
         []() { run_coefficient_identification_contract_test<2>(); }},
        {"pure_neumann",
         "nmopt.dealii.pure_neumann",
         {"dealii", "compiler"},
         60,
         []() { run_pure_neumann_contract_test<2>(); }},
        {"general_scalar_robin",
         "nmopt.dealii.general_scalar_robin",
         {"dealii", "compiler"},
         90,
         []() { run_general_scalar_robin_contract_test<2>(); }},
        {"compiler_diagnostics",
         "nmopt.dealii.compiler_diagnostics",
         {"dealii", "compiler", "contract"},
         60,
         []() { run_compiler_diagnostics_contract_test<2>(); }},
        {"compiler_session",
         "nmopt.dealii.compiler_session",
         {"dealii", "compiler", "contract"},
         60,
         []() { run_compiler_session_contract_test<2>(); }},
        {"supplied_otd_owned_session_lifetime",
         "nmopt.dealii.supplied_otd_owned_session_lifetime",
         {"dealii", "compiler", "formulation", "supplied-otd", "ownership"},
         60,
         []() { run_supplied_otd_owned_session_lifetime_contract_test<2>(); }},
      };
      const auto result = test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "deal.II diffusion DTO contract scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "deal.II diffusion DTO contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}


