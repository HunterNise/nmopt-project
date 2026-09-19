# External deal.II Step-4 integration fixture. This module is included only
# after the root establishes deal.II setup, nmopt_dealii_contract, the nmopt
# target/test helpers, and scenario discovery.

# Pinned upstream, stripped baseline, and adapted tutorial targets
# This target intentionally compiles the pinned upstream tutorial source
# without linking the nmopt contract or adding nmopt includes.
add_executable(
  nmopt_external_tutorial_step_4
  apps/external-dealii/step-4/source/upstream/step-4.cc)
deal_ii_setup_target(nmopt_external_tutorial_step_4)

# This target compiles the mechanically stripped baseline without linking
# the nmopt contract or adding nmopt includes.
add_executable(
  nmopt_external_tutorial_step_4_stripped
  apps/external-dealii/step-4/source/baseline/step-4-stripped.cc)
deal_ii_setup_target(nmopt_external_tutorial_step_4_stripped)

# This target compiles the adapted tutorial seams without linking the nmopt
# contract or adding nmopt includes.
add_executable(
  nmopt_external_tutorial_step_4_adapted
  apps/external-dealii/step-4/source/adapted/step-4.cc)
deal_ii_setup_target(nmopt_external_tutorial_step_4_adapted)

# Minimal nmopt consumers
nmopt_add_dealii_executable(
  nmopt_external_step4_minimal_problem_a
  apps/external-dealii/step-4/minimal/problem_a.cc)

nmopt_add_dealii_executable(
  nmopt_external_step4_minimal_problem_b
  apps/external-dealii/step-4/minimal/problem_b.cc)

if(BUILD_TESTING)
  # Step-4 integration contract targets
  # This test includes only the adapted tutorial source and the standard
  # scenario-discovery helper. It intentionally has no nmopt library link.
  add_executable(
    nmopt_external_step4_native_contract_test
    tests/dealii/external_step4_native_contract.cc)
  # Keep the plain signature for compatibility with deal_ii_setup_target().
  target_link_libraries(nmopt_external_step4_native_contract_test
    nmopt_build_flags)
  deal_ii_setup_target(nmopt_external_step4_native_contract_test)
  scenario_discovery_register(nmopt_external_step4_native_contract_test)

  nmopt_add_dealii_test(
    nmopt_external_step4_minimal_problem_a_contract_test
    tests/application/external_step4_minimal_problem_a_contract.cc)
  add_dependencies(nmopt_external_step4_minimal_problem_a_contract_test
                   nmopt_external_step4_minimal_problem_a)

  nmopt_add_dealii_test(
    nmopt_external_step4_minimal_problem_b_contract_test
    tests/application/external_step4_minimal_problem_b_contract.cc)
  add_dependencies(nmopt_external_step4_minimal_problem_b_contract_test
                   nmopt_external_step4_minimal_problem_b)

  nmopt_add_dealii_test(
    nmopt_external_step4_nmopt_contract_test
    tests/application/external_step4_nmopt_contract.cc)

  nmopt_add_dealii_test(
    nmopt_external_step4_problem_b_nmopt_contract_test
    tests/application/external_step4_problem_b_nmopt_contract.cc)

  nmopt_add_dealii_test(
    nmopt_external_step4_problem_b_optimization_contract_test
    tests/application/external_step4_problem_b_optimization_contract.cc)

  nmopt_add_dealii_test(
    nmopt_external_step4_optimization_contract_test
    tests/application/external_step4_optimization_contract.cc)

  # Reproduction / forward-comparison registrations
  # Ordinary Python tooling contracts in the root remain optional. The Step-4
  # comparison tests require Python when deal.II testing is enabled.
  find_package(Python3 COMPONENTS Interpreter REQUIRED)

  set(nmopt_external_tutorial_forward_root
      "${CMAKE_SOURCE_DIR}/runs/external-dealii/step-4/forward-comparison")
  set(nmopt_external_tutorial_upstream_run_dir
      "${nmopt_external_tutorial_forward_root}/ctest/upstream")
  set(nmopt_external_tutorial_stripped_run_dir
      "${nmopt_external_tutorial_forward_root}/ctest/stripped")
  set(nmopt_external_tutorial_adapted_run_dir
      "${nmopt_external_tutorial_forward_root}/ctest/adapted")
  set(nmopt_external_tutorial_adapted_comparison_root
      "${nmopt_external_tutorial_forward_root}/adapted")
  file(MAKE_DIRECTORY
    "${nmopt_external_tutorial_upstream_run_dir}"
    "${nmopt_external_tutorial_stripped_run_dir}"
    "${nmopt_external_tutorial_adapted_run_dir}")
  add_test(
    NAME nmopt.external_tutorial_step_4.forward
    COMMAND nmopt_external_tutorial_step_4)
  set_tests_properties(
    nmopt.external_tutorial_step_4.forward
    PROPERTIES
      LABELS "dealii;application;external;tutorial;reproduction"
      TIMEOUT 30
      WORKING_DIRECTORY "${nmopt_external_tutorial_upstream_run_dir}")

  add_test(
    NAME nmopt.external_tutorial_step_4.stripped_forward
    COMMAND nmopt_external_tutorial_step_4_stripped)
  set_tests_properties(
    nmopt.external_tutorial_step_4.stripped_forward
    PROPERTIES
      LABELS "dealii;application;external;tutorial;baseline;reproduction"
      TIMEOUT 30
      WORKING_DIRECTORY "${nmopt_external_tutorial_stripped_run_dir}")

  add_test(
    NAME nmopt.external_tutorial_step_4.adapted_forward
    COMMAND nmopt_external_tutorial_step_4_adapted)
  set_tests_properties(
    nmopt.external_tutorial_step_4.adapted_forward
    PROPERTIES
      LABELS "dealii;application;external;tutorial;adapted;reproduction"
      TIMEOUT 30
      WORKING_DIRECTORY "${nmopt_external_tutorial_adapted_run_dir}")

  add_test(
    NAME nmopt.external_tutorial_step_4.forward_comparison
    COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_SOURCE_DIR}/tools/external_dealii/check_forward.py"
      --upstream-executable $<TARGET_FILE:nmopt_external_tutorial_step_4>
      --stripped-executable $<TARGET_FILE:nmopt_external_tutorial_step_4_stripped>
      --output-root "${nmopt_external_tutorial_forward_root}"
      --file solution-2d.vtk
      --file solution-3d.vtk)
  set_tests_properties(
    nmopt.external_tutorial_step_4.forward_comparison
    PROPERTIES
      LABELS "dealii;application;external;tutorial;baseline;reproduction"
      TIMEOUT 60
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

  add_test(
    NAME nmopt.external_tutorial_step_4.adapted_forward_comparison
    COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_SOURCE_DIR}/tools/external_dealii/check_forward.py"
      --upstream-executable $<TARGET_FILE:nmopt_external_tutorial_step_4>
      --stripped-executable $<TARGET_FILE:nmopt_external_tutorial_step_4_adapted>
      --stripped-label adapted
      --output-root "${nmopt_external_tutorial_adapted_comparison_root}"
      --file solution-2d.vtk
      --file solution-3d.vtk)
  set_tests_properties(
    nmopt.external_tutorial_step_4.adapted_forward_comparison
    PROPERTIES
      LABELS "dealii;application;external;tutorial;adapted;reproduction"
      TIMEOUT 60
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

  add_test(
    NAME nmopt.external_tutorial_step_4.forward_comparator_contract
    COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_SOURCE_DIR}/tests/tools/external_dealii_forward_contract.py")
  set_tests_properties(
    nmopt.external_tutorial_step_4.forward_comparator_contract
    PROPERTIES
      LABELS "dealii;application;external;tutorial;tool"
      TIMEOUT 30
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()
