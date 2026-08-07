cmake_minimum_required(VERSION 3.20)

project(nmopt LANGUAGES CXX)

add_library(nmopt_contract INTERFACE)
target_include_directories(nmopt_contract INTERFACE include)
target_compile_features(nmopt_contract INTERFACE cxx_std_17)

include(CTest)

function(nmopt_add_scenario_test name target scenario labels)
  add_test(NAME ${name} COMMAND ${target} ${scenario})
  set_tests_properties(${name} PROPERTIES LABELS "${labels}")
endfunction()

if(BUILD_TESTING)
  add_executable(nmopt_contract_tests tests/reduced_dto_contract.cc)
  target_link_libraries(nmopt_contract_tests PRIVATE nmopt_contract)
  nmopt_add_scenario_test(
    nmopt.contract.v0
    nmopt_contract_tests
    v0_contract
    "backend-neutral;contract")
  nmopt_add_scenario_test(
    nmopt.contract.backend_parameterisation
    nmopt_contract_tests
    backend_parameterisation
    "backend-neutral;contract")

  add_executable(nmopt_semantic_v1_contract_test tests/semantic_v1_contract.cc)
  target_link_libraries(nmopt_semantic_v1_contract_test PRIVATE nmopt_contract)
  nmopt_add_scenario_test(
    nmopt.semantic.v1_validation
    nmopt_semantic_v1_contract_test
    validation
    "backend-neutral;semantic")
endif()

option(NMOPT_ENABLE_DEAL_II "Build the serial deal.II reference lowerer" ON)

if(NMOPT_ENABLE_DEAL_II)
  find_package(deal.II CONFIG QUIET)

  if(deal.II_FOUND)
    deal_ii_initialize_cached_variables()
    add_library(nmopt_dealii_contract INTERFACE)
    target_link_libraries(nmopt_dealii_contract INTERFACE nmopt_contract)
    set(NMOPT_USING_DEAL_II_CMAKE_PACKAGE TRUE)
  else()
    find_path(NMOPT_DEAL_II_INCLUDE_DIR deal.II/base/config.h)
    find_library(NMOPT_DEAL_II_LIBRARY NAMES deal.ii)

    if(NMOPT_DEAL_II_INCLUDE_DIR AND NMOPT_DEAL_II_LIBRARY)
      add_library(nmopt_dealii_contract INTERFACE)
      target_link_libraries(nmopt_dealii_contract
        INTERFACE
          nmopt_contract
          "${NMOPT_DEAL_II_LIBRARY}")
      target_include_directories(nmopt_dealii_contract
        INTERFACE
          "${NMOPT_DEAL_II_INCLUDE_DIR}")
    endif()
  endif()

  if(TARGET nmopt_dealii_contract AND BUILD_TESTING)
    add_executable(nmopt_dealii_diffusion_contract_test
      tests/dealii_diffusion_contract.cc)
    target_link_libraries(nmopt_dealii_diffusion_contract_test
      nmopt_dealii_contract)
    if(NMOPT_USING_DEAL_II_CMAKE_PACKAGE)
      deal_ii_setup_target(nmopt_dealii_diffusion_contract_test)
    endif()
    nmopt_add_scenario_test(
      nmopt.dealii.canonical_volume_control
      nmopt_dealii_diffusion_contract_test
      canonical_volume_control
      "dealii;compiler")
    nmopt_add_scenario_test(
      nmopt.dealii.fixed_dirichlet
      nmopt_dealii_diffusion_contract_test
      fixed_dirichlet
      "dealii;compiler")
    nmopt_add_scenario_test(
      nmopt.dealii.dirichlet_control
      nmopt_dealii_diffusion_contract_test
      dirichlet_control
      "dealii;compiler")
    nmopt_add_scenario_test(
      nmopt.dealii.subdomain_observation
      nmopt_dealii_diffusion_contract_test
      subdomain_observation
      "dealii;compiler")
    nmopt_add_scenario_test(
      nmopt.dealii.neumann_boundary
      nmopt_dealii_diffusion_contract_test
      neumann_boundary
      "dealii;compiler")
    nmopt_add_scenario_test(
      nmopt.dealii.h1_control
      nmopt_dealii_diffusion_contract_test
      h1_control
      "dealii;compiler")
    nmopt_add_scenario_test(
      nmopt.dealii.coefficient_identification
      nmopt_dealii_diffusion_contract_test
      coefficient_identification
      "dealii;compiler")
    nmopt_add_scenario_test(
      nmopt.dealii.pure_neumann
      nmopt_dealii_diffusion_contract_test
      pure_neumann
      "dealii;compiler")
  endif()
endif()
