cmake_minimum_required(VERSION 3.20)

project(nmopt LANGUAGES CXX)

option(NMOPT_ENABLE_DEAL_II "Build the serial deal.II reference lowerer" ON)
option(NMOPT_ENABLE_WARNINGS "Enable project-owned compiler warnings" OFF)
option(NMOPT_ENABLE_STRICT_WARNINGS
  "Enable additional conversion and shadow warnings" OFF)
option(NMOPT_WARNINGS_AS_ERRORS "Treat project warnings as errors" OFF)
option(NMOPT_ENABLE_SANITIZERS
  "Enable address and undefined-behavior sanitizers" OFF)

if(NMOPT_ENABLE_SANITIZERS AND NMOPT_ENABLE_DEAL_II)
  message(FATAL_ERROR
    "NMOPT_ENABLE_SANITIZERS is a backend-neutral profile. "
    "Configure with NMOPT_ENABLE_DEAL_II=OFF or use the sanitize-neutral preset.")
endif()

add_library(nmopt_project_options INTERFACE)
if(NMOPT_ENABLE_WARNINGS)
  if(MSVC)
    target_compile_options(nmopt_project_options INTERFACE /W4)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(nmopt_project_options
      INTERFACE
        -Wall
        -Wextra
        -Wpedantic)
    if(NMOPT_ENABLE_STRICT_WARNINGS)
      target_compile_options(nmopt_project_options
        INTERFACE
          -Wconversion
          -Wshadow)
    endif()
  else()
    message(WARNING
      "NMOPT_ENABLE_WARNINGS is ON, but compiler-specific warning flags are not defined for ${CMAKE_CXX_COMPILER_ID}.")
  endif()
endif()

if(NMOPT_WARNINGS_AS_ERRORS)
  if(MSVC)
    target_compile_options(nmopt_project_options INTERFACE /WX)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(nmopt_project_options INTERFACE -Werror)
  endif()
endif()

if(NMOPT_ENABLE_SANITIZERS)
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    message(FATAL_ERROR
      "NMOPT_ENABLE_SANITIZERS currently supports GNU and Clang compilers only.")
  endif()
  target_compile_options(nmopt_project_options
    INTERFACE
      -fsanitize=address,undefined
      -fno-omit-frame-pointer)
  target_link_options(nmopt_project_options
    INTERFACE
      -fsanitize=address,undefined
      -fno-omit-frame-pointer)
endif()

add_library(nmopt_contract INTERFACE)
target_include_directories(nmopt_contract INTERFACE include)
target_compile_features(nmopt_contract INTERFACE cxx_std_17)

include(CTest)

function(nmopt_add_scenario_test name target scenario labels timeout)
  add_test(NAME ${name} COMMAND ${target} ${scenario})
  set_tests_properties(${name}
    PROPERTIES
      LABELS "${labels}"
      TIMEOUT "${timeout}")
endfunction()

function(nmopt_add_scenario_tests target prefix name_prefix labels timeout)
  foreach(scenario IN LISTS ARGN)
    nmopt_add_scenario_test(
      "${prefix}.${name_prefix}${scenario}"
      "${target}"
      "${scenario}"
      "${labels}"
      "${timeout}")
  endforeach()
endfunction()

if(BUILD_TESTING)
  add_executable(nmopt_contract_tests tests/reduced_dto_contract.cc)
  target_link_libraries(nmopt_contract_tests
    PRIVATE
      nmopt_contract
      nmopt_project_options)
  nmopt_add_scenario_tests(
    nmopt_contract_tests
    nmopt.contract
    ""
    "backend-neutral;contract"
    30
    layout_invariant
    backend_parameterisation)
  nmopt_add_scenario_test(
    nmopt.contract.v0
    nmopt_contract_tests
    v0_contract
    "backend-neutral;contract"
    30)
  nmopt_add_scenario_test(
    nmopt.contract.projection_compatibility
    nmopt_contract_tests
    projection_compatibility
    "backend-neutral;contract;constraint"
    30)

  add_executable(nmopt_semantic_v1_contract_test tests/semantic_v1_contract.cc)
  target_link_libraries(nmopt_semantic_v1_contract_test
    PRIVATE
      nmopt_contract
      nmopt_project_options)
  nmopt_add_scenario_tests(
    nmopt_semantic_v1_contract_test
    nmopt.semantic
    v1_
    "backend-neutral;semantic"
    30
    validation
    graph_closure
    pairing_compatibility
    incomplete_components
    reference_delta_stability)
endif()

if(NMOPT_ENABLE_DEAL_II)
  find_package(deal.II CONFIG QUIET)

  if(NOT deal.II_FOUND)
    message(FATAL_ERROR
      "NMOPT_ENABLE_DEAL_II=ON requires the official deal.II CMake package. "
      "Set deal.II_DIR to the directory containing deal.IIConfig.cmake, or "
      "configure with NMOPT_ENABLE_DEAL_II=OFF for backend-neutral tests.")
  endif()

  deal_ii_initialize_cached_variables()
  add_library(nmopt_dealii_contract INTERFACE)
  target_link_libraries(nmopt_dealii_contract INTERFACE nmopt_contract)

  if(BUILD_TESTING)
    add_executable(nmopt_dealii_diffusion_contract_test
      tests/dealii_diffusion_contract.cc)
    target_link_libraries(nmopt_dealii_diffusion_contract_test
      nmopt_dealii_contract
      nmopt_project_options)
    deal_ii_setup_target(nmopt_dealii_diffusion_contract_test)
    nmopt_add_scenario_test(
      nmopt.dealii.canonical_volume_control
      nmopt_dealii_diffusion_contract_test
      canonical_volume_control
      "dealii;compiler"
      180)
    nmopt_add_scenario_tests(
      nmopt_dealii_diffusion_contract_test
      nmopt.dealii
      ""
      "dealii;compiler"
      60
      fixed_dirichlet
      dirichlet_control
      subdomain_observation
      neumann_boundary
      h1_control
      coefficient_identification
      pure_neumann)
    nmopt_add_scenario_test(
      nmopt.dealii.projection_compatibility
      nmopt_dealii_diffusion_contract_test
      projection_compatibility
      "dealii;contract;constraint"
      60)
    nmopt_add_scenario_test(
      nmopt.dealii.backend_size_conversion
      nmopt_dealii_diffusion_contract_test
      backend_size_conversion
      "dealii;contract;adapter"
      30)
  endif()
endif()

message(STATUS "nmopt configuration summary")
message(STATUS "  generator: ${CMAKE_GENERATOR}")
message(STATUS "  build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  testing: ${BUILD_TESTING}")
message(STATUS "  deal.II backend: ${NMOPT_ENABLE_DEAL_II}")
message(STATUS "  project warnings: ${NMOPT_ENABLE_WARNINGS}")
message(STATUS "  strict warnings: ${NMOPT_ENABLE_STRICT_WARNINGS}")
message(STATUS "  warnings as errors: ${NMOPT_WARNINGS_AS_ERRORS}")
message(STATUS "  sanitizers: ${NMOPT_ENABLE_SANITIZERS}")
if(BUILD_TESTING)
  message(STATUS "  backend-neutral tests: enabled")
  if(NMOPT_ENABLE_DEAL_II)
    message(STATUS "  deal.II tests: enabled")
  else()
    message(STATUS "  deal.II tests: disabled by NMOPT_ENABLE_DEAL_II=OFF")
  endif()
endif()
