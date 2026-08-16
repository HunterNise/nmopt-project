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

add_library(nmopt_build_flags INTERFACE)
if(NMOPT_ENABLE_WARNINGS)
  if(MSVC)
    target_compile_options(nmopt_build_flags INTERFACE /W4)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(nmopt_build_flags
      INTERFACE
        -Wall
        -Wextra
        -Wpedantic)
    if(NMOPT_ENABLE_STRICT_WARNINGS)
      target_compile_options(nmopt_build_flags
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
    target_compile_options(nmopt_build_flags INTERFACE /WX)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(nmopt_build_flags INTERFACE -Werror)
  endif()
endif()

if(NMOPT_ENABLE_SANITIZERS)
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    message(FATAL_ERROR
      "NMOPT_ENABLE_SANITIZERS currently supports GNU and Clang compilers only.")
  endif()
  target_compile_options(nmopt_build_flags
    INTERFACE
      -fsanitize=address,undefined
      -fno-omit-frame-pointer)
  target_link_options(nmopt_build_flags
    INTERFACE
      -fsanitize=address,undefined
      -fno-omit-frame-pointer)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
file(MAKE_DIRECTORY
  "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
  "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
  "${CMAKE_BINARY_DIR}/generated/scenarios")

add_library(nmopt_contract INTERFACE)
target_include_directories(nmopt_contract INTERFACE include)
target_compile_features(nmopt_contract INTERFACE cxx_std_17)

include(CTest)
include(cmake/ScenarioDiscovery.cmake)

function(nmopt_add_neutral_executable target source_file)
  add_executable("${target}" "${source_file}")
  target_link_libraries("${target}"
    PRIVATE
      nmopt_contract
      nmopt_build_flags)
endfunction()

function(nmopt_add_neutral_test target source_file)
  nmopt_add_neutral_executable("${target}" "${source_file}")
  scenario_discovery_register("${target}")
endfunction()

function(nmopt_add_dealii_executable target source_file)
  add_executable("${target}" "${source_file}")
  target_link_libraries("${target}"
    nmopt_dealii_contract
    nmopt_build_flags)
  deal_ii_setup_target("${target}")
endfunction()

function(nmopt_add_dealii_test target source_file)
  nmopt_add_dealii_executable("${target}" "${source_file}")
  scenario_discovery_register("${target}")
endfunction()

if(BUILD_TESTING)
  nmopt_add_neutral_test(
    nmopt_reduced_dto_contract_test
    tests/reduced_dto_contract.cc)

  nmopt_add_neutral_test(
    nmopt_supplied_otd_contract_test
    tests/supplied_otd_contract.cc)

  nmopt_add_neutral_test(
    nmopt_supplied_otd_reference_contract_test
    tests/supplied_linear_quadratic_otd_contract.cc)

  nmopt_add_neutral_test(
    nmopt_quadratic_kkt_contract_test
    tests/quadratic_kkt_contract.cc)

  nmopt_add_neutral_test(
    nmopt_quadratic_kkt_reference_contract_test
    tests/quadratic_kkt_reference_contract.cc)

  nmopt_add_neutral_test(
    nmopt_complementarity_contract_test
    tests/complementarity_contract.cc)

  nmopt_add_neutral_test(
    nmopt_semantic_v1_contract_test
    tests/semantic_v1_contract.cc)
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
    nmopt_add_dealii_test(
      nmopt_dealii_diffusion_contract_test
      tests/dealii_diffusion_contract.cc)

    nmopt_add_dealii_test(
      nmopt_dealii_kkt_contract_test
      tests/dealii_kkt_contract.cc)

    nmopt_add_dealii_test(
      nmopt_dealii_trace_hhalf_metric_contract_test
      tests/dealii_trace_hhalf_metric_contract.cc)
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
