cmake_minimum_required(VERSION 3.20)

project(nmopt LANGUAGES CXX)

option(NMOPT_ENABLE_DEAL_II "Build the serial deal.II reference lowerer" ON)
option(NMOPT_ENABLE_WARNINGS "Enable project-owned compiler warnings" OFF)
option(NMOPT_ENABLE_STRICT_WARNINGS
  "Enable additional conversion and shadow warnings" OFF)
option(NMOPT_WARNINGS_AS_ERRORS "Treat project warnings as errors" OFF)
option(NMOPT_ENABLE_SANITIZERS
  "Enable address and undefined-behavior sanitizers" OFF)

if(CMAKE_BUILD_TYPE)
  string(TOLOWER "${CMAKE_BUILD_TYPE}" NMOPT_BUILD_CONFIGURATION)
else()
  set(NMOPT_BUILD_CONFIGURATION "unknown")
endif()
if(NMOPT_ENABLE_DEAL_II)
  set(NMOPT_BUILD_PROFILE_NAME "${NMOPT_BUILD_CONFIGURATION}-dealii")
else()
  set(NMOPT_BUILD_PROFILE_NAME "${NMOPT_BUILD_CONFIGURATION}-neutral")
endif()

if(NMOPT_ENABLE_SANITIZERS AND NMOPT_ENABLE_DEAL_II)
  message(FATAL_ERROR
    "NMOPT_ENABLE_SANITIZERS is a backend-neutral profile. "
    "Configure with NMOPT_ENABLE_DEAL_II=OFF or use the sanitize-neutral preset.")
endif()

add_library(nmopt_build_flags INTERFACE)
target_compile_definitions(nmopt_build_flags
  INTERFACE
    NMOPT_COMPILED_BUILD_PROFILE=\"${NMOPT_BUILD_PROFILE_NAME}\")
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

function(nmopt_add_python_test test_name script labels timeout)
  add_test(
    NAME "${test_name}"
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/${script}")
  set_tests_properties(
    "${test_name}"
    PROPERTIES
      LABELS "${labels}"
      TIMEOUT "${timeout}"
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
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
    tests/contract/reduced_dto_contract.cc)

  nmopt_add_neutral_test(
    nmopt_executable_model_contract_test
    tests/contract/executable_model_contract.cc)

  nmopt_add_neutral_test(
    nmopt_supplied_otd_contract_test
    tests/contract/supplied_otd_contract.cc)

  nmopt_add_neutral_test(
    nmopt_supplied_otd_reference_contract_test
    tests/contract/supplied_linear_quadratic_otd_contract.cc)

  nmopt_add_neutral_test(
    nmopt_quadratic_kkt_contract_test
    tests/contract/quadratic_kkt_contract.cc)

  nmopt_add_neutral_test(
    nmopt_quadratic_kkt_reference_contract_test
    tests/contract/quadratic_kkt_reference_contract.cc)

  nmopt_add_neutral_test(
    nmopt_complementarity_contract_test
    tests/contract/complementarity_contract.cc)

  nmopt_add_neutral_test(
    nmopt_pdas_contract_test
    tests/contract/pdas_contract.cc)

  nmopt_add_neutral_test(
    nmopt_semantic_v1_contract_test
    tests/semantic/semantic_v1_contract.cc)

  nmopt_add_neutral_test(
    nmopt_application_boundary_contract_test
    tests/application/application_boundary_contract.cc)

  nmopt_add_neutral_test(
    nmopt_runner_contract_test
    tests/application/runner_contract.cc)

  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    nmopt_add_python_test(
      nmopt.tools.chapter6_report
      tests/tools/chapter6_report_contract.py
      "python;tool;report"
      30)
    nmopt_add_python_test(
      nmopt.tools.artifact_records
      tests/tools/artifact_records_contract.py
      "python;tool;artifact"
      30)

    execute_process(
      COMMAND "${Python3_EXECUTABLE}" -c "import numpy, matplotlib, meshio"
      RESULT_VARIABLE nmopt_plotting_import_result
      OUTPUT_QUIET
      ERROR_QUIET)
    if(nmopt_plotting_import_result EQUAL 0)
      nmopt_add_python_test(
        nmopt.tools.postprocess_configuration
        tests/tools/postprocess_configuration_contract.py
        "python;tool;postprocess"
        60)
      nmopt_add_python_test(
        nmopt.tools.postprocess_comparison_plan
        tests/tools/postprocess_comparison_plan_contract.py
        "python;tool;postprocess"
        60)
      nmopt_add_python_test(
        nmopt.tools.postprocess_history
        tests/tools/postprocess_history_contract.py
        "python;tool;postprocess"
        60)
      nmopt_add_python_test(
        nmopt.tools.postprocess_boundary_field
        tests/tools/postprocess_boundary_field_contract.py
        "python;tool;postprocess"
        60)
      nmopt_add_python_test(
        nmopt.tools.postprocess_render_policy
        tests/tools/postprocess_render_policy_contract.py
        "python;tool;postprocess"
        60)
      set(nmopt_python_tool_test_status "stdlib tooling + postprocessing")
    else()
      set(nmopt_python_tool_test_status
          "stdlib tooling only (plotting dependencies unavailable)")
    endif()
  else()
    set(nmopt_python_tool_test_status
        "disabled (Python interpreter not found)")
  endif()
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

  nmopt_add_dealii_executable(
    nmopt_runner
    apps/nmopt-runner/main.cc)

  nmopt_add_dealii_executable(
    nmopt_external_poisson_forward
    tests/dealii/external_poisson_forward.cc)

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

  nmopt_add_dealii_executable(
    nmopt_external_step4_minimal_problem_a
    apps/external-dealii/step-4/minimal/problem_a.cc)

  nmopt_add_dealii_executable(
    nmopt_external_step4_minimal_problem_b
    apps/external-dealii/step-4/minimal/problem_b.cc)

  if(BUILD_TESTING)
    # This test includes only the adapted tutorial source and the standard
    # scenario-discovery helper. It intentionally has no nmopt library link.
    add_executable(
      nmopt_external_step4_native_contract_test
      tests/dealii/external_step4_native_contract.cc)
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

    add_test(
      NAME nmopt.external_application.standalone_forward
      COMMAND nmopt_external_poisson_forward)
    set_tests_properties(
      nmopt.external_application.standalone_forward
      PROPERTIES LABELS "dealii;application;external" TIMEOUT 30)

    nmopt_add_dealii_test(
      nmopt_external_application_dealii_contract_test
      tests/application/external_application_dealii_contract.cc)

    nmopt_add_dealii_test(
      nmopt_dealii_diffusion_contract_test
      tests/compiler/dealii/dealii_compiler_contract.cc)

    nmopt_add_dealii_test(
      nmopt_dealii_backend_contract_test
      tests/dealii/dealii_backend_contract.cc)

    nmopt_add_dealii_test(
      nmopt_dealii_trace_hhalf_metric_contract_test
      tests/dealii/dealii_trace_hhalf_metric_contract.cc)

    nmopt_add_dealii_test(
      nmopt_application_chapter6_b1_dealii_contract_test
      tests/application/chapter6_b1_dealii_contract.cc)

    nmopt_add_dealii_test(
      nmopt_compiled_application_view_dealii_contract_test
      tests/application/compiled_application_view_dealii_contract.cc)

    nmopt_add_dealii_test(
    nmopt_application_chapter6_b2_dealii_contract_test
      tests/application/chapter6_b2_dealii_contract.cc)

    nmopt_add_dealii_test(
      nmopt_application_parameter_files_dealii_contract_test
      tests/application/parameter_files_dealii_contract.cc)
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
  message(STATUS "  Python tool tests: ${nmopt_python_tool_test_status}")
  message(STATUS "  backend-neutral tests: enabled")
  if(NMOPT_ENABLE_DEAL_II)
    message(STATUS "  deal.II tests: enabled")
  else()
    message(STATUS "  deal.II tests: disabled by NMOPT_ENABLE_DEAL_II=OFF")
  endif()
endif()
