include_guard(GLOBAL)

set(_SCENARIO_DISCOVERY_SCRIPT
  "${CMAKE_CURRENT_LIST_DIR}/WriteDiscoveredScenarios.cmake")

function(scenario_discovery_register target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR
      "scenario_discovery_register requires an existing executable target; "
      "'${target}' is not a target.")
  endif()

  get_target_property(target_type "${target}" TYPE)
  if(NOT target_type STREQUAL "EXECUTABLE")
    message(FATAL_ERROR
      "scenario_discovery_register requires an executable target; "
      "'${target}' has type '${target_type}'.")
  endif()

  set(scenario_output_dir
    "${CMAKE_CURRENT_BINARY_DIR}/generated/scenarios/${target}")
  file(MAKE_DIRECTORY "${scenario_output_dir}")

  set(ctest_file "${scenario_output_dir}/discovered_scenarios.cmake")
  set(ctest_include_file "${scenario_output_dir}/include.cmake")

  add_custom_command(
    TARGET "${target}"
    POST_BUILD
    BYPRODUCTS "${ctest_file}"
    COMMAND "${CMAKE_COMMAND}"
      -D "SCENARIO_DISCOVERY_TEST_EXECUTABLE=$<TARGET_FILE:${target}>"
      -D "SCENARIO_DISCOVERY_CTEST_FILE=${ctest_file}"
      -D "SCENARIO_DISCOVERY_TIMEOUT=30"
      -P "${_SCENARIO_DISCOVERY_SCRIPT}"
    VERBATIM)

  set_property(TARGET "${target}" APPEND PROPERTY
    LINK_DEPENDS "${_SCENARIO_DISCOVERY_SCRIPT}")

  file(WRITE "${ctest_include_file}"
    "# Generated scenario-registration wrapper for ${target}.\n"
    "if(EXISTS \"${ctest_file}\")\n"
    "  include(\"${ctest_file}\")\n"
    "else()\n"
    "  add_test(${target}_NOT_BUILT ${target}_NOT_BUILT)\n"
    "endif()\n")

  set_property(DIRECTORY APPEND PROPERTY
    TEST_INCLUDE_FILES "${ctest_include_file}")
endfunction()
