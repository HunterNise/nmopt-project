include_guard(GLOBAL)

set(_NMOPT_SCENARIO_DISCOVERY_SCRIPT
  "${CMAKE_CURRENT_LIST_DIR}/NmoptWriteDiscoveredScenarios.cmake")

function(nmopt_discover_scenarios target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR
      "nmopt_discover_scenarios requires an existing executable target; "
      "'${target}' is not a target.")
  endif()

  get_target_property(target_type "${target}" TYPE)
  if(NOT target_type STREQUAL "EXECUTABLE")
    message(FATAL_ERROR
      "nmopt_discover_scenarios requires an executable target; "
      "'${target}' has type '${target_type}'.")
  endif()

  set(ctest_file
    "${CMAKE_CURRENT_BINARY_DIR}/${target}_discovered_scenarios.cmake")
  set(ctest_include_file
    "${CMAKE_CURRENT_BINARY_DIR}/${target}_include.cmake")

  add_custom_command(
    TARGET "${target}"
    POST_BUILD
    BYPRODUCTS "${ctest_file}"
    COMMAND "${CMAKE_COMMAND}"
      -D "NMOPT_TEST_EXECUTABLE=$<TARGET_FILE:${target}>"
      -D "NMOPT_CTEST_FILE=${ctest_file}"
      -D "NMOPT_DISCOVERY_TIMEOUT=30"
      -P "${_NMOPT_SCENARIO_DISCOVERY_SCRIPT}"
    VERBATIM)

  set_property(TARGET "${target}" APPEND PROPERTY
    LINK_DEPENDS "${_NMOPT_SCENARIO_DISCOVERY_SCRIPT}")

  file(WRITE "${ctest_include_file}"
    "if(EXISTS \"${ctest_file}\")\n"
    "  include(\"${ctest_file}\")\n"
    "else()\n"
    "  add_test(${target}_NOT_BUILT ${target}_NOT_BUILT)\n"
    "endif()\n")

  set_property(DIRECTORY APPEND PROPERTY
    TEST_INCLUDE_FILES "${ctest_include_file}")
endfunction()
