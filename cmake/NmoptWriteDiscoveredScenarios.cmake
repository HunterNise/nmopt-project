cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED NMOPT_TEST_EXECUTABLE OR
   NOT DEFINED NMOPT_CTEST_FILE OR
   NOT DEFINED NMOPT_DISCOVERY_TIMEOUT)
  message(FATAL_ERROR
    "Scenario discovery requires NMOPT_TEST_EXECUTABLE, NMOPT_CTEST_FILE, "
    "and NMOPT_DISCOVERY_TIMEOUT.")
endif()

# Discovery reads static metadata and is not a leak-checking execution. This
# also keeps instrumented discovery usable in ptrace-based environments.
set(ENV{ASAN_OPTIONS} "detect_leaks=0")

execute_process(
  COMMAND "${NMOPT_TEST_EXECUTABLE}" --list-scenarios
  TIMEOUT "${NMOPT_DISCOVERY_TIMEOUT}"
  RESULT_VARIABLE discovery_result
  OUTPUT_VARIABLE scenario_manifest
  ERROR_VARIABLE discovery_error
  OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT discovery_result EQUAL 0)
  message(FATAL_ERROR
    "Failed to discover scenarios from '${NMOPT_TEST_EXECUTABLE}'.\n"
    "Exit result: ${discovery_result}\n"
    "Standard error: ${discovery_error}")
endif()

if(scenario_manifest STREQUAL "")
  message(FATAL_ERROR
    "Scenario discovery from '${NMOPT_TEST_EXECUTABLE}' returned no tests.")
endif()

string(REPLACE "\r\n" "\n" scenario_manifest "${scenario_manifest}")
string(REPLACE "\r" "\n" scenario_manifest "${scenario_manifest}")
string(REPLACE "\n" ";" manifest_lines "${scenario_manifest}")

set(ctest_content "")
set(discovered_names "")
set(scenario_count 0)

foreach(manifest_line IN LISTS manifest_lines)
  string(REPLACE "|" ";" fields "${manifest_line}")
  list(LENGTH fields field_count)
  if(NOT field_count EQUAL 5)
    message(FATAL_ERROR
      "Malformed scenario manifest line from '${NMOPT_TEST_EXECUTABLE}': "
      "'${manifest_line}'.")
  endif()

  list(GET fields 0 marker)
  list(GET fields 1 ctest_name)
  list(GET fields 2 scenario_name)
  list(GET fields 3 timeout)
  list(GET fields 4 label_field)

  if(NOT marker STREQUAL "NMOPT_SCENARIO" OR
     NOT ctest_name MATCHES "^[A-Za-z0-9_.-]+$" OR
     NOT scenario_name MATCHES "^[A-Za-z0-9_.-]+$" OR
     NOT timeout MATCHES "^[1-9][0-9]*$" OR
     label_field STREQUAL "")
    message(FATAL_ERROR
      "Invalid scenario metadata from '${NMOPT_TEST_EXECUTABLE}': "
      "'${manifest_line}'.")
  endif()

  if("${ctest_name}" IN_LIST discovered_names)
    message(FATAL_ERROR
      "Duplicate discovered CTest name '${ctest_name}' from "
      "'${NMOPT_TEST_EXECUTABLE}'.")
  endif()
  list(APPEND discovered_names "${ctest_name}")

  string(REPLACE "," ";" labels "${label_field}")
  foreach(label IN LISTS labels)
    if(NOT label MATCHES "^[A-Za-z0-9_.-]+$")
      message(FATAL_ERROR
        "Invalid label '${label}' for discovered test '${ctest_name}'.")
    endif()
  endforeach()

  string(APPEND ctest_content
    "add_test([=[${ctest_name}]=] "
    "[=[${NMOPT_TEST_EXECUTABLE}]=] [=[${scenario_name}]=])\n"
    "set_tests_properties([=[${ctest_name}]=] PROPERTIES "
    "LABELS [=[${labels}]=] TIMEOUT [=[${timeout}]=])\n")
  math(EXPR scenario_count "${scenario_count} + 1")
endforeach()

if(scenario_count EQUAL 0)
  message(FATAL_ERROR
    "Scenario discovery from '${NMOPT_TEST_EXECUTABLE}' returned no tests.")
endif()

file(WRITE "${NMOPT_CTEST_FILE}" "${ctest_content}")
