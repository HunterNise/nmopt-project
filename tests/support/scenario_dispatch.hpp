#pragma once

#include <cstddef>
#include <exception>
#include <functional>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace nmopt::test_support
{
  struct Scenario
  {
    std::string              name;
    std::string              ctest_name;
    std::vector<std::string> labels;
    unsigned int             timeout_seconds;
    std::function<void()>    run;
  };

  struct ScenarioRunResult
  {
    bool        listed = false;
    std::string executed;
  };

  inline bool
  is_manifest_token(const std::string &value)
  {
    if (value.empty())
      return false;

    for (const char character : value)
      if (!((character >= 'a' && character <= 'z') ||
            (character >= 'A' && character <= 'Z') ||
            (character >= '0' && character <= '9') || character == '_' ||
            character == '-' || character == '.'))
        return false;
    return true;
  }

  inline void
  validate_scenarios(const std::vector<Scenario> &scenarios)
  {
    if (scenarios.empty())
      throw std::invalid_argument("scenario registry is empty");

    std::set<std::string> scenario_names;
    std::set<std::string> ctest_names;
    for (const auto &scenario : scenarios)
      {
        if (!is_manifest_token(scenario.name))
          throw std::invalid_argument("invalid scenario name '" +
                                      scenario.name + "'");
        if (!is_manifest_token(scenario.ctest_name))
          throw std::invalid_argument("invalid CTest name '" +
                                      scenario.ctest_name + "'");
        if (scenario.timeout_seconds == 0)
          throw std::invalid_argument("scenario '" + scenario.name +
                                      "' has a zero timeout");
        if (scenario.labels.empty())
          throw std::invalid_argument("scenario '" + scenario.name +
                                      "' has no labels");
        for (const auto &label : scenario.labels)
          if (!is_manifest_token(label))
            throw std::invalid_argument("scenario '" + scenario.name +
                                        "' has invalid label '" + label +
                                        "'");
        if (!scenario_names.insert(scenario.name).second)
          throw std::invalid_argument("duplicate scenario name '" +
                                      scenario.name + "'");
        if (!ctest_names.insert(scenario.ctest_name).second)
          throw std::invalid_argument("duplicate CTest name '" +
                                      scenario.ctest_name + "'");
      }
  }

  inline void
  write_scenario_manifest(std::ostream &                 output,
                          const std::vector<Scenario> &scenarios)
  {
    for (const auto &scenario : scenarios)
      {
        output << "SCENARIO_DISCOVERY|" << scenario.ctest_name << '|'
               << scenario.name << '|' << scenario.timeout_seconds << '|';
        for (std::size_t index = 0; index < scenario.labels.size(); ++index)
          {
            if (index != 0)
              output << ',';
            output << scenario.labels[index];
          }
        output << '\n';
      }
  }

  inline void
  run_scenario(const Scenario &scenario)
  {
    try
      {
        scenario.run();
      }
    catch (const std::exception &exception)
      {
        throw std::runtime_error("scenario '" + scenario.name +
                                 "' failed: " + exception.what());
      }
    catch (...)
      {
        throw std::runtime_error("scenario '" + scenario.name +
                                 "' failed with a non-standard exception");
      }
  }

  inline ScenarioRunResult
  run_requested_scenarios(const int argc,
                          char **argv,
                          const std::vector<Scenario> &scenarios,
                          std::ostream &               output)
  {
    validate_scenarios(scenarios);

    if (argc == 2 && std::string(argv[1]) == "--list-scenarios")
      {
        write_scenario_manifest(output, scenarios);
        return {true, {}};
      }

    if (argc == 1)
      {
        for (const auto &scenario : scenarios)
          run_scenario(scenario);
        return {false, "all"};
      }

    if (argc == 2)
      for (const auto &scenario : scenarios)
        if (scenario.name == argv[1])
          {
            run_scenario(scenario);
            return {false, scenario.name};
          }

    std::ostringstream message;
    message << "usage: " << argv[0]
            << " [scenario|--list-scenarios]; available scenarios:";
    for (const auto &scenario : scenarios)
      message << ' ' << scenario.name;
    throw std::invalid_argument(message.str());
  }
} // namespace nmopt::test_support
