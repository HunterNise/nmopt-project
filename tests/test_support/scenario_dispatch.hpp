#pragma once

#include <exception>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace nmopt::test_support
{
  struct Scenario
  {
    std::string           name;
    std::function<void()> run;
  };

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

  inline std::string
  run_requested_scenarios(const int argc,
                          char **argv,
                          const std::vector<Scenario> &scenarios)
  {
    if (argc == 1)
      {
        for (const auto &scenario : scenarios)
          run_scenario(scenario);
        return "all";
      }

    if (argc == 2)
      for (const auto &scenario : scenarios)
        if (scenario.name == argv[1])
          {
            run_scenario(scenario);
            return scenario.name;
          }

    std::ostringstream message;
    message << "usage: " << argv[0] << " [scenario]; available scenarios:";
    for (const auto &scenario : scenarios)
      message << ' ' << scenario.name;
    throw std::invalid_argument(message.str());
  }
} // namespace nmopt::test_support
