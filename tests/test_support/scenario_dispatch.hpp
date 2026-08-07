#pragma once

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

  inline std::string
  run_requested_scenarios(const int argc,
                          char **argv,
                          const std::vector<Scenario> &scenarios)
  {
    if (argc == 1)
      {
        for (const auto &scenario : scenarios)
          scenario.run();
        return "all";
      }

    if (argc == 2)
      for (const auto &scenario : scenarios)
        if (scenario.name == argv[1])
          {
            scenario.run();
            return scenario.name;
          }

    std::ostringstream message;
    message << "usage: " << argv[0] << " [scenario]; available scenarios:";
    for (const auto &scenario : scenarios)
      message << ' ' << scenario.name;
    throw std::invalid_argument(message.str());
  }
} // namespace nmopt::test_support
