#pragma once

#include "run_set_plan.hpp"

#include <array>
#include <functional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace nmopt::application::runner
{
  struct ParameterFile;
  struct ResolvedRunConfiguration;
  class RunSetManifest;

  struct BenchmarkRegistration
  {
    std::string_view id;
    std::string_view parameter_benchmark_id;
    std::string_view default_parameter_file;
  };

  inline const std::array<BenchmarkRegistration, 2> &
  benchmark_registrations()
  {
    static const std::array<BenchmarkRegistration, 2> registrations{{
      {"b1",
       "chapter-6.b1.distributed-laplace",
       "parameters/chapter-6/b1/authoritative.prm"},
      {"b2",
       "chapter-6.b2.graetz-flow",
       "parameters/chapter-6/b2/authoritative.prm"}}};
    return registrations;
  }

  inline const BenchmarkRegistration *
  find_benchmark_registration(const std::string_view id)
  {
    for (const auto &registration : benchmark_registrations())
      if (registration.id == id)
        return &registration;
    return nullptr;
  }

  inline const BenchmarkRegistration *
  find_benchmark_registration_for_parameter_id(
    const std::string_view parameter_benchmark_id)
  {
    for (const auto &registration : benchmark_registrations())
      if (registration.parameter_benchmark_id == parameter_benchmark_id)
        return &registration;
    return nullptr;
  }

  inline std::string
  registered_benchmark_ids()
  {
    std::string result;
    for (const auto &registration : benchmark_registrations())
      {
        if (!result.empty())
          result += ", ";
        result += registration.id;
      }
    return result;
  }

  using BenchmarkSelectionFilters =
    std::vector<std::pair<std::string, std::string>>;
  using BenchmarkArtifactPlanner = std::function<
    std::vector<std::string>(const RunSetPlan &)>;
  using BenchmarkExecutionCallback = std::function<
    bool(const ResolvedRunConfiguration &,
         const std::vector<std::string> &,
         const RunSetPlan &,
         const ParameterFile &,
         RunSetManifest &)>;

  struct BenchmarkExecutionRegistration
  {
    const BenchmarkRegistration *metadata = nullptr;
    BenchmarkArtifactPlanner      artifact_planner;
    BenchmarkExecutionCallback    execute;
  };
} // namespace nmopt::application::runner
