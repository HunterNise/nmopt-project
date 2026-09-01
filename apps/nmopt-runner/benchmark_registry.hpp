#pragma once

#include <array>
#include <string>
#include <string_view>

namespace nmopt::application::runner
{
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
} // namespace nmopt::application::runner
