#pragma once

#include "nmopt/application/chapter6.hpp"
#include "nmopt/semantic/v1/validation.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::application::benchmark
{
  struct BenchmarkIdentity
  {
    unsigned int             schema_version = 1;
    std::string              scenario_id;
    std::string              recipe_id;
    std::string              output_id;
    std::string              source_reference;
    std::string              source_revision;
    std::string              build_profile;
    std::string              artifact_directory;
    bool                     deterministic = true;
    std::vector<std::string> requirements;
  };

  struct BenchmarkMeasurements
  {
    bool        timing_collected = false;
    double      wall_seconds = 0.0;
    double      cpu_seconds = 0.0;
    bool        memory_collected = false;
    std::size_t peak_memory_bytes = 0;
  };

  inline void
  validate_benchmark_identity(const BenchmarkIdentity &identity)
  {
    if (identity.schema_version == 0)
      throw std::invalid_argument(
        "benchmark identity needs a positive schema version");
    if (identity.scenario_id.empty())
      throw std::invalid_argument("benchmark identity needs a scenario id");
    if (identity.recipe_id.empty())
      throw std::invalid_argument("benchmark identity needs a recipe id");
    if (identity.output_id.empty())
      throw std::invalid_argument("benchmark identity needs an output id");
    if (identity.source_reference.empty())
      throw std::invalid_argument(
        "benchmark identity needs a source reference");
    if (identity.source_revision.empty())
      throw std::invalid_argument(
        "benchmark identity needs a source revision");
    if (identity.build_profile.empty())
      throw std::invalid_argument(
        "benchmark identity needs a build profile");
    if (identity.artifact_directory.empty())
      throw std::invalid_argument(
        "benchmark identity needs an artifact directory");
  }

  inline void
  validate_benchmark_measurements(const BenchmarkMeasurements &measurements)
  {
    if (measurements.timing_collected &&
        (!std::isfinite(measurements.wall_seconds) ||
         !std::isfinite(measurements.cpu_seconds) ||
         measurements.wall_seconds < 0.0 ||
         measurements.cpu_seconds < 0.0))
      throw std::invalid_argument(
        "collected benchmark timings must be finite and nonnegative");
  }

  template <typename Envelope>
  class BenchmarkArtifactT final
  {
  public:
    using envelope_type = Envelope;

    BenchmarkArtifactT(BenchmarkIdentity                 identity,
                       semantic::v1::ValidationReport   diagnostics,
                       Envelope                         envelope,
                       BenchmarkMeasurements            measurements,
                       std::vector<std::string>          selected_fields = {})
      : identity_(std::move(identity))
      , diagnostics_(std::move(diagnostics))
      , envelope_(std::move(envelope))
      , measurements_(std::move(measurements))
      , selected_fields_(std::move(selected_fields))
    {
      validate_benchmark_identity(identity_);
      validate_benchmark_measurements(measurements_);
    }

    const BenchmarkIdentity &
    identity() const noexcept
    {
      return identity_;
    }

    const semantic::v1::ValidationReport &
    diagnostics() const noexcept
    {
      return diagnostics_;
    }

    const Envelope &
    envelope() const noexcept
    {
      return envelope_;
    }

    const BenchmarkMeasurements &
    measurements() const noexcept
    {
      return measurements_;
    }

    const std::vector<std::string> &
    selected_fields() const noexcept
    {
      return selected_fields_;
    }

  private:
    BenchmarkIdentity                 identity_;
    semantic::v1::ValidationReport    diagnostics_;
    Envelope                          envelope_;
    BenchmarkMeasurements             measurements_;
    std::vector<std::string>          selected_fields_;
  };

  // This helper intentionally accepts the concrete Chapter 6 scenario shape.
  // It is the one place where scenario metadata is projected into the B0
  // artifact identity; serializers and backend runners consume the result.
  template <typename Scenario>
  inline BenchmarkIdentity
  make_benchmark_identity(const Scenario &scenario)
  {
    scenario.validate();
    return {1,
            scenario.metadata.id,
            scenario.metadata.recipe_id,
            scenario.experiment.scenario_output_id,
            scenario.experiment.source_reference,
            scenario.experiment.source_revision,
            scenario.experiment.build_profile,
            scenario.experiment.harness.artifact_directory,
            scenario.experiment.harness.deterministic,
            scenario.metadata.requirements};
  }

  template <typename Scenario>
  class BenchmarkHarnessT final
  {
  public:
    using scenario_type = Scenario;

    explicit BenchmarkHarnessT(const Scenario &scenario)
      : identity_(make_benchmark_identity(scenario))
    {
      validate_benchmark_identity(identity_);
    }

    const BenchmarkIdentity &
    identity() const noexcept
    {
      return identity_;
    }

    template <typename Envelope>
    BenchmarkArtifactT<Envelope>
    finalize(Envelope                         envelope,
             semantic::v1::ValidationReport   diagnostics = {},
             BenchmarkMeasurements            measurements = {},
             std::vector<std::string>          selected_fields = {}) const
    {
      return BenchmarkArtifactT<Envelope>{identity_,
                                          std::move(diagnostics),
                                          std::move(envelope),
                                          std::move(measurements),
                                          std::move(selected_fields)};
    }

  private:
    BenchmarkIdentity identity_;
  };
} // namespace nmopt::application::benchmark
