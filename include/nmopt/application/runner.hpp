#pragma once

#include "nmopt/application/artifact_writer.hpp"

#include <chrono>
#include <functional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace nmopt::application::benchmark
{
  template <typename Envelope>
  struct BenchmarkExecutionEvidenceT
  {
    using envelope_type = Envelope;

    Envelope                          envelope;
    semantic::v1::ValidationReport    diagnostics;
    BenchmarkMeasurements             measurements;
    std::vector<std::string>          selected_fields;
    std::vector<ArtifactField>        fields;
  };

  template <typename Envelope>
  struct BenchmarkRunResultT
  {
    using envelope_type = Envelope;

    BenchmarkArtifactT<Envelope> artifact;
    std::string                  document;
  };

  // This adapter owns orchestration only. The problem builder and execution
  // callback are supplied by the application/backend integration; no compiler
  // or optimization algorithm is reproduced here.
  template <typename Scenario>
  class HeadlessBenchmarkRunnerT final
  {
  public:
    using scenario_type = Scenario;

    explicit HeadlessBenchmarkRunnerT(const Scenario &scenario)
      : scenario_(&scenario)
      , harness_(scenario)
    {}

    const BenchmarkIdentity &
    identity() const noexcept
    {
      return harness_.identity();
    }

    template <typename ProblemBuilder, typename ExecutionAdapter>
    auto
    run(ProblemBuilder &&   build_problem,
        ExecutionAdapter && execute) const
      -> BenchmarkRunResultT<typename std::decay_t<
        std::invoke_result_t<ExecutionAdapter,
                             decltype(std::invoke(
                               std::declval<ProblemBuilder>(),
                               std::declval<const typename Scenario::problem_parameters_type &>())),
                             const Scenario &>>::envelope_type>
    {
      const auto start = std::chrono::steady_clock::now();
      auto       specification = std::invoke(
        std::forward<ProblemBuilder>(build_problem), scenario_->problem);
      auto evidence = std::invoke(std::forward<ExecutionAdapter>(execute),
                                  specification,
                                  *scenario_);
      const auto finish = std::chrono::steady_clock::now();

      if (scenario_->experiment.harness.measure_timings)
        {
          evidence.measurements.timing_collected = true;
          evidence.measurements.wall_seconds =
            std::chrono::duration<double>(finish - start).count();
        }

      auto artifact = harness_.finalize(std::move(evidence.envelope),
                                        std::move(evidence.diagnostics),
                                        std::move(evidence.measurements),
                                        std::move(evidence.selected_fields));
      const auto document = writer_.render(artifact, std::move(evidence.fields));
      return {std::move(artifact), document};
    }

  private:
    const Scenario *             scenario_ = nullptr;
    BenchmarkHarnessT<Scenario>  harness_;
    BenchmarkArtifactWriter      writer_;
  };
} // namespace nmopt::application::benchmark
