#pragma once

#include "../../apps/external-dealii/step-4/diagnostics/instrumentation.hpp"

#include <chrono>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace external_dealii_step4_test
{
  using Instrumentation = external_dealii_step4::Instrumentation;

  struct InstrumentationView
  {
    const char *             path;
    const Instrumentation *  value;
  };

  class EvidenceGuard;
  inline thread_local EvidenceGuard *active_guard = nullptr;

  inline const char *
  solve_role_name(const external_dealii_step4::SolveRole role)
  {
    return role == external_dealii_step4::SolveRole::state ? "state" :
                                                               "adjoint";
  }

  inline std::filesystem::path
  create_unique_artifact_root(const std::filesystem::path &parent,
                              const std::string            &prefix)
  {
    std::filesystem::create_directories(parent);
    const auto now = std::chrono::system_clock::now().time_since_epoch();
    const auto run_id = std::to_string(
      std::chrono::duration_cast<std::chrono::microseconds>(now).count());
    for (std::size_t attempt = 0;; ++attempt)
      {
        const auto suffix = attempt == 0 ? std::string{} :
                                           "-" + std::to_string(attempt);
        const auto root = parent / (prefix + "-" + run_id + suffix);
        if (std::filesystem::create_directories(root))
          return root;
      }
  }

  inline void
  begin_artifact(const std::filesystem::path &root, const char *const scenario)
  {
    std::filesystem::create_directories(root);
    std::ofstream status(root / "status.txt");
    if (!status)
      throw std::runtime_error("could not open Step-4 evidence status");
    status << "status running\nscenario " << scenario << '\n';
  }

  inline void
  write_counters(const std::filesystem::path &                  root,
                 const std::vector<InstrumentationView> &views)
  {
    const std::pair<const char *, std::size_t Instrumentation::*> counters[] = {
      {"assembly_calls", &Instrumentation::assembly_calls},
      {"state_solve_calls", &Instrumentation::state_solve_calls},
      {"adjoint_solve_calls", &Instrumentation::adjoint_solve_calls},
      {"solve_failures", &Instrumentation::solve_failures},
      {"objective_calls", &Instrumentation::objective_calls},
      {"objective_derivative_calls",
       &Instrumentation::objective_derivative_calls},
      {"residual_calls", &Instrumentation::residual_calls},
      {"residual_jvp_calls", &Instrumentation::residual_jvp_calls},
      {"residual_vjp_calls", &Instrumentation::residual_vjp_calls},
      {"control_vjp_calls", &Instrumentation::control_vjp_calls},
      {"explicit_matrix_vmult_calls",
       &Instrumentation::explicit_matrix_vmult_calls},
      {"explicit_matrix_tvmult_calls",
       &Instrumentation::explicit_matrix_tvmult_calls},
      {"value_evaluations", &Instrumentation::value_evaluations},
      {"derivative_augmentations",
       &Instrumentation::derivative_augmentations},
      {"output_calls", &Instrumentation::output_calls},
      {"metric_apply_calls", &Instrumentation::metric_apply_calls},
      {"metric_inverse_apply_calls",
       &Instrumentation::metric_inverse_apply_calls}};
    std::ofstream output(root / "counters.csv");
    if (!output)
      throw std::runtime_error("could not open Step-4 evidence counters");
    output << "path,counter,value\n";
    for (const auto &view : views)
      for (const auto &[name, member] : counters)
        output << view.path << ',' << name << ',' << view.value->*member
               << '\n';
  }

  inline void
  write_solve_records(const std::filesystem::path &                  root,
                      const std::vector<InstrumentationView> &views)
  {
    std::ofstream output(root / "solve-records.csv");
    if (!output)
      throw std::runtime_error("could not open Step-4 evidence solve records");
    output << "path,status,role,iterations,initial_residual,final_residual\n"
           << std::setprecision(17);
    for (const auto &view : views)
      {
        const auto &instrumentation = *view.value;
        for (const auto &record : instrumentation.solve_records)
          output << view.path << ",success," << solve_role_name(record.role)
                 << ',' << record.iterations << ',' << record.initial_residual
                 << ',' << record.final_residual << '\n';
        for (const auto &record : instrumentation.solve_failure_records)
          output << view.path << ",failure," << solve_role_name(record.role)
                 << ',' << record.iterations << ",," << record.final_residual
                 << '\n';
      }
  }

  inline void
  write_problem_b_matrix_actions(
    const std::filesystem::path &                  root,
    const std::vector<InstrumentationView> &views)
  {
    std::ofstream output(root / "matrix-actions.csv");
    if (!output)
      throw std::runtime_error(
        "could not open Step-4 evidence matrix actions");
    output << "path,action,purpose\n";
    for (const auto &view : views)
      for (const auto &record : view.value->problem_b_matrix_actions)
        output << view.path << ','
               << external_dealii_step4::problem_b_matrix_action_name(
                    record.action)
               << ','
               << external_dealii_step4::problem_b_matrix_purpose_name(
                    record.purpose)
               << '\n';
  }

  inline void
  write_metric_solve_records(const std::filesystem::path &                  root,
                             const std::vector<InstrumentationView> &views)
  {
    std::ofstream output(root / "metric-solve-records.csv");
    if (!output)
      throw std::runtime_error(
        "could not open Step-4 evidence metric solve records");
    output << "path,status,purpose,iterations,initial_residual,final_residual\n"
           << std::setprecision(17);
    for (const auto &view : views)
      for (const auto &record : view.value->metric_solve_records)
        output << view.path << ',' << (record.converged ? "success" : "failure")
               << ','
               << external_dealii_step4::metric_solve_purpose_name(
                    record.purpose)
               << ',' << record.iterations << ',' << record.initial_residual
               << ',' << record.final_residual << '\n';
  }

  inline void
  write_complete(const std::filesystem::path &root,
                 const std::vector<InstrumentationView> &views)
  {
    write_counters(root, views);
    write_solve_records(root, views);
    write_problem_b_matrix_actions(root, views);
    write_metric_solve_records(root, views);
    std::ofstream status(root / "status.txt");
    if (!status)
      throw std::runtime_error("could not update Step-4 evidence status");
    status << "status complete\n";
  }

  inline void
  write_failure(const std::filesystem::path &root,
                const char *const              scenario,
                const std::string &             message,
                const std::vector<InstrumentationView> &views) noexcept
  {
    try
      {
        std::ofstream failure(root / "failure.txt");
        if (failure)
          failure << "status failed\nscenario " << scenario << '\n'
                  << "error " << message << '\n';
        write_counters(root, views);
        write_solve_records(root, views);
        write_problem_b_matrix_actions(root, views);
        write_metric_solve_records(root, views);
        std::ofstream status(root / "status.txt");
        if (status)
          status << "status failed\nscenario " << scenario << '\n';
      }
    catch (...)
      {}
  }

  inline std::string
  current_exception_message()
  {
    try
      {
        if (const auto exception = std::current_exception(); exception)
          std::rethrow_exception(exception);
      }
    catch (const std::exception &exception)
      {
        return exception.what();
      }
    catch (...)
      {
        return "non-standard exception";
      }
    return "scenario terminated before completion";
  }

  class EvidenceGuard final
  {
  public:
    EvidenceGuard(const std::filesystem::path &root,
                  const char *const              scenario,
                  const std::vector<InstrumentationView> &views)
      : root_(root)
      , scenario_(scenario)
      , views_(views)
      , previous_(active_guard)
    {
      begin_artifact(root_, scenario_);
      active_guard = this;
    }

    ~EvidenceGuard() noexcept
    {
      if (!completed_ && !failure_recorded_)
        write_failure(root_,
                      scenario_,
                      failure_message_.empty() ? current_exception_message() :
                                                 failure_message_,
                      views_);
      active_guard = previous_;
    }

    void
    complete()
    {
      write_complete(root_, views_);
      completed_ = true;
    }

    void
    fail(const std::string &message) noexcept
    {
      failure_message_ = message;
      failure_recorded_ = true;
      write_failure(root_, scenario_, failure_message_, views_);
    }

    void
    fail_current_exception() noexcept
    {
      fail(current_exception_message());
    }

  private:
    std::filesystem::path              root_;
    const char *                        scenario_;
    std::vector<InstrumentationView>   views_;
    EvidenceGuard *                    previous_;
    std::string                        failure_message_;
    bool                               completed_ = false;
    bool                               failure_recorded_ = false;
  };

  inline void
  note_failure(const std::string &message) noexcept
  {
    if (active_guard != nullptr)
      active_guard->fail(message);
  }
} // namespace external_dealii_step4_test
