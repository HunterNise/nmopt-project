#pragma once

#include "nmopt/contract/linalg.hpp"
#include "nmopt/semantic/v1/validation.hpp"

#include <cstddef>
#include <sstream>
#include <string>

namespace nmopt::test_support
{
  inline const char *
  diagnostic_category_name(
    const semantic::v1::DiagnosticCategory category)
  {
    using semantic::v1::DiagnosticCategory;
    switch (category)
      {
        case DiagnosticCategory::structural:
          return "structural";
        case DiagnosticCategory::analytical_policy:
          return "analytical_policy";
        case DiagnosticCategory::lowerability:
          return "lowerability";
        case DiagnosticCategory::formulation_capability:
          return "formulation_capability";
      }
    return "unknown";
  }

  inline void
  require_exact_diagnostic(
    const semantic::v1::ValidationReport &report,
    const semantic::v1::DiagnosticCategory category,
    const std::string &component_id,
    const std::string &capability,
    const std::string &description)
  {
    std::size_t matches = 0;
    for (const auto &diagnostic : report.diagnostics())
      if (diagnostic.category == category &&
          diagnostic.component_id == component_id &&
          diagnostic.capability == capability)
        ++matches;
    if (matches == 1)
      return;

    std::ostringstream message;
    message << description << ": expected exactly one diagnostic {"
            << diagnostic_category_name(category) << ", " << component_id
            << ", " << capability << "}, found " << matches
            << "; actual diagnostics:";
    for (const auto &diagnostic : report.diagnostics())
      message << " {" << diagnostic_category_name(diagnostic.category) << ", "
              << diagnostic.component_id << ", " << diagnostic.capability
              << "}";
    if (report.diagnostics().empty())
      message << " none";
    throw contract::ContractError(message.str());
  }
} // namespace nmopt::test_support
