#pragma once

#include "nmopt/application/harness.hpp"

#include <algorithm>
#include <iomanip>
#include <locale>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::application::benchmark
{
  struct ArtifactField
  {
    std::string key;
    std::string value;
  };

  inline const char *
  diagnostic_category_name(const semantic::v1::DiagnosticCategory category)
  {
    switch (category)
      {
        case semantic::v1::DiagnosticCategory::structural:
          return "structural";
        case semantic::v1::DiagnosticCategory::analytical_policy:
          return "analytical_policy";
        case semantic::v1::DiagnosticCategory::lowerability:
          return "lowerability";
        case semantic::v1::DiagnosticCategory::formulation_capability:
          return "formulation_capability";
      }
    return "unknown";
  }

  inline std::string
  escape_artifact_value(const std::string &value)
  {
    std::string escaped;
    escaped.reserve(value.size());
    for (const char character : value)
      switch (character)
        {
          case '\\':
            escaped += "\\\\";
            break;
          case '\n':
            escaped += "\\n";
            break;
          case '\r':
            escaped += "\\r";
            break;
          case '\t':
            escaped += "\\t";
            break;
          case '=':
            escaped += "\\=";
            break;
          default:
            escaped += character;
            break;
        }
    return escaped;
  }

  inline void
  validate_artifact_fields(const std::vector<ArtifactField> &fields)
  {
    for (std::size_t index = 0; index < fields.size(); ++index)
      {
        if (fields[index].key.empty())
          throw std::invalid_argument("artifact fields need nonempty keys");
        if (fields[index].key.find_first_of("=\\\n\r\t") !=
            std::string::npos)
          throw std::invalid_argument(
            "artifact field keys cannot contain separators");
        if (index != 0 && fields[index - 1].key >= fields[index].key)
          throw std::invalid_argument(
            "artifact fields must be sorted by unique key");
      }
  }

  class BenchmarkArtifactWriter final
  {
  public:
    template <typename Envelope>
    std::string
    render(const BenchmarkArtifactT<Envelope> &artifact,
           std::vector<ArtifactField>          fields = {}) const
    {
      std::sort(fields.begin(), fields.end(), [](const ArtifactField &left,
                                                 const ArtifactField &right) {
        return left.key < right.key;
      });
      validate_artifact_fields(fields);

      std::ostringstream output;
      output.imbue(std::locale::classic());
      output << std::setprecision(17);
      write_line(output, "artifact.schema", "nmopt-benchmark-v1");
      write_line(output,
                 "identity.schema_version",
                 number(artifact.identity().schema_version));
      write_line(output, "identity.scenario_id", artifact.identity().scenario_id);
      write_line(output, "identity.recipe_id", artifact.identity().recipe_id);
      write_line(output, "identity.output_id", artifact.identity().output_id);
      write_line(output,
                 "identity.source_reference",
                 artifact.identity().source_reference);
      write_line(output,
                 "identity.source_revision",
                 artifact.identity().source_revision);
      write_line(output,
                 "identity.build_profile",
                 artifact.identity().build_profile);
      write_line(output,
                 "identity.artifact_directory",
                 artifact.identity().artifact_directory);
      write_line(output,
                 "identity.deterministic",
                 artifact.identity().deterministic ? "true" : "false");
      write_line(output,
                 "identity.requirement_count",
                 number(artifact.identity().requirements.size()));
      for (std::size_t index = 0;
           index < artifact.identity().requirements.size();
           ++index)
        write_line(output,
                   "identity.requirement[" + number(index) + "]",
                   artifact.identity().requirements[index]);

      write_line(output,
                 "diagnostics.valid",
                 artifact.diagnostics().valid() ? "true" : "false");
      write_line(output,
                 "diagnostics.count",
                 number(artifact.diagnostics().diagnostics().size()));
      for (std::size_t index = 0;
           index < artifact.diagnostics().diagnostics().size();
           ++index)
        {
          const auto &diagnostic = artifact.diagnostics().diagnostics()[index];
          const auto  prefix = "diagnostic[" + number(index) + "]";
          write_line(output,
                     prefix + ".category",
                     diagnostic_category_name(diagnostic.category));
          write_line(output, prefix + ".component_id", diagnostic.component_id);
          write_line(output, prefix + ".capability", diagnostic.capability);
          write_line(output, prefix + ".remedy", diagnostic.remedy);
        }

      write_line(output,
                 "measurements.timing_collected",
                 artifact.measurements().timing_collected ? "true" : "false");
      if (artifact.measurements().timing_collected)
        {
          write_line(output,
                     "measurements.wall_seconds",
                     number(artifact.measurements().wall_seconds));
          write_line(output,
                     "measurements.cpu_seconds",
                     number(artifact.measurements().cpu_seconds));
        }
      write_line(output,
                 "measurements.memory_collected",
                 artifact.measurements().memory_collected ? "true" : "false");
      if (artifact.measurements().memory_collected)
        write_line(output,
                   "measurements.peak_memory_bytes",
                   number(artifact.measurements().peak_memory_bytes));

      write_line(output,
                 "selected_field_count",
                 number(artifact.selected_fields().size()));
      for (std::size_t index = 0; index < artifact.selected_fields().size(); ++index)
        write_line(output,
                   "selected_field[" + number(index) + "]",
                   artifact.selected_fields()[index]);

      for (const auto &field : fields)
        write_line(output, field.key, field.value);
      return output.str();
    }

    template <typename Envelope>
    void
    write(std::ostream &                 output,
          const BenchmarkArtifactT<Envelope> &artifact,
          std::vector<ArtifactField>     fields = {}) const
    {
      output << render(artifact, std::move(fields));
      if (!output)
        throw std::runtime_error("could not write benchmark artifact");
    }

  private:
    static void
    write_line(std::ostream &output,
               const std::string &key,
               const std::string &value)
    {
      output << key << '=' << escape_artifact_value(value) << '\n';
    }

    template <typename Value>
    static std::string
    number(const Value value)
    {
      std::ostringstream output;
      output.imbue(std::locale::classic());
      output << std::setprecision(17) << value;
      return output.str();
    }
  };
} // namespace nmopt::application::benchmark
