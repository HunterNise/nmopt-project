#pragma once

#include "nmopt/application/chapter6.hpp"
#include "nmopt/semantic/v1/types.hpp"

#include <deal.II/base/parameter_handler.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace nmopt::application::runner
{
  struct ParameterAxis
  {
    std::string              id;
    std::vector<std::string> values;
  };

  struct ParameterCombination
  {
    std::map<std::string, std::string> values;
  };

  struct ParameterFile
  {
    std::filesystem::path              path;
    std::string                        content_hash;
    std::map<std::string, std::string> values;
    std::vector<ParameterAxis>         matrix;
    std::map<std::string, std::string> selection;
    std::vector<ParameterCombination>  excluded_combinations;

    const std::string &
    value(const std::string_view key) const
    {
      const auto found = values.find(std::string(key));
      if (found == values.end())
        throw std::invalid_argument("parameter file has no entry '" +
                                    std::string(key) + "'");
      return found->second;
    }

    std::string
    optional_value(const std::string_view key,
                   const std::string_view fallback = {}) const
    {
      const auto found = values.find(std::string(key));
      return found == values.end() ? std::string(fallback) : found->second;
    }

    std::vector<ParameterCombination>
    combinations(const std::vector<std::pair<std::string, std::string>> &cli_filters = {}) const;
  };

  struct ResolvedParameterValue
  {
    std::string key;
    std::string value;
  };

  enum class ParameterPresence
  {
    optional,
    required
  };

  enum class ParameterOwnership
  {
    consumed,
    locked_profile,
    provenance_only
  };

  struct ParameterSchemaEntry
  {
    std::string path;
    std::string default_value;
    std::shared_ptr<const ::dealii::Patterns::PatternBase> pattern;
    ParameterPresence presence;
    ParameterOwnership ownership;
  };

  inline semantic::v1::TransportBoundaryForm
  b2_transport_boundary_form(const ParameterFile &file)
  {
    const auto &boundary_form =
      file.value("Boundary/transport boundary form");
    const auto &conormal_form = file.value("Boundary/conormal form");
    if (boundary_form == "ordinary-normal-minus-transport")
      {
        if (conormal_form != "unspecified")
          throw std::invalid_argument(
            "B2 ordinary-normal-minus-transport requires an unspecified conormal form");
        return semantic::v1::TransportBoundaryForm::
          ordinary_normal_minus_transport;
      }
    if (boundary_form == "total-conormal")
      {
        if (conormal_form != "diffusion-minus-transport")
          throw std::invalid_argument(
            "B2 total-conormal requires the diffusion-minus-transport conormal form");
        return semantic::v1::TransportBoundaryForm::total_conormal;
      }
    throw std::invalid_argument("unknown B2 transport boundary form '" +
                                boundary_form + "'");
  }

  inline semantic::v1::NeumannControlDiscretisation
  b2_neumann_control_discretisation(const ParameterFile &file)
  {
    const auto &representation = file.value("Problem/control representation");
    if (representation == "facewise-constant")
      return semantic::v1::NeumannControlDiscretisation::facewise_constant;
    if (representation == "continuous-nodal-trace")
      return semantic::v1::NeumannControlDiscretisation::continuous_nodal_trace;
    throw std::invalid_argument("unknown B2 control representation '" +
                                representation + "'");
  }

  inline chapter6::ReducedGlobalization
  reduced_globalization(const ParameterFile &file)
  {
    const auto selection = file.optional_value("Solver/globalization", "armijo");
    if (selection.empty() || selection == "armijo")
      return chapter6::ReducedGlobalization::armijo;
    if (selection == "fixed-step")
      return chapter6::ReducedGlobalization::fixed_step;
    throw std::invalid_argument("unknown reduced globalization '" + selection +
                                "'");
  }

  inline ResolvedParameterValue
  resolve_method_parameter(const ParameterFile &file,
                           const std::string_view method_id,
                           const std::string_view entry)
  {
    const auto policy_key = "Solver/method policy " + std::string(method_id) +
                            "/" + std::string(entry);
    const auto policy_value = file.optional_value(policy_key);
    if (!policy_value.empty())
      return {policy_key, policy_value};
    const auto global_key = "Solver/" + std::string(entry);
    return {global_key, file.optional_value(global_key)};
  }

  namespace detail
  {
    inline constexpr std::array<const char *, 16> method_policy_entries = {
      {"maximum iterations",
       "maximum line search trials",
       "maximum backtracking reductions",
       "gradient tolerance",
       "stopping criterion",
       "relative gradient tolerance",
       "objective change tolerance",
       "step tolerance",
       "initial step length",
       "Armijo fraction",
       "backtracking factor",
       "minimum step length",
       "memory",
       "curvature tolerance",
       "initial inverse Hessian scaling",
       "declared minimum step length"}};

    inline std::string hash_text(const std::string &text);
  }

  inline std::string
  parameter_file_text(const std::filesystem::path &path)
  {
    std::ifstream input(path);
    if (!input)
      throw std::invalid_argument("could not open configuration file '" +
                                  path.string() + "'");
    return {std::istreambuf_iterator<char>(input),
            std::istreambuf_iterator<char>()};
  }

  inline std::string
  parameter_file_hash(const std::filesystem::path &path)
  {
    return detail::hash_text(parameter_file_text(path));
  }

  inline std::filesystem::path
  find_file_from_current_or_parent(const std::filesystem::path &path)
  {
    if (std::filesystem::exists(path))
      return path;
    auto directory = std::filesystem::current_path();
    for (;;)
      {
        const auto candidate = directory / path;
        if (std::filesystem::exists(candidate))
          return candidate;
        const auto parent = directory.parent_path();
        if (parent == directory)
          break;
        directory = parent;
      }
    throw std::invalid_argument("could not locate configuration file '" +
                                path.string() + "'");
  }

  inline double
  parameter_double(const ParameterFile &file, const std::string_view key)
  {
    const auto text = file.value(key);
    std::size_t consumed = 0;
    double      result = 0.0;
    try
      {
        result = std::stod(text, &consumed);
      }
    catch (const std::exception &)
      {
        throw std::invalid_argument("parameter '" + std::string(key) +
                                    "' needs a finite number");
      }
    if (consumed != text.size())
      throw std::invalid_argument("parameter '" + std::string(key) +
                                  "' needs a finite number");
    return result;
  }

  inline unsigned int
  parameter_unsigned(const ParameterFile &file, const std::string_view key)
  {
    const auto value = parameter_double(file, key);
    if (value < 0.0 || value != static_cast<unsigned int>(value))
      throw std::invalid_argument("parameter '" + std::string(key) +
                                  "' needs a nonnegative integer");
    return static_cast<unsigned int>(value);
  }

  inline chapter6::VolumeObservationOptions
  b2_volume_observation_options(const ParameterFile &file)
  {
    chapter6::VolumeObservationOptions options;
    options.quadrature_order = parameter_unsigned(
      file, "Compile/volume observation quadrature order");
    if (options.quadrature_order == 0)
      throw std::invalid_argument(
        "B2 volume-observation quadrature order must be positive");
    const auto &target_realisation =
      file.value("Compile/volume observation target realisation");
    if (target_realisation == "analytic-quadrature")
      options.target_realisation = chapter6::
        VolumeObservationTargetRealisation::analytic_quadrature;
    else if (target_realisation == "state-fe-interpolation")
      options.target_realisation = chapter6::
        VolumeObservationTargetRealisation::state_fe_interpolation;
    else
      throw std::invalid_argument(
        "unknown B2 volume-observation target realisation '" +
        target_realisation + "'");
    return options;
  }

  inline bool
  parameter_bool(const ParameterFile &file, const std::string_view key)
  {
    const auto value = file.value(key);
    if (value == "true")
      return true;
    if (value == "false")
      return false;
    throw std::invalid_argument("parameter '" + std::string(key) +
                                "' needs true or false");
  }

  inline ScalarFunctionDefinition
  parameter_scalar_function_definition_at(const ParameterFile &file,
                                          const std::string &path,
                                          const std::string &id_path)
  {
    ScalarFunctionDefinition definition;
    definition.id = file.value(id_path);
    definition.kind =
      scalar_function_kind_from_name(file.value(path + "/kind"));
    definition.provenance = file.value(path + "/provenance");
    const auto value = file.optional_value(path + "/value");
    const auto expression = file.optional_value(path + "/expression");
    if (definition.kind == ScalarFunctionKind::constant)
      {
        if (value.empty())
          throw std::invalid_argument("parameter '" + path +
                                      "/value' needs a finite number");
        if (!expression.empty())
          throw std::invalid_argument("parameter '" + path +
                                      "' cannot set both value and expression");
        definition.value = parameter_double(file, path + "/value");
      }
    else if (definition.kind == ScalarFunctionKind::expression)
      {
        if (!value.empty())
          throw std::invalid_argument("parameter '" + path +
                                      "' cannot set both value and expression");
        definition.expression = expression;
      }
    else if (!value.empty() || !expression.empty())
      throw std::invalid_argument("parameter '" + path +
                                  "' zero kind cannot set value or expression");
    validate_scalar_function_definition(definition, path);
    return definition;
  }

  inline ScalarFunctionDefinition
  parameter_scalar_function_definition(const ParameterFile &file,
                                       const std::string_view key)
  {
    const auto path = std::string(key);
    return parameter_scalar_function_definition_at(file, path, path);
  }

  inline ScalarFunctionDefinition
  parameter_scalar_function_section_definition(
    const ParameterFile &file,
    const std::string_view section)
  {
    const auto path = "Functions/target definitions/" + std::string(section);
    return parameter_scalar_function_definition_at(file, path, path + "/id");
  }

  inline chapter6::B2TargetParameters
  b2_target_parameters(const ParameterFile &file)
  {
    chapter6::B2TargetParameters parameters;
    parameters.constant =
      parameter_scalar_function_section_definition(file, "constant");
    parameters.parabolic =
      parameter_scalar_function_section_definition(file, "parabolic");
    chapter6::validate_b2_target_parameters(parameters);
    return parameters;
  }

  namespace detail
  {
    inline std::string
    trim(std::string value)
    {
      const auto first = value.find_first_not_of(" \t\r\n");
      if (first == std::string::npos)
        return {};
      const auto last = value.find_last_not_of(" \t\r\n");
      return value.substr(first, last - first + 1);
    }

    inline std::vector<std::string>
    split_list(const std::string &value)
    {
      std::vector<std::string> result;
      std::size_t              begin = 0;
      while (begin <= value.size())
        {
          const auto end = value.find(',', begin);
          const auto item = trim(value.substr(
            begin, end == std::string::npos ? std::string::npos : end - begin));
          if (!item.empty())
            result.push_back(item);
          if (end == std::string::npos)
            break;
          begin = end + 1;
        }
      return result;
    }

    inline std::vector<ParameterCombination>
    parse_excluded_combinations(const std::string                &value,
                                const std::vector<ParameterAxis> &matrix)
    {
      if (trim(value).empty())
        return {};

      std::vector<ParameterCombination> result;
      std::size_t                       begin = 0;
      while (begin <= value.size())
        {
          const auto end = value.find(';', begin);
          const auto text = trim(value.substr(
            begin, end == std::string::npos ? std::string::npos : end - begin));
          if (text.size() < 2 || text.front() != '[' || text.back() != ']')
            throw std::invalid_argument("excluded matrix combination '" +
                                        text +
                                        "' must be enclosed in brackets");

          ParameterCombination combination;
          const auto           assignments =
            text.substr(1, text.size() - 2);
          std::size_t assignment_begin = 0;
          while (assignment_begin <= assignments.size())
            {
              const auto assignment_end =
                assignments.find(',', assignment_begin);
              const auto assignment = trim(assignments.substr(
                assignment_begin,
                assignment_end == std::string::npos ?
                  std::string::npos :
                  assignment_end - assignment_begin));
              const auto separator = assignment.find('=');
              if (separator == std::string::npos)
                throw std::invalid_argument(
                  "excluded matrix combination '" + text +
                  "' must contain axis=value assignments");
              const auto axis = trim(assignment.substr(0, separator));
              const auto selected = trim(assignment.substr(separator + 1));
              if (axis.empty() || selected.empty())
                throw std::invalid_argument(
                  "excluded matrix combination '" + text +
                  "' contains an empty axis or value");
              if (!combination.values.emplace(axis, selected).second)
                throw std::invalid_argument(
                  "excluded matrix combination '" + text +
                  "' repeats axis '" + axis + "'");
              if (assignment_end == std::string::npos)
                break;
              assignment_begin = assignment_end + 1;
            }

          for (const auto &[axis, selected] : combination.values)
            {
              const auto declared = std::find_if(
                matrix.begin(), matrix.end(), [&](const auto &candidate) {
                  return candidate.id == axis;
                });
              if (declared == matrix.end())
                throw std::invalid_argument(
                  "excluded matrix combination '" + text +
                  "' names undeclared matrix axis '" + axis + "'");
              if (std::find(declared->values.begin(),
                            declared->values.end(),
                            selected) == declared->values.end())
                throw std::invalid_argument(
                  "excluded matrix combination '" + text +
                  "' selects unknown value '" + selected +
                  "' on matrix axis '" + axis + "'");
            }
          for (const auto &axis : matrix)
            if (combination.values.find(axis.id) == combination.values.end())
              throw std::invalid_argument(
                "excluded matrix combination '" + text +
                "' is missing matrix axis '" + axis.id + "'");
          if (std::find_if(result.begin(), result.end(), [&](const auto &entry) {
                return entry.values == combination.values;
              }) != result.end())
            throw std::invalid_argument("excluded matrix combination '" +
                                        text + "' is duplicated");
          result.push_back(std::move(combination));

          if (end == std::string::npos)
            break;
          begin = end + 1;
        }
      return result;
    }

    inline std::string
    hash_text(const std::string &text)
    {
      // A deterministic digest is enough to detect configuration drift. It is
      // deliberately labelled because this is provenance, not authentication.
      std::uint64_t hash = 1469598103934665603ULL;
      for (const unsigned char character : text)
        {
          hash ^= character;
          hash *= 1099511628211ULL;
        }
      std::ostringstream output;
      output << "fnv1a64:" << std::hex << std::setfill('0') << std::setw(16)
             << hash;
      return output.str();
    }

    class EmptyOrPattern final : public ::dealii::Patterns::PatternBase
    {
    public:
      explicit EmptyOrPattern(
        std::shared_ptr<const ::dealii::Patterns::PatternBase> pattern)
        : pattern(std::move(pattern))
      {}

      bool
      match(const std::string &test_string) const override
      {
        return test_string.empty() || pattern->match(test_string);
      }

      std::string
      description(const OutputStyle style = Machine) const override
      {
        return "[empty or " + pattern->description(style) + "]";
      }

      std::unique_ptr<PatternBase>
      clone() const override
      {
        return std::make_unique<EmptyOrPattern>(*this);
      }

    private:
      std::shared_ptr<const ::dealii::Patterns::PatternBase> pattern;
    };

    inline std::shared_ptr<const ::dealii::Patterns::PatternBase>
    anything_pattern()
    {
      return std::make_shared<::dealii::Patterns::Anything>();
    }

    inline std::shared_ptr<const ::dealii::Patterns::PatternBase>
    empty_or(std::shared_ptr<const ::dealii::Patterns::PatternBase> pattern)
    {
      return std::make_shared<EmptyOrPattern>(std::move(pattern));
    }

    inline bool
    ends_with(const std::string &value, const std::string &suffix)
    {
      return value.size() >= suffix.size() &&
             value.compare(value.size() - suffix.size(), suffix.size(), suffix) ==
               0;
    }

    inline std::shared_ptr<const ::dealii::Patterns::PatternBase>
    pattern_for(const std::string &path)
    {
      const auto entry = path.substr(path.rfind('/') + 1);

      if (path.rfind("Matrix/", 0) == 0 ||
          path.rfind("Selection/", 0) == 0)
        return anything_pattern();

      if (path == "Solver/globalization")
        return std::make_shared<::dealii::Patterns::Selection>(
          "armijo|fixed-step");

      if (path == "Compile/volume observation target realisation")
        return std::make_shared<::dealii::Patterns::MultipleSelection>(
          "analytic-quadrature|state-fe-interpolation");

      if (path == "Problem/control representation")
        return std::make_shared<::dealii::Patterns::MultipleSelection>(
          "cellwise-volume|continuous-volume-homogeneous-dirichlet|"
          "facewise-constant|continuous-nodal-trace");

      if (entry == "cellwise box constraint" ||
          entry == "facewise box constraint" || entry == "owned session" ||
          entry == "deterministic" ||
          entry == "serialize artifacts" || entry == "measure timings" ||
          entry == "measure memory" || entry == "retain fields")
        return empty_or(std::make_shared<::dealii::Patterns::Bool>());

      if (entry == "dimension" || entry == "refinement" ||
          entry == "subdivisions" || entry == "centroid splits" ||
          entry == "selection seed" || entry == "state degree" ||
          entry == "volume observation quadrature order" ||
          entry == "observed material id" || entry == "material id" ||
          entry == "fixed boundary id" || entry == "control boundary id" ||
          entry == "outflow boundary id" || entry == "maximum iterations" ||
          entry == "maximum line search trials" ||
          entry == "maximum backtracking reductions" || entry == "memory")
        return empty_or(std::make_shared<::dealii::Patterns::Integer>());

      if (entry == "diffusion" || entry == "reaction" ||
          entry == "regularisation" || entry == "upstream transition x" ||
          entry == "outflow x" || entry == "value" ||
          ends_with(entry, "tolerance") || entry == "initial step length" ||
          entry == "Armijo fraction" || entry == "backtracking factor" ||
          entry == "minimum step length")
        return empty_or(std::make_shared<::dealii::Patterns::Double>());

      return anything_pattern();
    }

    inline ParameterOwnership
    ownership_for(const std::string &path)
    {
      if (path == "Benchmark/source reference" ||
          path == "Benchmark/source revision" ||
          path == "Solver/declared minimum step length" ||
          ends_with(path, "/provenance"))
        return ParameterOwnership::provenance_only;

      if (path.rfind("Boundary/", 0) == 0 || path == "Mesh/geometry" ||
          path == "Mesh/generator" ||
          path == "Functions/conservative transport" ||
          path.rfind("Functions/graetz/", 0) == 0)
        return ParameterOwnership::locked_profile;

      return ParameterOwnership::consumed;
    }

    inline void
    append_schema_entry(std::vector<ParameterSchemaEntry> &registry,
                        const std::string                      &path,
                        const std::string                      &default_value = {},
                        const ParameterPresence presence =
                          ParameterPresence::optional)
    {
      if (std::any_of(registry.begin(), registry.end(), [&](const auto &entry) {
            return entry.path == path;
          }))
        throw std::logic_error("duplicate parameter schema path '" + path +
                              "'");
      registry.push_back({path,
                          default_value,
                          pattern_for(path),
                          presence,
                          ownership_for(path)});
    }

    inline const std::vector<ParameterSchemaEntry> &
    parameter_schema_registry()
    {
      static const auto registry = [] {
        std::vector<ParameterSchemaEntry> result;
        const auto add_section = [&](const std::string &section,
                                     const std::initializer_list<const char *> entries,
                                     const ParameterPresence presence =
                                       ParameterPresence::optional) {
          for (const auto *entry : entries)
            append_schema_entry(result, section + "/" + entry, {}, presence);
        };

        append_schema_entry(result,
                            "Benchmark/id",
                            {},
                            ParameterPresence::required);
        add_section("Benchmark", {"label"});
        append_schema_entry(result,
                            "Benchmark/recipe",
                            {},
                            ParameterPresence::required);
        add_section("Benchmark", {"source reference", "source revision"});
        add_section("Matrix",
                    {"method",
                     "regularisation",
                     "case",
                     "forcing",
                     "observation-region",
                     "target-profile"});
        add_section("Selection",
                    {"method",
                     "regularisation",
                     "case",
                     "forcing",
                     "observation-region",
                     "target-profile",
                     "exclude combinations"});
        add_section("Problem",
                    {"control representation",
                     "cellwise box constraint",
                     "observation",
                     "observed material id",
                     "facewise box constraint",
                     "initial control"});
        add_section("Observation", {"active region", "material id"});
        for (const auto *region : {"wings", "full"})
          append_schema_entry(result,
                              "Observation/region " + std::string(region) +
                                "/geometry");

        add_section("Functions",
                    {"forcing",
                     "desired state",
                     "fixed Dirichlet data",
                     "conservative transport"});
        for (const auto *section : {"forcing",
                                    "desired state",
                                    "fixed-temperature",
                                    "graetz"})
          for (const auto *entry : {"kind", "expression", "provenance", "value"})
            append_schema_entry(result,
                                "Functions/" + std::string(section) + "/" +
                                  entry);
        for (const auto *profile : {"constant", "parabolic"})
          for (const auto *entry : {"id", "kind", "expression", "provenance", "value"})
            append_schema_entry(result,
                                "Functions/target definitions/" +
                                  std::string(profile) + "/" + entry);
        for (const auto *forcing : {"zero", "constant-one", "constant-two"})
          for (const auto *entry : {"kind", "provenance", "value"})
            append_schema_entry(result,
                                "Functions/forcing definition " +
                                  std::string(forcing) + "/" + entry);

        add_section("Runtime", {"diffusion", "reaction", "regularisation"});
        add_section("Boundary",
                    {"fixed region",
                     "fixed boundary id",
                     "control region",
                     "control boundary id",
                     "outflow region",
                     "outflow boundary id",
                     "upstream transition x",
                     "outflow x",
                     "transport boundary form",
                     "conormal form",
                     "normal orientation",
                     "trace evaluation",
                     "face quadrature"});
        add_section("Mesh",
                    {"dimension",
                     "geometry",
                     "lower",
                     "upper",
                     "refinement",
                     "provenance"});
        append_schema_entry(result, "Mesh/generator", "framework-native");
        append_schema_entry(result, "Mesh/subdivisions", "0");
        append_schema_entry(result, "Mesh/axis subdivisions");
        append_schema_entry(result, "Mesh/centroid splits", "0");
        append_schema_entry(result, "Mesh/selection seed", "0");

        add_section("Compile",
                    {"state degree",
                     "execution",
                     "product",
                     "owned session",
                     "stabilization"});
        append_schema_entry(result, "Compile/volume observation quadrature order");
        append_schema_entry(result, "Compile/volume observation target realisation");
        append_schema_entry(result, "Compile/state solve maximum iterations", "0");
        append_schema_entry(result, "Compile/state solve relative tolerance", "1e-12");
        append_schema_entry(result, "Compile/state solve absolute tolerance", "1e-14");
        append_schema_entry(result, "Compile/adjoint solve maximum iterations", "0");
        append_schema_entry(result, "Compile/adjoint solve relative tolerance", "1e-12");
        append_schema_entry(result, "Compile/adjoint solve absolute tolerance", "1e-14");
        append_schema_entry(result,
                            "Compile/control metric solve maximum iterations",
                            "1000");
        append_schema_entry(result,
                            "Compile/control metric solve relative tolerance",
                            "1e-12");
        append_schema_entry(result,
                            "Compile/control metric solve absolute tolerance",
                            "1e-14");

        append_schema_entry(result, "Solver/globalization", "armijo");
        add_section("Solver",
                    {"method",
                     "initial control",
                     "maximum iterations",
                     "maximum line search trials",
                     "maximum backtracking reductions",
                     "gradient tolerance",
                     "stopping criterion",
                     "relative gradient tolerance",
                     "objective change tolerance",
                     "step tolerance",
                     "objective target",
                     "objective target policy",
                     "objective target reference method",
                     "initial step length",
                     "Armijo fraction",
                     "backtracking factor",
                     "minimum step length",
                     "declared minimum step length"});
        for (const auto *method : {"steepest-descent", "l-bfgs", "bfgs"})
          for (const auto *entry : method_policy_entries)
            append_schema_entry(result,
                                "Solver/method policy " + std::string(method) +
                                  "/" + entry);

        add_section("Run",
                    {"kind",
                     "build profile",
                     "output root",
                     "deterministic",
                     "serialize artifacts",
                     "measure timings",
                     "measure memory"});
        add_section("Output", {"retain fields", "selected fields"});
        add_section("Postprocessing",
                    {"style profile",
                     "comparison rows",
                     "comparison columns",
                     "comparison group by",
                     "output formats"});
        return result;
      }();
      return registry;
    }

    inline std::vector<std::string>
    schema_sections(const std::string &path)
    {
      const auto last_separator = path.rfind('/');
      if (last_separator == std::string::npos || last_separator == 0 ||
          last_separator + 1 == path.size())
        throw std::logic_error("parameter schema path must name a section and entry: '" +
                              path + "'");

      std::vector<std::string> sections;
      std::size_t              begin = 0;
      while (begin < last_separator)
        {
          const auto separator = path.find('/', begin);
          const auto end = separator == std::string::npos ||
                                   separator > last_separator ?
                                 last_separator :
                                 separator;
          sections.push_back(path.substr(begin, end - begin));
          begin = end + 1;
        }
      return sections;
    }

    inline std::string
    schema_entry_name(const std::string &path)
    {
      return path.substr(path.rfind('/') + 1);
    }

    inline void
    declare(::dealii::ParameterHandler &handler,
            const ParameterSchemaEntry   &schema_entry)
    {
      const auto sections = schema_sections(schema_entry.path);
      for (const auto &section : sections)
        handler.enter_subsection(section);
      handler.declare_entry(schema_entry_name(schema_entry.path),
                            schema_entry.default_value,
                            *schema_entry.pattern);
      for (std::size_t index = 0; index < sections.size(); ++index)
        handler.leave_subsection();
    }

    inline std::string
    get(::dealii::ParameterHandler &handler,
        const ParameterSchemaEntry   &schema_entry)
    {
      const auto sections = schema_sections(schema_entry.path);
      for (const auto &section : sections)
        handler.enter_subsection(section);
      const auto value = handler.get(schema_entry_name(schema_entry.path));
      for (std::size_t index = 0; index < sections.size(); ++index)
        handler.leave_subsection();
      return value;
    }

    inline void
    declare_schema(::dealii::ParameterHandler &handler)
    {
      for (const auto &schema_entry : parameter_schema_registry())
        declare(handler, schema_entry);
    }

    inline std::map<std::string, std::string>
    read_values(::dealii::ParameterHandler &handler)
    {
      std::map<std::string, std::string> values;
      for (const auto &schema_entry : parameter_schema_registry())
        values.emplace(schema_entry.path, get(handler, schema_entry));
      return values;
    }

    inline bool
    matches(const std::string &value, const std::string &filter)
    {
      const auto candidates = split_list(filter);
      return std::find(candidates.begin(), candidates.end(), value) !=
             candidates.end();
    }
  } // namespace detail

  inline std::vector<unsigned int>
  parameter_positive_unsigned_list(const ParameterFile &file,
                                   const std::string_view key)
  {
    const auto text = detail::trim(file.value(key));
    if (text.empty())
      return {};

    std::vector<unsigned int> result;
    std::size_t               begin = 0;
    while (begin <= text.size())
      {
        const auto end = text.find(',', begin);
        const auto item = detail::trim(text.substr(
          begin, end == std::string::npos ? std::string::npos : end - begin));
        if (item.empty())
          throw std::invalid_argument(
            "parameter '" + std::string(key) +
            "' needs comma-separated positive integers");

        std::size_t        consumed = 0;
        unsigned long long parsed = 0;
        try
          {
            parsed = std::stoull(item, &consumed);
          }
        catch (const std::exception &)
          {
            throw std::invalid_argument(
              "parameter '" + std::string(key) +
              "' needs comma-separated positive integers");
          }
        if (item.front() == '-' || consumed != item.size() || parsed == 0 ||
            parsed > std::numeric_limits<unsigned int>::max())
          throw std::invalid_argument(
            "parameter '" + std::string(key) +
            "' needs comma-separated positive integers");
        result.push_back(static_cast<unsigned int>(parsed));

        if (end == std::string::npos)
          break;
        begin = end + 1;
      }
    return result;
  }

  inline ParameterFile
  read_parameter_file(const std::filesystem::path &path)
  {
    std::ifstream input(path);
    if (!input)
      throw std::invalid_argument("could not open parameter file '" +
                                  path.string() + "'");
    const std::string content((std::istreambuf_iterator<char>(input)),
                              std::istreambuf_iterator<char>());

    ::dealii::ParameterHandler handler;
    detail::declare_schema(handler);
    std::istringstream input_stream(content);
    handler.parse_input(input_stream, path.string());

    ParameterFile result;
    result.path = path;
    result.content_hash = detail::hash_text(content);
    result.values = detail::read_values(handler);
    for (const auto &schema_entry : detail::parameter_schema_registry())
      if (schema_entry.path.rfind("Matrix/", 0) == 0)
      {
        const auto axis_id = schema_entry.path.substr(std::string("Matrix/").size());
        const auto values = detail::split_list(result.value(
          std::string("Matrix/") + axis_id));
        if (!values.empty())
          {
            for (std::size_t index = 0; index < values.size(); ++index)
              if (std::find(values.begin() + static_cast<std::ptrdiff_t>(index + 1),
                            values.end(),
                            values[index]) != values.end())
                throw std::invalid_argument("matrix axis '" + std::string(axis_id) +
                                            "' contains duplicate value '" +
                                            values[index] + "'");
            result.matrix.push_back({axis_id, values});
          }
      }
    for (const auto &axis : result.matrix)
      {
        const auto selected = detail::split_list(result.value(
          std::string("Selection/") + axis.id));
        if (!selected.empty())
          result.selection.emplace(axis.id, result.value(
            std::string("Selection/") + axis.id));
      }
    if (result.matrix.empty())
      throw std::invalid_argument("parameter file must declare at least one matrix axis");
    result.excluded_combinations = detail::parse_excluded_combinations(
      result.value("Selection/exclude combinations"), result.matrix);
    if (result.value("Benchmark/id").empty() ||
        result.value("Benchmark/recipe").empty())
      throw std::invalid_argument("parameter file needs Benchmark id and recipe");
    return result;
  }

  inline std::vector<ParameterCombination>
  ParameterFile::combinations(
    const std::vector<std::pair<std::string, std::string>> &cli_filters) const
  {
    std::map<std::string, std::string> filters = selection;
    for (const auto &filter : cli_filters)
      {
        const auto axis = std::find_if(matrix.begin(), matrix.end(), [&](const auto &candidate) {
          return candidate.id == filter.first;
        });
        if (axis == matrix.end())
          throw std::invalid_argument("selection names undeclared matrix axis '" +
                                      filter.first + "'");
        filters[filter.first] = filter.second;
      }

    std::vector<ParameterCombination> result(1);
    for (const auto &axis : matrix)
      {
        std::vector<std::string> allowed = axis.values;
        const auto filter = filters.find(axis.id);
          if (filter != filters.end())
          {
            allowed = detail::split_list(filter->second);
            if (allowed.empty())
              throw std::invalid_argument("selection for matrix axis '" +
                                          axis.id + "' is empty");
            for (const auto &value : allowed)
              if (std::find(axis.values.begin(), axis.values.end(), value) ==
                  axis.values.end())
                throw std::invalid_argument("selection value '" + value +
                                            "' is not declared on matrix axis '" +
                                            axis.id + "'");
            for (std::size_t index = 0; index < allowed.size(); ++index)
              if (std::find(allowed.begin() +
                              static_cast<std::ptrdiff_t>(index + 1),
                            allowed.end(),
                            allowed[index]) != allowed.end())
                throw std::invalid_argument("selection for matrix axis '" +
                                            axis.id + "' contains duplicate value '" +
                                            allowed[index] + "'");
          }

        std::vector<ParameterCombination> expanded;
        for (const auto &partial : result)
          for (const auto &value : allowed)
            {
              auto combination = partial;
              combination.values[axis.id] = value;
              expanded.push_back(std::move(combination));
            }
        result = std::move(expanded);
      }
    result.erase(
      std::remove_if(
        result.begin(),
        result.end(),
        [&](const auto &combination) {
          return std::find_if(
                   excluded_combinations.begin(),
                   excluded_combinations.end(),
                   [&](const auto &excluded) {
                     return excluded.values == combination.values;
                   }) != excluded_combinations.end();
        }),
      result.end());
    if (result.empty())
      throw std::invalid_argument("parameter selection resolves to no combinations");
    return result;
  }
} // namespace nmopt::application::runner
