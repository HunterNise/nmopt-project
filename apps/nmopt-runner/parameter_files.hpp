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
#include <iterator>
#include <limits>
#include <map>
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

    inline void
    declare(::dealii::ParameterHandler &handler,
            const std::vector<std::string> &sections,
            const std::string              &entry,
            const std::string              &default_value = {})
    {
      for (const auto &section : sections)
        handler.enter_subsection(section);
      handler.declare_entry(entry,
                            default_value,
                            ::dealii::Patterns::Anything());
      for (std::size_t index = 0; index < sections.size(); ++index)
        handler.leave_subsection();
    }

    inline std::string
    get(::dealii::ParameterHandler &handler,
        const std::vector<std::string> &sections,
        const std::string              &entry)
    {
      for (const auto &section : sections)
        handler.enter_subsection(section);
      const auto value = handler.get(entry);
      for (std::size_t index = 0; index < sections.size(); ++index)
        handler.leave_subsection();
      return value;
    }

    inline void
    declare_schema(::dealii::ParameterHandler &handler)
    {
      const auto declare_section = [&](const std::string &section,
                                       const std::string &entry,
                                       const std::string &value = {}) {
        declare(handler, {section}, entry, value);
      };

      for (const auto &entry : {"id",
                                "label",
                                "recipe",
                                "source reference",
                                "source revision"})
        declare_section("Benchmark", entry);

      // ParameterHandler has no enumeration API. Matrix and selection keys
      // are therefore registered once here, while their values remain fully
      // data-driven and benchmark adapters only consume resolved combinations.
      for (const auto &entry : {"method",
                                "regularisation",
                                "case",
                                "forcing",
                                "observation-region",
                                "target-profile"})
        {
          declare_section("Matrix", entry);
          declare_section("Selection", entry);
        }
      declare_section("Selection", "exclude combinations");

      for (const auto &entry : {"control representation",
                                "cellwise box constraint",
                                "observation",
                                "observed material id",
                                "facewise box constraint",
                                "initial control"})
        declare_section("Problem", entry);

      declare_section("Observation", "active region");
      declare_section("Observation", "material id");
      for (const auto &region : {"wings", "full"})
        declare(handler, {"Observation", "region " + std::string(region)},
                "geometry");

      for (const auto &entry : {"forcing",
                                "desired state",
                                "fixed Dirichlet data",
                                "conservative transport"})
        declare_section("Functions", entry);
      for (const auto &entry : {"kind", "expression", "provenance", "value"})
        {
          declare(handler, {"Functions", "forcing"}, entry);
          declare(handler, {"Functions", "desired state"}, entry);
          declare(handler, {"Functions", "fixed-temperature"}, entry);
          declare(handler, {"Functions", "graetz"}, entry);
        }
      for (const auto &profile : {"constant", "parabolic"})
        for (const auto &entry : {"id", "kind", "expression", "provenance", "value"})
          declare(handler,
                  {"Functions", "target definitions", profile},
                  entry);
      for (const auto &forcing : {"zero", "constant-one", "constant-two"})
        for (const auto &entry : {"kind", "provenance", "value"})
          declare(handler,
                  {"Functions", "forcing definition " + std::string(forcing)},
                  entry);

      for (const auto &entry : {"diffusion", "reaction", "regularisation"})
        declare_section("Runtime", entry);

      for (const auto &entry : {"fixed region",
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
                                "face quadrature"})
        declare_section("Boundary", entry);

      for (const auto &entry : {"dimension",
                                "geometry",
                                "lower",
                                "upper",
                                "refinement",
                                "provenance"})
        declare_section("Mesh", entry);
      declare_section("Mesh", "generator", "framework-native");
      declare_section("Mesh", "subdivisions", "0");
      declare_section("Mesh", "axis subdivisions");
      declare_section("Mesh", "centroid splits", "0");
      declare_section("Mesh", "selection seed", "0");
      for (const auto &entry : {"state degree",
                                "execution",
                                "product",
                                "owned session",
                                "stabilization"})
        declare_section("Compile", entry);
      declare_section("Compile", "volume observation quadrature order");
      declare_section("Compile", "volume observation target realisation");
      declare_section("Compile", "state solve maximum iterations", "0");
      declare_section("Compile", "state solve relative tolerance", "1e-12");
      declare_section("Compile", "state solve absolute tolerance", "1e-14");
      declare_section("Compile", "adjoint solve maximum iterations", "0");
      declare_section("Compile", "adjoint solve relative tolerance", "1e-12");
      declare_section("Compile", "adjoint solve absolute tolerance", "1e-14");
      declare_section("Compile",
                      "control metric solve maximum iterations",
                      "1000");
      declare_section("Compile",
                      "control metric solve relative tolerance",
                      "1e-12");
      declare_section("Compile",
                      "control metric solve absolute tolerance",
                      "1e-14");

      declare_section("Solver", "globalization", "armijo");
      for (const auto &entry : {"method",
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
                                "declared minimum step length"})
        declare_section("Solver", entry);
      for (const auto &method : {"steepest-descent", "l-bfgs", "bfgs"})
        for (const auto *entry : method_policy_entries)
          declare(handler,
                  {"Solver", "method policy " + std::string(method)},
                  entry);

      for (const auto &entry : {"kind",
                                "build profile",
                                "output root",
                                "deterministic",
                                "serialize artifacts",
                                "measure timings",
                                "measure memory"})
        declare_section("Run", entry);
      for (const auto &entry : {"retain fields", "selected fields"})
        declare_section("Output", entry);
      for (const auto &entry : {"style profile",
                                "comparison rows",
                                "comparison columns",
                                "comparison group by",
                                "output formats"})
        declare_section("Postprocessing", entry);
    }

    inline std::map<std::string, std::string>
    read_values(::dealii::ParameterHandler &handler)
    {
      std::map<std::string, std::string> values;
      const auto remember = [&](const std::string &section,
                                const std::string &entry) {
        values.emplace(section + "/" + entry,
                       get(handler, {section}, entry));
      };

      for (const auto &entry : {"id",
                                "label",
                                "recipe",
                                "source reference",
                                "source revision"})
        remember("Benchmark", entry);
      for (const auto &entry : {"method",
                                "regularisation",
                                "case",
                                "forcing",
                                "observation-region",
                                "target-profile"})
        remember("Matrix", entry);
      for (const auto &entry : {"method",
                                "regularisation",
                                "case",
                                "forcing",
                                "observation-region",
                                "target-profile"})
        remember("Selection", entry);
      remember("Selection", "exclude combinations");
      for (const auto &entry : {"control representation",
                                "cellwise box constraint",
                                "observation",
                                "observed material id",
                                "facewise box constraint",
                                "initial control"})
        remember("Problem", entry);
      for (const auto &entry : {"active region", "material id"})
        remember("Observation", entry);
      for (const auto &region : {"wings", "full"})
        values.emplace("Observation/region " + std::string(region) + "/geometry",
                       get(handler,
                           {"Observation", "region " + std::string(region)},
                           "geometry"));
      for (const auto &entry : {"forcing",
                                "desired state",
                                "fixed Dirichlet data",
                                "conservative transport"})
        remember("Functions", entry);
      for (const auto &entry : {"kind", "expression", "provenance", "value"})
        {
          for (const auto &section : {"forcing",
                                      "desired state",
                                      "fixed-temperature",
                                      "graetz"})
            values.emplace("Functions/" + std::string(section) + "/" + entry,
                           get(handler, {"Functions", section}, entry));
        }
      for (const auto &profile : {"constant", "parabolic"})
        for (const auto &entry : {"id", "kind", "expression", "provenance", "value"})
          values.emplace("Functions/target definitions/" +
                           std::string(profile) + "/" + entry,
                         get(handler,
                             {"Functions", "target definitions", profile},
                             entry));
      for (const auto &forcing : {"zero", "constant-one", "constant-two"})
        for (const auto &entry : {"kind", "provenance", "value"})
          values.emplace("Functions/forcing definition " +
                           std::string(forcing) + "/" + entry,
                         get(handler,
                             {"Functions",
                              "forcing definition " + std::string(forcing)},
                             entry));

      for (const auto &section : {"Runtime", "Boundary", "Mesh", "Compile",
                                  "Solver", "Run", "Output", "Postprocessing"})
        {
          const std::vector<std::string> entries =
            section == "Runtime"
              ? std::vector<std::string>{"diffusion", "reaction", "regularisation"}
              : section == "Boundary"
                  ? std::vector<std::string>{"fixed region",
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
                                             "face quadrature"}
                  : section == "Mesh"
                      ? std::vector<std::string>{"dimension",
                                                 "geometry",
                                                 "lower",
                                                 "upper",
                                                 "refinement",
                                                 "generator",
                                                 "subdivisions",
                                                 "axis subdivisions",
                                                 "centroid splits",
                                                 "selection seed",
                                                 "provenance"}
                      : section == "Compile"
                          ? std::vector<std::string>{"state degree",
                                                     "execution",
                                                     "product",
                                                     "owned session",
                                                     "stabilization",
                                                     "volume observation quadrature order",
                                                     "volume observation target realisation",
                                                     "state solve maximum iterations",
                                                     "state solve relative tolerance",
                                                     "state solve absolute tolerance",
                                                     "adjoint solve maximum iterations",
                                                     "adjoint solve relative tolerance",
                                                     "adjoint solve absolute tolerance",
                                                     "control metric solve maximum iterations",
                                                     "control metric solve relative tolerance",
                                                     "control metric solve absolute tolerance"}
                          : section == "Solver"
                              ? std::vector<std::string>{"method",
                                                         "globalization",
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
                                                         "declared minimum step length"}
                              : section == "Run"
                                  ? std::vector<std::string>{"kind",
                                                             "build profile",
                                                             "output root",
                                                             "deterministic",
                                                             "serialize artifacts",
                                                             "measure timings",
                                                             "measure memory"}
                                  : section == "Output"
                                      ? std::vector<std::string>{"retain fields",
                                                                 "selected fields"}
                                      : std::vector<std::string>{"style profile",
                                                                 "comparison rows",
                                                                 "comparison columns",
                                                                 "comparison group by",
                                                                 "output formats"};
          for (const auto &entry : entries)
            remember(section, entry);
        }
      for (const auto &method : {"steepest-descent", "l-bfgs", "bfgs"})
        for (const auto *entry : method_policy_entries)
          values.emplace("Solver/method policy " + std::string(method) + "/" +
                           entry,
                         get(handler,
                             {"Solver",
                              "method policy " + std::string(method)},
                             entry));
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
    for (const auto &axis_id : {"method",
                                "regularisation",
                                "case",
                                "forcing",
                                "observation-region",
                                "target-profile"})
      {
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
