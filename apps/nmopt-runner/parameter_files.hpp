#pragma once

#include <deal.II/base/parameter_handler.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
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

  namespace detail
  {
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
      for (const auto &entry : {"constant", "parabolic", "provenance"})
        declare(handler, {"Functions", "target definitions"}, entry);
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
      for (const auto &entry : {"state degree",
                                "execution",
                                "product",
                                "owned session",
                                "stabilization"})
        declare_section("Compile", entry);

      for (const auto &entry : {"method",
                                "initial control",
                                "maximum iterations",
                                "maximum line search trials",
                                "gradient tolerance",
                                "stopping criterion",
                                "relative gradient tolerance",
                                "objective change tolerance",
                                "step tolerance",
                                "initial step length",
                                "Armijo fraction",
                                "backtracking factor",
                                "declared minimum step length"})
        declare_section("Solver", entry);
      for (const auto &method : {"steepest-descent", "l-bfgs", "bfgs"})
        for (const auto &entry : {"gradient tolerance",
                                  "declared minimum step length"})
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
      for (const auto &entry : {"constant", "parabolic", "provenance"})
        values.emplace("Functions/target definitions/" + std::string(entry),
                       get(handler, {"Functions", "target definitions"}, entry));
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
                                                 "provenance"}
                      : section == "Compile"
                          ? std::vector<std::string>{"state degree",
                                                     "execution",
                                                     "product",
                                                     "owned session",
                                                     "stabilization"}
                          : section == "Solver"
                              ? std::vector<std::string>{"method",
                                                         "initial control",
                                                         "maximum iterations",
                                                         "maximum line search trials",
                                                         "gradient tolerance",
                                                         "stopping criterion",
                                                         "relative gradient tolerance",
                                                         "objective change tolerance",
                                                         "step tolerance",
                                                         "initial step length",
                                                         "Armijo fraction",
                                                         "backtracking factor",
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
        for (const auto &entry : {"gradient tolerance",
                                  "declared minimum step length"})
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
    if (result.empty())
      throw std::invalid_argument("parameter selection resolves to no combinations");
    return result;
  }
} // namespace nmopt::application::runner
