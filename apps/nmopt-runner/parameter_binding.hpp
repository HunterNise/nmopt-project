#pragma once

#include "nmopt/application/chapter6.hpp"
#include "parameter_files.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

namespace nmopt::application::runner::binding
{
  inline const std::string &
  combination_value(const ParameterCombination &combination,
                    const std::string_view   axis)
  {
    const auto found = combination.values.find(std::string(axis));
    if (found == combination.values.end())
      throw std::invalid_argument("resolved combination has no matrix axis '" +
                                  std::string(axis) + "'");
    return found->second;
  }

  inline nmopt::application::chapter6::ReducedMethod
  parse_method(const std::string &value)
  {
    using nmopt::application::chapter6::ReducedMethod;
    if (value == "steepest-descent")
      return ReducedMethod::steepest_descent;
    if (value == "l-bfgs")
      return ReducedMethod::limited_memory_bfgs;
    if (value == "bfgs")
      return ReducedMethod::bfgs;
    throw std::invalid_argument("unknown solver method '" + value + "'");
  }

  inline nmopt::application::chapter6::MeshGeneration
  parse_mesh_generation(const std::string &value)
  {
    using nmopt::application::chapter6::MeshGeneration;
    if (value == "framework-native")
      return MeshGeneration::framework_native;
    if (value == "structured-simplex")
      return MeshGeneration::structured_simplex;
    if (value == "centroid-split-simplex")
      return MeshGeneration::centroid_split_simplex;
    throw std::invalid_argument("unknown mesh generator '" + value + "'");
  }

  inline double
  parse_number_text(const std::string &text, const std::string &key)
  {
    if (text.empty() || text == "none")
      throw std::invalid_argument("parameter '" + key + "' needs a number");
    std::size_t consumed = 0;
    double      value = 0.0;
    try
      {
        value = std::stod(text, &consumed);
      }
    catch (const std::exception &)
      {
        throw std::invalid_argument("parameter '" + key + "' needs a number");
      }
    if (consumed != text.size() || !std::isfinite(value))
      throw std::invalid_argument("parameter '" + key + "' needs a finite number");
    return value;
  }

  inline unsigned int
  parse_unsigned_text(const std::string &text, const std::string &key)
  {
    const auto value = parse_number_text(text, key);
    if (value < 0.0 || value > std::numeric_limits<unsigned int>::max() ||
        value != std::floor(value))
      throw std::invalid_argument("parameter '" + key +
                                  "' needs a nonnegative integer");
    return static_cast<unsigned int>(value);
  }

  inline nmopt::solvers::ReducedStoppingCriterion
  parse_stopping_criterion(const std::string &value)
  {
    using Criterion = nmopt::solvers::ReducedStoppingCriterion;
    if (value == "automatic")
      return Criterion::automatic;
    if (value == "gradient-norm")
      return Criterion::gradient_norm;
    if (value == "relative-gradient-norm")
      return Criterion::relative_gradient_norm;
    if (value == "objective-change")
      return Criterion::objective_change;
    if (value == "step-norm")
      return Criterion::step_norm;
    throw std::invalid_argument("unknown reduced stopping criterion '" +
                                value + "'");
  }

  inline void
  require_parameter(const ParameterFile    &file,
                    const std::string_view key,
                    const std::string_view expected)
  {
    if (file.value(key) != expected)
      throw std::invalid_argument("parameter '" + std::string(key) +
                                  "' is '" + file.value(key) +
                                  "', but this adapter requires '" +
                                  std::string(expected) + "'");
  }
} // namespace nmopt::application::runner::binding
