#pragma once

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>
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

  struct RunSetComparisonCoordinates
  {
    std::string rows;
    std::string columns;
    std::string group_by;
  };

  struct RunSetParameterProvenance
  {
    std::filesystem::path file;
    std::string            content_hash;
  };

  struct RunSetArtifactCoordinateComponent
  {
    std::string axis;
    std::string value;
  };

  struct RunSetCombination
  {
    ParameterCombination                       values;
    std::vector<RunSetArtifactCoordinateComponent> artifact_coordinates;
  };

  struct RunSetPlan
  {
    std::string                                  benchmark_id;
    std::vector<ParameterAxis>                   matrix_axes;
    std::map<std::string, std::string>           selection;
    std::vector<ParameterCombination>            excluded_combinations;
    std::vector<RunSetCombination>                resolved_combinations;
    RunSetComparisonCoordinates                   comparison;
    RunSetParameterProvenance                     parameter_provenance;
  };

  inline std::vector<RunSetArtifactCoordinateComponent>
  default_artifact_coordinate_components(
    const std::vector<ParameterAxis> &axes,
    const ParameterCombination       &combination)
  {
    std::vector<RunSetArtifactCoordinateComponent> result;
    result.reserve(axes.size());
    for (const auto &axis : axes)
      {
        const auto selected = combination.values.find(axis.id);
        if (selected == combination.values.end())
          throw std::invalid_argument(
            "run-set combination is missing matrix axis '" + axis.id + "'");
        result.push_back({axis.id, selected->second});
      }
    return result;
  }

  inline void
  validate_run_set_plan(const RunSetPlan &plan)
  {
    if (plan.benchmark_id.empty())
      throw std::invalid_argument("run-set plan needs a benchmark identifier");
    if (plan.matrix_axes.empty())
      throw std::invalid_argument("run-set plan needs at least one matrix axis");

    for (std::size_t axis_index = 0; axis_index < plan.matrix_axes.size();
         ++axis_index)
      {
        const auto &axis = plan.matrix_axes[axis_index];
        if (axis.id.empty())
          throw std::invalid_argument("run-set matrix axes need nonempty IDs");
        if (std::any_of(plan.matrix_axes.begin(),
                        plan.matrix_axes.begin() +
                          static_cast<std::ptrdiff_t>(axis_index),
                        [&](const auto &candidate) {
                          return candidate.id == axis.id;
                        }))
          throw std::invalid_argument("run-set matrix axes must have unique IDs");
        if (axis.values.empty())
          throw std::invalid_argument("run-set matrix axes need values");
        for (std::size_t value_index = 0; value_index < axis.values.size();
             ++value_index)
          {
            if (axis.values[value_index].empty())
              throw std::invalid_argument(
                "run-set matrix axis values must be nonempty");
            if (std::find(axis.values.begin(),
                          axis.values.begin() +
                            static_cast<std::ptrdiff_t>(value_index),
                          axis.values[value_index]) !=
                axis.values.begin() + static_cast<std::ptrdiff_t>(value_index))
              throw std::invalid_argument(
                "run-set matrix axes must not repeat values");
          }
      }

    const auto validate_combination = [&](const ParameterCombination &combination,
                                          const char                  *kind) {
      if (combination.values.size() != plan.matrix_axes.size())
        throw std::invalid_argument(std::string("run-set ") + kind +
                                    " must name every matrix axis exactly once");
      for (const auto &axis : plan.matrix_axes)
        {
          const auto selected = combination.values.find(axis.id);
          if (selected == combination.values.end())
            throw std::invalid_argument(std::string("run-set ") + kind +
                                        " is missing matrix axis '" + axis.id +
                                        "'");
          if (std::find(axis.values.begin(), axis.values.end(), selected->second) ==
              axis.values.end())
            throw std::invalid_argument(
              std::string("run-set ") + kind + " selects undeclared value '" +
              selected->second + "' on matrix axis '" + axis.id + "'");
        }
    };

    for (const auto &combination : plan.excluded_combinations)
      validate_combination(combination, "exclusion");

    if (plan.resolved_combinations.empty())
      throw std::invalid_argument("run-set plan resolves to no combinations");
    for (const auto &combination : plan.resolved_combinations)
      {
        validate_combination(combination.values, "combination");
        if (combination.artifact_coordinates.size() != plan.matrix_axes.size())
          throw std::invalid_argument(
            "run-set combinations need one artifact coordinate per matrix axis");
        for (std::size_t index = 0; index < plan.matrix_axes.size(); ++index)
          if (combination.artifact_coordinates[index].axis !=
                plan.matrix_axes[index].id ||
              combination.artifact_coordinates[index].value !=
                combination.values.values.at(plan.matrix_axes[index].id))
            throw std::invalid_argument(
              "run-set artifact coordinates must follow matrix axis order");
      }
  }
} // namespace nmopt::application::runner
