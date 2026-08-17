#pragma once

#include <stdexcept>
#include <string>
#include <vector>

namespace nmopt::application
{
  struct RecipeMetadata
  {
    std::string              id;
    std::string              label;
    std::string              description;
    std::string              chapter;
    std::vector<std::string> requirements;
  };

  struct ScenarioMetadata
  {
    std::string              id;
    std::string              label;
    std::string              description;
    std::string              chapter;
    std::string              recipe_id;
    std::vector<std::string> requirements;
  };

  namespace detail
  {
    template <typename Metadata>
    inline void
    validate_common_metadata(const Metadata &metadata)
    {
      if (metadata.id.empty())
        throw std::invalid_argument("application metadata needs an id");
      if (metadata.label.empty())
        throw std::invalid_argument("application metadata needs a label");
      if (metadata.chapter.empty())
        throw std::invalid_argument("application metadata needs a chapter");
    }

    inline void
    validate_recipe_metadata(const RecipeMetadata &metadata)
    {
      validate_common_metadata(metadata);
    }

    inline void
    validate_scenario_metadata(const ScenarioMetadata &metadata)
    {
      validate_common_metadata(metadata);
      if (metadata.recipe_id.empty())
        throw std::invalid_argument(
          "scenario metadata needs the id of its problem recipe");
    }
  } // namespace detail
} // namespace nmopt::application
