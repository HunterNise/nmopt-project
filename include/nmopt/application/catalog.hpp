#pragma once

#include "nmopt/application/metadata.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::application
{
  enum class CatalogEntryKind
  {
    recipe,
    scenario
  };

  struct CatalogEntry
  {
    CatalogEntryKind         kind = CatalogEntryKind::recipe;
    std::string              id;
    std::string              label;
    std::string              description;
    std::string              chapter;
    std::string              recipe_id;
    std::vector<std::string> requirements;
  };

  inline CatalogEntry
  make_catalog_entry(const RecipeMetadata &metadata)
  {
    detail::validate_recipe_metadata(metadata);
    return {CatalogEntryKind::recipe,
            metadata.id,
            metadata.label,
            metadata.description,
            metadata.chapter,
            metadata.id,
            metadata.requirements};
  }

  inline CatalogEntry
  make_catalog_entry(const ScenarioMetadata &metadata)
  {
    detail::validate_scenario_metadata(metadata);
    return {CatalogEntryKind::scenario,
            metadata.id,
            metadata.label,
            metadata.description,
            metadata.chapter,
            metadata.recipe_id,
            metadata.requirements};
  }

  // The catalog is intentionally metadata-only. It is suitable for discovery
  // without erasing or owning backend-specific recipe builders.
  class ApplicationCatalog final
  {
  public:
    void
    add(CatalogEntry entry)
    {
      if (entry.id.empty())
        throw std::invalid_argument("catalog entries need an id");
      if (entry.label.empty())
        throw std::invalid_argument("catalog entries need a label");
      if (entry.chapter.empty())
        throw std::invalid_argument("catalog entries need a chapter");
      if (entry.kind == CatalogEntryKind::scenario && entry.recipe_id.empty())
        throw std::invalid_argument(
          "scenario catalog entries need a recipe id");
      if (find(entry.id) != nullptr)
        throw std::invalid_argument("duplicate application catalog id '" +
                                    entry.id + "'");
      entries_.push_back(std::move(entry));
    }

    void
    add(const RecipeMetadata &metadata)
    {
      add(make_catalog_entry(metadata));
    }

    void
    add(const ScenarioMetadata &metadata)
    {
      add(make_catalog_entry(metadata));
    }

    const CatalogEntry *
    find(const std::string &id) const noexcept
    {
      const auto entry = std::find_if(
        entries_.begin(), entries_.end(), [&id](const CatalogEntry &candidate) {
          return candidate.id == id;
        });
      return entry == entries_.end() ? nullptr : &*entry;
    }

    const std::vector<CatalogEntry> &
    entries() const noexcept
    {
      return entries_;
    }

  private:
    std::vector<CatalogEntry> entries_;
  };
} // namespace nmopt::application
