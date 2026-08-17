#pragma once

#include "nmopt/application/metadata.hpp"
#include "nmopt/semantic/v1/types.hpp"

#include <functional>
#include <utility>

namespace nmopt::application
{
  // A recipe is a typed builder for a semantic problem graph. It has no
  // backend, mesh, solver, or experiment ownership.
  template <typename Parameters>
  class ProblemRecipeT final
  {
  public:
    using parameters_type = Parameters;
    using builder_type =
      std::function<semantic::v1::ProblemSpec(const Parameters &)>;

    ProblemRecipeT(RecipeMetadata metadata, builder_type builder)
      : metadata_(std::move(metadata))
      , builder_(std::move(builder))
    {
      detail::validate_recipe_metadata(metadata_);
      if (!builder_)
        throw std::invalid_argument("a problem recipe needs a builder");
    }

    const RecipeMetadata &
    metadata() const noexcept
    {
      return metadata_;
    }

    semantic::v1::ProblemSpec
    build(const Parameters &parameters) const
    {
      return builder_(parameters);
    }

    semantic::v1::ProblemSpec
    operator()(const Parameters &parameters) const
    {
      return build(parameters);
    }

  private:
    RecipeMetadata metadata_;
    builder_type    builder_;
  };
} // namespace nmopt::application
