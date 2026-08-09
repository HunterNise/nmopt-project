#pragma once

#include "nmopt/semantic/v1/validation.hpp"

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace nmopt::semantic::v1
{
  // A validated, closed semantic graph with every stable ID resolved once.
  // The view borrows its ProblemSpec and is intentionally confined to one
  // validation/compilation operation.
  class ResolvedProblemView final
  {
  public:
    const ProblemSpec &
    specification() const
    {
      return *specification_;
    }

    const RegionSpec &
    region(const std::string &id) const
    {
      return *regions_.at(id);
    }

    const SpaceSpec &
    space(const std::string &id) const
    {
      return *spaces_.at(id);
    }

    const PairingSpec &
    pairing(const std::string &id) const
    {
      return *pairings_.at(id);
    }

    const VariableSpec &
    variable(const std::string &id) const
    {
      return *variables_.at(id);
    }

    const DataSpec &
    datum(const std::string &id) const
    {
      return *data_.at(id);
    }

    const TransformationSpec &
    transformation(const std::string &id) const
    {
      return *transformations_.at(id);
    }

    const ResidualTermSpec &
    residual_term(const std::string &id) const
    {
      return *residual_terms_.at(id);
    }

    const EquationBlockSpec &
    equation(const std::string &id) const
    {
      return *equations_.at(id);
    }

    const ObservationSpec &
    observation(const std::string &id) const
    {
      return *observations_.at(id);
    }

    const LossSpec &
    loss(const std::string &id) const
    {
      return *losses_.at(id);
    }

    const MetricSpec &
    metric(const std::string &id) const
    {
      return *metrics_.at(id);
    }

    const ConstraintSpec &
    constraint(const std::string &id) const
    {
      return *constraints_.at(id);
    }

    const RequirementPolicySpec &
    requirement(const std::string &id) const
    {
      return *requirements_.at(id);
    }

  private:
    template <typename Component>
    using Index = std::unordered_map<std::string, const Component *>;

    explicit ResolvedProblemView(const ProblemSpec &specification)
      : specification_(&specification)
      , regions_(index(specification.regions))
      , spaces_(index(specification.spaces))
      , pairings_(index(specification.pairings))
      , variables_(index(specification.variables))
      , data_(index(specification.data))
      , transformations_(index(specification.transformations))
      , residual_terms_(index(specification.residual_terms))
      , equations_(index(specification.equations))
      , observations_(index(specification.observations))
      , losses_(index(specification.losses))
      , metrics_(index(specification.metrics))
      , constraints_(index(specification.constraints))
      , requirements_(index(specification.requirement_policies))
    {}

    template <typename Component>
    static Index<Component>
    index(const std::vector<Component> &components)
    {
      Index<Component> result;
      result.reserve(components.size());
      for (const auto &component : components)
        result.emplace(component.id, &component);
      return result;
    }

    const ProblemSpec *                  specification_;
    Index<RegionSpec>                    regions_;
    Index<SpaceSpec>                     spaces_;
    Index<PairingSpec>                   pairings_;
    Index<VariableSpec>                  variables_;
    Index<DataSpec>                      data_;
    Index<TransformationSpec>            transformations_;
    Index<ResidualTermSpec>              residual_terms_;
    Index<EquationBlockSpec>             equations_;
    Index<ObservationSpec>               observations_;
    Index<LossSpec>                      losses_;
    Index<MetricSpec>                    metrics_;
    Index<ConstraintSpec>                constraints_;
    Index<RequirementPolicySpec>         requirements_;

    friend class SemanticResolver;
  };

  struct SemanticResolution
  {
    ValidationReport                   diagnostics;
    std::optional<ResolvedProblemView> problem;

    bool
    succeeded() const
    {
      return diagnostics.valid() && problem.has_value();
    }
  };

  class SemanticResolver final
  {
  public:
    SemanticResolution
    resolve(const ProblemSpec &specification) const
    {
      SemanticResolution result;
      result.diagnostics = SemanticValidator().validate(specification);
      if (result.diagnostics.valid())
        result.problem = ResolvedProblemView(specification);
      return result;
    }
  };
} // namespace nmopt::semantic::v1
