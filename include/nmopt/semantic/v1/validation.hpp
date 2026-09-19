#pragma once

#include "nmopt/semantic/v1/types.hpp"

#include <algorithm>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace nmopt::semantic::v1
{
  class ValidationReport
  {
  public:
    bool
    valid() const
    {
      return diagnostics_.empty();
    }

    const std::vector<Diagnostic> &
    diagnostics() const
    {
      return diagnostics_;
    }

    bool
    has_category(const DiagnosticCategory category) const
    {
      return std::any_of(diagnostics_.begin(),
                         diagnostics_.end(),
                         [category](const Diagnostic &diagnostic) {
                           return diagnostic.category == category;
                         });
    }

    void
    add(DiagnosticCategory category,
        std::string        component_id,
        std::string        capability,
        std::string        remedy)
    {
      diagnostics_.push_back(
        {category,
         std::move(component_id),
         std::move(capability),
         std::move(remedy)});
    }

  private:
    std::vector<Diagnostic> diagnostics_;
  };

  // Validation here is deliberately limited to semantic structure and stated
  // policies. Backend lowerability and formulation capabilities belong to the
  // compiler and are appended to the same diagnostic report there.
  class SemanticValidator final
  {
  public:
    ValidationReport
    validate(const ProblemSpec &specification) const;

  private:
    template <typename Component>
    using Index = std::unordered_map<std::string, const Component *>;

    template <typename Component>
    static Index<Component>
    index(const std::vector<Component> &components,
          ValidationReport &           report,
          const char *                 component_name);

    template <typename Component>
    static bool
    contains(const Index<Component> &index, const std::string &id);

    template <typename Component>
    static void
    validate_labels(const std::vector<Component> &components,
                    ValidationReport &           report,
                    const char *                 component_name);

    static void
    require_specified(const bool         specified,
                      const std::string &component_id,
                      const char *       fallback_id,
                      const char *       capability,
                      ValidationReport & report);

    static void
    validate_required_enum_fields(const ProblemSpec &specification,
                                  ValidationReport & report);

    static bool
    pairing_matches_space(const PairingSpec &pairing,
                          const std::string &space_id);

    static void
    validate_regions(const ProblemSpec &specification, ValidationReport &report);

    static void
    validate_spaces(const ProblemSpec &      specification,
                    const Index<RegionSpec> &regions,
                    ValidationReport &       report);

    static void
    validate_pairings(const ProblemSpec &     specification,
                      const Index<SpaceSpec> &spaces,
                      ValidationReport &      report);

    static void
    validate_variables(const ProblemSpec &              specification,
                       const Index<SpaceSpec> &         spaces,
                       const Index<TransformationSpec> &transformations,
                       ValidationReport &               report);

    static void
    validate_data(const ProblemSpec &     specification,
                  const Index<SpaceSpec> &spaces,
                  ValidationReport &      report);

    static void
    validate_general_scalar_data_spaces(
      const ProblemSpec &          specification,
      const Index<SpaceSpec> &     spaces,
      const Index<RegionSpec> &    regions,
      const Index<DataSpec> &      data,
      ValidationReport &           report);

    static void
    validate_transformations(const ProblemSpec &              specification,
                             const Index<VariableSpec> &     variables,
                             const Index<DataSpec> &          data,
                             const Index<SpaceSpec> &         spaces,
                             ValidationReport &               report);

    static void
    validate_equations(const ProblemSpec &             specification,
                       const Index<SpaceSpec> &        spaces,
                       const Index<PairingSpec> &      pairings,
                       const Index<ResidualTermSpec> & terms,
                       ValidationReport &              report);

    static void
    validate_terms(const ProblemSpec &             specification,
                   const Index<VariableSpec> &     variables,
                   const Index<DataSpec> &         data,
                   const Index<EquationBlockSpec> &equations,
                   const Index<SpaceSpec> &        spaces,
                   const Index<RegionSpec> &       regions,
                   ValidationReport &               report);

    static void
    validate_natural_boundary_sources(
      const ProblemSpec &      specification,
      const Index<DataSpec> &  data,
      const Index<SpaceSpec> & spaces,
      const Index<RegionSpec> &regions,
      ValidationReport &       report);

    static void
    validate_observations(const ProblemSpec &          specification,
                          const Index<VariableSpec> &  variables,
                          const Index<DataSpec> &      data,
                          const Index<RegionSpec> &    regions,
                          const Index<SpaceSpec> &     spaces,
                          const Index<PairingSpec> &   pairings,
                          ValidationReport &            report);

    static void
    validate_losses(const ProblemSpec &             specification,
                    const Index<ObservationSpec> &  observations,
                    const Index<SpaceSpec> &        spaces,
                    const Index<PairingSpec> &      pairings,
                    const Index<DataSpec> &         data,
                    ValidationReport &               report);

    static void
    validate_metrics(const ProblemSpec &            specification,
                     const Index<VariableSpec> &    variables,
                     const Index<PairingSpec> &     pairings,
                     ValidationReport &              report);

    static void
    validate_constraints(const ProblemSpec &          specification,
                         const Index<VariableSpec> &  variables,
                         const Index<DataSpec> &      data,
                         const Index<SpaceSpec> &     spaces,
                         const Index<RegionSpec> &    regions,
                         ValidationReport &            report);

    static void
    validate_policy_regions(const ProblemSpec &      specification,
                            const Index<RegionSpec> &regions,
                            ValidationReport &       report);

    static bool
    has_role(const std::vector<std::string> &ids,
             const Index<DataSpec> &         data,
             const DataRole                  role);

    static bool
    has_variable_role(const std::vector<std::string> &ids,
                      const Index<VariableSpec> &     variables,
                      const VariableRole              role);

    static void
    validate_term_signature(const ResidualTermSpec &    term,
                            const Index<VariableSpec> & variables,
                            const Index<DataSpec> &     data,
                            ValidationReport &          report);

    static void
    validate_loss_signature(const LossSpec &        loss,
                            const Index<DataSpec> & data,
                            ValidationReport &      report);

    static void
    validate_supplied_otd_block(const SuppliedOTDBlockSpec & block,
                                const SuppliedOTDBlockRole    expected_role,
                                const Index<SpaceSpec> &      spaces,
                                const Index<PairingSpec> &    pairings,
                                ValidationReport &            report);

    static void
    validate_supplied_otd_declaration(
      const ProblemSpec &             specification,
      const Index<VariableSpec> &     variables,
      const Index<SpaceSpec> &        spaces,
      const Index<PairingSpec> &      pairings,
      const Index<EquationBlockSpec> &equations,
      ValidationReport &               report);

    static void
    validate_formulation(const FormulationSpec &         formulation,
                         const Index<VariableSpec> &     variables,
                         const Index<EquationBlockSpec> &equations,
                         const Index<MetricSpec> &        metrics,
                         const Index<ConstraintSpec> &    constraints,
                         ValidationReport &               report);

    static void
    validate_policies(const ProblemSpec &specification,
                      ValidationReport & report);
  };
} // namespace nmopt::semantic::v1

#include "nmopt/semantic/v1/detail/validation_detail.hpp"
