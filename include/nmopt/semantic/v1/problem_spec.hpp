#pragma once

#include "nmopt/contract/linalg.hpp"

#include <algorithm>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// The v1 semantic graph deliberately contains no discretisation-backend
// objects. Values such as deal.II Functions and FE choices are bound by a
// compiler after this graph has been validated.
namespace nmopt::semantic::v1
{
  enum class RegionKind
  {
    volume,
    boundary
  };

  enum class SpaceTopology
  {
    h1,
    l2
  };

  enum class SpaceRole
  {
    state,
    test,
    control,
    observation,
    data
  };

  enum class VariableRole
  {
    state,
    control
  };

  enum class DataKind
  {
    function,
    scalar_constant,
    cellwise_bound
  };

  enum class DataRole
  {
    forcing,
    desired_state,
    diffusion,
    reaction,
    regularisation_weight,
    lower_bound,
    upper_bound
  };

  enum class ResidualTermKind
  {
    diffusion_reaction,
    volume_source,
    volume_control
  };

  enum class ObservationKind
  {
    volume_restriction
  };

  enum class LossKind
  {
    quadratic_tracking,
    quadratic_control_regularisation
  };

  enum class MetricKind
  {
    l2
  };

  enum class ConstraintKind
  {
    cellwise_box
  };

  enum class RequirementKind
  {
    fixed_dirichlet,
    discrete_cellwise_bounds
  };

  enum class RequirementStatus
  {
    provided,
    user_assumed,
    selected_discrete_realisation
  };

  enum class RequirementScope
  {
    continuous_semantics,
    discrete_compilation,
    both
  };

  enum class FormulationKind
  {
    reduced_dto,
    all_at_once
  };

  struct RegionSpec
  {
    std::string               id;
    std::string               label;
    RegionKind                kind;
    bool                      is_full_domain = false;
    std::vector<unsigned int> boundary_ids;
  };

  struct SpaceSpec
  {
    std::string   id;
    std::string   label;
    std::string   region_id;
    SpaceTopology topology;
    SpaceRole     role;
    bool          is_scalar = true;
  };

  struct PairingSpec
  {
    std::string id;
    std::string label;
    std::string primal_space_id;
    std::string covector_space_id;
  };

  struct VariableSpec
  {
    std::string  id;
    std::string  label;
    VariableRole role;
    std::string  space_id;
  };

  struct DataSpec
  {
    std::string id;
    std::string label;
    DataKind    kind;
    DataRole    role;
    std::string space_id;
  };

  struct ResidualTermSpec
  {
    std::string              id;
    std::string              label;
    ResidualTermKind         kind;
    std::string              equation_id;
    std::vector<std::string> variable_ids;
    std::vector<std::string> data_ids;
  };

  struct EquationBlockSpec
  {
    std::string              id;
    std::string              label;
    std::string              test_space_id;
    std::string              test_pairing_id;
    std::vector<std::string> residual_term_ids;
  };

  struct ObservationSpec
  {
    std::string     id;
    std::string     label;
    ObservationKind kind;
    std::string     input_variable_id;
    std::string     region_id;
    std::string     output_space_id;
    std::string     output_pairing_id;
  };

  struct LossSpec
  {
    std::string  id;
    std::string  label;
    LossKind     kind;
    std::string  source_observation_id;
    std::string  data_id;
    std::string  pairing_id;
  };

  struct MetricSpec
  {
    std::string id;
    std::string label;
    MetricKind  kind;
    std::string variable_id;
    std::string pairing_id;
  };

  struct ConstraintSpec
  {
    std::string    id;
    std::string    label;
    ConstraintKind kind;
    std::string    variable_id;
    std::string    lower_bound_data_id;
    std::string    upper_bound_data_id;
  };

  struct RequirementPolicySpec
  {
    std::string       id;
    std::string       subject_id;
    RequirementKind   kind;
    RequirementStatus status;
    RequirementScope  scope;
    std::string       selected_policy;
    std::string       region_id;
  };

  struct ReducedFormulationSpec
  {
    std::string     id;
    FormulationKind kind;
    std::string     state_variable_id;
    std::string     control_variable_id;
    std::string     equation_id;
    std::string     metric_id;
    std::string     constraint_id;
  };

  // This is a composition root, not a PDE-model class. The v1 compiler only
  // accepts the explicitly represented narrow stationary volume-control
  // graph; extensions add semantic node kinds and registered lowerers.
  struct ProblemSpec
  {
    std::string                        id;
    std::string                        label;
    std::vector<RegionSpec>            regions;
    std::vector<SpaceSpec>             spaces;
    std::vector<PairingSpec>           pairings;
    std::vector<VariableSpec>          variables;
    std::vector<DataSpec>              data;
    std::vector<ResidualTermSpec>      residual_terms;
    std::vector<EquationBlockSpec>     equations;
    std::vector<ObservationSpec>       observations;
    std::vector<LossSpec>              losses;
    std::vector<MetricSpec>            metrics;
    std::vector<ConstraintSpec>        constraints;
    std::vector<RequirementPolicySpec> requirement_policies;
    ReducedFormulationSpec              formulation;
  };

  enum class DiagnosticCategory
  {
    structural,
    analytical_policy,
    lowerability,
    formulation_capability
  };

  struct Diagnostic
  {
    DiagnosticCategory category;
    std::string        component_id;
    std::string        capability;
    std::string        remedy;
  };

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
    validate(const ProblemSpec &specification) const
    {
      ValidationReport report;
      if (specification.id.empty())
        report.add(DiagnosticCategory::structural,
                   "problem",
                   "stable_problem_identity",
                   "Give ProblemSpec a non-empty stable identifier.");

      const auto regions = index(specification.regions, report, "region");
      const auto spaces = index(specification.spaces, report, "space");
      const auto pairings = index(specification.pairings, report, "pairing");
      const auto variables = index(specification.variables, report, "variable");
      const auto data = index(specification.data, report, "data");
      const auto terms = index(specification.residual_terms, report,
                               "residual term");
      const auto equations = index(specification.equations, report,
                                   "equation");
      const auto observations = index(specification.observations, report,
                                      "observation");
      const auto metrics = index(specification.metrics, report, "metric");
      const auto constraints = index(specification.constraints, report,
                                     "constraint");
      index(specification.requirement_policies, report, "requirement policy");

      for (const auto &region : specification.regions)
        {
          if (region.kind == RegionKind::volume && !region.boundary_ids.empty())
            report.add(DiagnosticCategory::structural,
                       region.id,
                       "volume_region_has_no_boundary_ids",
                       "Declare boundary ids only on a boundary region.");
          if (region.kind == RegionKind::boundary && region.is_full_domain)
            report.add(DiagnosticCategory::structural,
                       region.id,
                       "boundary_region_is_not_full_domain",
                       "Use a volume region for a full-domain declaration.");
          if (region.kind == RegionKind::boundary && region.boundary_ids.empty())
            report.add(DiagnosticCategory::structural,
                       region.id,
                       "boundary_region_ids",
                       "Declare at least one boundary id for a boundary region.");
        }

      for (const auto &space : specification.spaces)
        {
          if (!contains(regions, space.region_id))
            report.add(DiagnosticCategory::structural,
                       space.id,
                       "space_region_port",
                       "Reference a declared base region.");
          if (!space.is_scalar)
            report.add(DiagnosticCategory::structural,
                       space.id,
                       "scalar_field_shape",
                       "The v1 semantic slice currently declares scalar spaces only.");
        }

      for (const auto &pairing : specification.pairings)
        {
          if (!contains(spaces, pairing.primal_space_id) ||
              !contains(spaces, pairing.covector_space_id))
            report.add(DiagnosticCategory::structural,
                       pairing.id,
                       "pairing_space_ports",
                       "Reference declared primal and covector spaces.");
        }

      for (const auto &variable : specification.variables)
        {
          const auto space = spaces.find(variable.space_id);
          if (space == spaces.end())
            report.add(DiagnosticCategory::structural,
                       variable.id,
                       "variable_space_port",
                       "Reference a declared semantic space.");
          else if ((variable.role == VariableRole::state &&
                    space->second->role != SpaceRole::state) ||
                   (variable.role == VariableRole::control &&
                    space->second->role != SpaceRole::control))
            report.add(DiagnosticCategory::structural,
                       variable.id,
                       "variable_space_role",
                       "Connect each variable to a matching state or control space.");
        }

      for (const auto &datum : specification.data)
        if (!datum.space_id.empty() && !contains(spaces, datum.space_id))
          report.add(DiagnosticCategory::structural,
                     datum.id,
                     "data_space_port",
                     "Reference a declared data space or leave scalar constants unspaced.");

      for (const auto &equation : specification.equations)
        {
          if (!contains(spaces, equation.test_space_id))
            report.add(DiagnosticCategory::structural,
                       equation.id,
                       "equation_test_space",
                       "Reference a declared scalar test space.");
          else if (spaces.at(equation.test_space_id)->role != SpaceRole::test)
            report.add(DiagnosticCategory::structural,
                       equation.id,
                       "equation_test_space_role",
                       "Reference a space declared with the test role.");
          const auto pairing = pairings.find(equation.test_pairing_id);
          if (pairing == pairings.end() ||
              pairing->second->primal_space_id != equation.test_space_id)
            report.add(DiagnosticCategory::structural,
                       equation.id,
                       "equation_test_pairing",
                       "Reference a pairing whose primal port is the equation test space.");
          for (const auto &term_id : equation.residual_term_ids)
            {
              const auto term = terms.find(term_id);
              if (term == terms.end())
                report.add(DiagnosticCategory::structural,
                           equation.id,
                           "equation_residual_term_port",
                           "Reference a declared residual term.");
              else if (term->second->equation_id != equation.id)
                report.add(DiagnosticCategory::structural,
                           term_id,
                           "residual_term_target_equation",
                           "Connect the term to the equation that owns it.");
            }
        }

      for (const auto &term : specification.residual_terms)
        {
          if (!contains(equations, term.equation_id))
            report.add(DiagnosticCategory::structural,
                       term.id,
                       "residual_term_target_equation",
                       "Reference a declared target equation block.");
          for (const auto &variable_id : term.variable_ids)
            if (!contains(variables, variable_id))
              report.add(DiagnosticCategory::structural,
                         term.id,
                         "residual_term_variable_port",
                         "Reference a declared variable input.");
          for (const auto &data_id : term.data_ids)
            if (!contains(data, data_id))
              report.add(DiagnosticCategory::structural,
                         term.id,
                         "residual_term_data_port",
                         "Reference declared immutable data.");
          validate_term_signature(term, variables, data, report);
        }

      for (const auto &observation : specification.observations)
        {
          if (!contains(variables, observation.input_variable_id))
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "observation_input_port",
                       "Reference a declared variable input.");
          if (!contains(regions, observation.region_id))
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "observation_region_port",
                       "Reference a declared volume region.");
          if (!contains(spaces, observation.output_space_id))
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "observation_output_space",
                       "Reference a declared observation space.");
          else if (spaces.at(observation.output_space_id)->role !=
                   SpaceRole::observation)
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "observation_output_space_role",
                       "Reference a space declared with the observation role.");
          const auto pairing = pairings.find(observation.output_pairing_id);
          if (pairing == pairings.end() ||
              pairing->second->primal_space_id != observation.output_space_id)
            report.add(DiagnosticCategory::structural,
                       observation.id,
                       "observation_output_pairing",
                       "Reference a pairing whose primal port is the observation output.");
        }

      for (const auto &loss : specification.losses)
        {
          if (!contains(observations, loss.source_observation_id))
            report.add(DiagnosticCategory::structural,
                       loss.id,
                       "loss_observation_port",
                       "Connect the loss to a declared observation.");
          if (!contains(data, loss.data_id))
            report.add(DiagnosticCategory::structural,
                       loss.id,
                       "loss_data_port",
                       "Reference declared target or regularisation data.");
          const auto observation = observations.find(loss.source_observation_id);
          const auto pairing = pairings.find(loss.pairing_id);
          if (observation == observations.end() || pairing == pairings.end() ||
              pairing->second->primal_space_id !=
                observation->second->output_space_id)
            report.add(DiagnosticCategory::structural,
                       loss.id,
                       "loss_pairing",
                       "Reference the pairing for the loss observation output.");
          validate_loss_signature(loss, data, report);
        }

      for (const auto &metric : specification.metrics)
        {
          const auto variable = variables.find(metric.variable_id);
          if (variable == variables.end())
            report.add(DiagnosticCategory::structural,
                       metric.id,
                       "metric_variable_port",
                       "Reference the primal variable identified by this metric.");
          const auto pairing = pairings.find(metric.pairing_id);
          if (variable == variables.end() || pairing == pairings.end() ||
              pairing->second->primal_space_id != variable->second->space_id)
            report.add(DiagnosticCategory::structural,
                       metric.id,
                       "metric_pairing",
                       "Reference the declared primal-dual pairing for the metric variable.");
          else if (variable->second->role != VariableRole::control)
            report.add(DiagnosticCategory::structural,
                       metric.id,
                       "metric_control_variable",
                       "The v1 L2 metric currently identifies the control derivative only.");
        }

      for (const auto &constraint : specification.constraints)
        {
          const auto variable = variables.find(constraint.variable_id);
          if (variable == variables.end())
            report.add(DiagnosticCategory::structural,
                       constraint.id,
                       "constraint_variable_port",
                       "Reference the constrained control variable.");
          else if (variable->second->role != VariableRole::control)
            report.add(DiagnosticCategory::structural,
                       constraint.id,
                       "constraint_control_variable",
                       "The v1 cellwise box can constrain the control variable only.");
          const auto lower = data.find(constraint.lower_bound_data_id);
          const auto upper = data.find(constraint.upper_bound_data_id);
          if (lower == data.end() || upper == data.end())
            report.add(DiagnosticCategory::structural,
                       constraint.id,
                       "constraint_bound_data_ports",
                       "Reference declared lower and upper bound data.");
          else if (lower->second->role != DataRole::lower_bound ||
                   upper->second->role != DataRole::upper_bound ||
                   lower->second->kind != DataKind::cellwise_bound ||
                   upper->second->kind != DataKind::cellwise_bound)
            report.add(DiagnosticCategory::structural,
                       constraint.id,
                       "constraint_cellwise_bound_data",
                       "Use declared cellwise lower and upper bound data.");
        }

      for (const auto &policy : specification.requirement_policies)
        if (!policy.region_id.empty() && !contains(regions, policy.region_id))
          report.add(DiagnosticCategory::structural,
                     policy.id,
                     "requirement_policy_region",
                     "Reference a declared region when a policy is region-specific.");

      validate_formulation(specification.formulation,
                           variables,
                           equations,
                           metrics,
                           constraints,
                           report);
      validate_policies(specification, report);
      return report;
    }

  private:
    template <typename Component>
    using Index = std::unordered_map<std::string, const Component *>;

    template <typename Component>
    static Index<Component>
    index(const std::vector<Component> &components,
          ValidationReport &           report,
          const char *                  component_name)
    {
      Index<Component> result;
      for (const auto &component : components)
        {
          if (component.id.empty())
            {
              report.add(DiagnosticCategory::structural,
                         component_name,
                         "stable_component_identity",
                         "Give every semantic component a non-empty identifier.");
              continue;
            }
          if (!result.emplace(component.id, &component).second)
            report.add(DiagnosticCategory::structural,
                       component.id,
                       "unique_component_identity",
                       "Use a unique identifier for this component kind.");
        }
      return result;
    }

    template <typename Component>
    static bool
    contains(const Index<Component> &index, const std::string &id)
    {
      return index.find(id) != index.end();
    }

    static bool
    has_role(const std::vector<std::string> &ids,
             const Index<DataSpec> &         data,
             const DataRole                  role)
    {
      return std::any_of(ids.begin(), ids.end(), [&data, role](const auto &id) {
        const auto entry = data.find(id);
        return entry != data.end() && entry->second->role == role;
      });
    }

    static bool
    has_variable_role(const std::vector<std::string> &ids,
                      const Index<VariableSpec> &     variables,
                      const VariableRole               role)
    {
      return std::any_of(ids.begin(), ids.end(), [&variables, role](const auto &id) {
        const auto entry = variables.find(id);
        return entry != variables.end() && entry->second->role == role;
      });
    }

    static void
    validate_term_signature(const ResidualTermSpec &    term,
                            const Index<VariableSpec> & variables,
                            const Index<DataSpec> &     data,
                            ValidationReport &          report)
    {
      bool valid = false;
      switch (term.kind)
        {
          case ResidualTermKind::diffusion_reaction:
            valid = term.variable_ids.size() == 1 && term.data_ids.size() == 2 &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::state) &&
                    has_role(term.data_ids, data, DataRole::diffusion) &&
                    has_role(term.data_ids, data, DataRole::reaction);
            break;
          case ResidualTermKind::volume_source:
            valid = term.variable_ids.empty() && term.data_ids.size() == 1 &&
                    has_role(term.data_ids, data, DataRole::forcing);
            break;
          case ResidualTermKind::volume_control:
            valid = term.variable_ids.size() == 1 &&
                    has_variable_role(term.variable_ids, variables,
                                      VariableRole::control) &&
                    term.data_ids.empty();
            break;
        }
      if (!valid)
        report.add(DiagnosticCategory::structural,
                   term.id,
                   "residual_term_signature",
                   "Supply the declared v1 term inputs and no undeclared ports.");
    }

    static void
    validate_loss_signature(const LossSpec &        loss,
                            const Index<DataSpec> & data,
                            ValidationReport &      report)
    {
      const auto datum = data.find(loss.data_id);
      if (datum == data.end())
        return;

      const bool valid =
        (loss.kind == LossKind::quadratic_tracking &&
         datum->second->role == DataRole::desired_state) ||
        (loss.kind == LossKind::quadratic_control_regularisation &&
         datum->second->role == DataRole::regularisation_weight);
      if (!valid)
        report.add(DiagnosticCategory::structural,
                   loss.id,
                   "loss_data_role",
                   "Use desired-state data for tracking and a scalar regularisation weight for control loss.");
    }

    static void
    validate_formulation(const ReducedFormulationSpec & formulation,
                         const Index<VariableSpec> &    variables,
                         const Index<EquationBlockSpec> & equations,
                         const Index<MetricSpec> &       metrics,
                         const Index<ConstraintSpec> &   constraints,
                         ValidationReport &              report)
    {
      const auto state = variables.find(formulation.state_variable_id);
      const auto control = variables.find(formulation.control_variable_id);
      if (state == variables.end() || state->second->role != VariableRole::state)
        report.add(DiagnosticCategory::structural,
                   formulation.id,
                   "formulation_state_variable",
                   "Select one declared state variable for elimination.");
      if (control == variables.end() ||
          control->second->role != VariableRole::control)
        report.add(DiagnosticCategory::structural,
                   formulation.id,
                   "formulation_control_variable",
                   "Select one declared control variable for optimisation.");
      if (!contains(equations, formulation.equation_id))
        report.add(DiagnosticCategory::structural,
                   formulation.id,
                   "formulation_equation",
                   "Reference the equation that defines the state.");
      if (!contains(metrics, formulation.metric_id))
        report.add(DiagnosticCategory::structural,
                   formulation.id,
                   "formulation_metric",
                   "Reference a declared search metric.");
      if (!formulation.constraint_id.empty() &&
          !contains(constraints, formulation.constraint_id))
        report.add(DiagnosticCategory::structural,
                   formulation.id,
                   "formulation_constraint",
                   "Reference a declared constraint or leave this port empty.");
    }

    static void
    validate_policies(const ProblemSpec &specification,
                      ValidationReport & report)
    {
      const auto has_policy = [&specification](const std::string &subject,
                                                const RequirementKind kind) {
        return std::any_of(
          specification.requirement_policies.begin(),
          specification.requirement_policies.end(),
          [&subject, kind](const RequirementPolicySpec &policy) {
            return policy.subject_id == subject && policy.kind == kind &&
                   policy.status ==
                     RequirementStatus::selected_discrete_realisation &&
                   !policy.selected_policy.empty();
          });
      };

      for (const auto &variable : specification.variables)
        if (variable.role == VariableRole::state &&
            !has_policy(variable.id, RequirementKind::fixed_dirichlet))
          report.add(DiagnosticCategory::analytical_policy,
                     variable.id,
                     "fixed_dirichlet_realisation",
                     "Declare the selected fixed-Dirichlet discrete policy.");

      for (const auto &constraint : specification.constraints)
        if (constraint.kind == ConstraintKind::cellwise_box &&
            !has_policy(constraint.id,
                        RequirementKind::discrete_cellwise_bounds))
          report.add(DiagnosticCategory::analytical_policy,
                     constraint.id,
                     "cellwise_bound_realisation",
                     "Declare the FE_DGQ(0) coefficientwise bound policy.");
    }
  };

  inline ProblemSpec
  make_scalar_diffusion_reaction_problem(const bool with_cellwise_box = false)
  {
    ProblemSpec specification;
    specification.id = "scalar_diffusion_reaction_volume_control";
    specification.label = "Scalar diffusion-reaction volume control";
    specification.regions = {
      {"domain", "Full volume domain", RegionKind::volume, true, {}},
      {"dirichlet_boundary", "Homogeneous Dirichlet boundary", RegionKind::boundary,
       false, {0}}};
    specification.spaces = {
      {"state_space", "State", "domain", SpaceTopology::h1, SpaceRole::state},
      {"state_test_space", "State test", "domain", SpaceTopology::h1,
       SpaceRole::test},
      {"control_space", "Cellwise control", "domain", SpaceTopology::l2,
       SpaceRole::control},
      {"state_observation_space", "State observation", "domain",
       SpaceTopology::l2, SpaceRole::observation},
      {"control_observation_space", "Control observation", "domain",
       SpaceTopology::l2, SpaceRole::observation}};
    specification.pairings = {
      {"state_pairing", "State coefficient pairing", "state_space", "state_space"},
      {"state_test_pairing", "State-test coefficient pairing",
       "state_test_space", "state_test_space"},
      {"control_pairing", "Control coefficient pairing", "control_space",
       "control_space"},
      {"state_observation_pairing", "State-observation coefficient pairing",
       "state_observation_space", "state_observation_space"},
      {"control_observation_pairing",
       "Control-observation coefficient pairing", "control_observation_space",
       "control_observation_space"}};
    specification.variables = {
      {"state", "State", VariableRole::state, "state_space"},
      {"control", "Control", VariableRole::control, "control_space"}};
    specification.data = {
      {"forcing", "Volume forcing", DataKind::function, DataRole::forcing,
       "state_test_space"},
      {"desired_state", "Desired state", DataKind::function,
       DataRole::desired_state, "state_observation_space"},
      {"diffusion", "Diffusion coefficient", DataKind::scalar_constant,
       DataRole::diffusion, ""},
      {"reaction", "Reaction coefficient", DataKind::scalar_constant,
       DataRole::reaction, ""},
      {"regularisation_weight", "Control regularisation", DataKind::scalar_constant,
       DataRole::regularisation_weight, ""}};
    specification.residual_terms = {
      {"diffusion_reaction", "Diffusion and reaction",
       ResidualTermKind::diffusion_reaction, "state_equation", {"state"},
       {"diffusion", "reaction"}},
      {"volume_source", "Volume source", ResidualTermKind::volume_source,
       "state_equation", {}, {"forcing"}},
      {"volume_control", "Volume control", ResidualTermKind::volume_control,
       "state_equation", {"control"}, {}}};
    specification.equations = {
      {"state_equation", "State residual", "state_test_space",
       "state_test_pairing",
       {"diffusion_reaction", "volume_source", "volume_control"}}};
    specification.observations = {
      {"state_observation", "Full-domain state restriction",
       ObservationKind::volume_restriction, "state", "domain",
       "state_observation_space", "state_observation_pairing"},
      {"control_observation", "Full-domain control restriction",
       ObservationKind::volume_restriction, "control", "domain",
       "control_observation_space", "control_observation_pairing"}};
    specification.losses = {
      {"state_tracking", "Quadratic tracking", LossKind::quadratic_tracking,
       "state_observation", "desired_state", "state_observation_pairing"},
      {"control_regularisation", "Quadratic control regularisation",
       LossKind::quadratic_control_regularisation, "control_observation",
       "regularisation_weight", "control_observation_pairing"}};
    specification.metrics = {
      {"control_l2_metric", "Cellwise L2 metric", MetricKind::l2, "control",
       "control_pairing"}};
    specification.requirement_policies = {
      {"state_fixed_dirichlet", "state", RequirementKind::fixed_dirichlet,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "homogeneous full-vector Dirichlet rows", "dirichlet_boundary"}};
    specification.formulation = {"reduced_dto", FormulationKind::reduced_dto,
                                  "state", "control", "state_equation",
                                  "control_l2_metric", ""};

    if (with_cellwise_box)
      {
        specification.data.push_back(
          {"lower_bound", "Control lower bound", DataKind::cellwise_bound,
           DataRole::lower_bound, "control_space"});
        specification.data.push_back(
          {"upper_bound", "Control upper bound", DataKind::cellwise_bound,
           DataRole::upper_bound, "control_space"});
        specification.constraints.push_back(
          {"control_box", "Cellwise L2 box", ConstraintKind::cellwise_box,
           "control", "lower_bound", "upper_bound"});
        specification.requirement_policies.push_back(
          {"control_box_policy", "control_box",
           RequirementKind::discrete_cellwise_bounds,
           RequirementStatus::selected_discrete_realisation,
           RequirementScope::discrete_compilation,
           "FE_DGQ(0) coefficientwise clipping in l2_cellwise", "domain"});
        specification.formulation.constraint_id = "control_box";
      }

    return specification;
  }
} // namespace nmopt::semantic::v1
