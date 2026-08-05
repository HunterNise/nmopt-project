#pragma once

#include "nmopt/semantic/v1/types.hpp"

#include <algorithm>

namespace nmopt::semantic::v1
{
  // This is the current reference graph used to compare direct v0 assembly
  // with semantic compilation. It is a factory, not a PDE problem class.
  inline ProblemSpec
  make_scalar_diffusion_reaction_problem(const bool with_cellwise_box = false)
  {
    ProblemSpec specification;
    specification.id = "scalar_diffusion_reaction_volume_control";
    specification.label = "Scalar diffusion-reaction volume control";
    specification.regions = {
      {"domain", "Full volume domain", RegionKind::volume, true, {}, {}},
      {"dirichlet_boundary", "Homogeneous Dirichlet boundary", RegionKind::boundary,
       false, {0}, {}}};
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
      {"state", "State", VariableRole::state, "state_space", ""},
      {"control", "Control", VariableRole::control, "control_space", ""}};
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
       {"diffusion", "reaction"}, ""},
      {"volume_source", "Volume source", ResidualTermKind::volume_source,
       "state_equation", {}, {"forcing"}, ""},
      {"volume_control", "Volume control", ResidualTermKind::volume_control,
       "state_equation", {"control"}, {}, ""}};
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
       "homogeneous full-vector Dirichlet rows", "dirichlet_boundary"},
      {"desired_state_quadrature_policy", "desired_state",
       RequirementKind::analytic_quadrature_evaluation,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "analytic Function evaluated at selected volume quadrature", "domain"}};
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

  inline ProblemSpec
  make_subdomain_tracking_scalar_diffusion_reaction_problem(
    const unsigned int observed_material_id,
    const bool         with_cellwise_box = false)
  {
    ProblemSpec specification =
      make_scalar_diffusion_reaction_problem(with_cellwise_box);
    specification.id = "scalar_diffusion_reaction_subdomain_tracking";
    specification.label = "Scalar diffusion-reaction with subdomain tracking";
    specification.regions.push_back(
      {"observation_subdomain", "Material subdomain observation region",
       RegionKind::volume, false, {}, {observed_material_id}});
    specification.spaces.at(3).region_id = "observation_subdomain";
    specification.observations.at(0).label = "Subdomain state restriction";
    specification.observations.at(0).region_id = "observation_subdomain";
    for (auto &policy : specification.requirement_policies)
      if (policy.subject_id == "desired_state" &&
          policy.kind == RequirementKind::analytic_quadrature_evaluation)
        policy.region_id = "observation_subdomain";
    return specification;
  }

  // P2.3's first half changes the objective, not the search geometry.  The
  // control has continuous FE_Q coordinates and the declared loss is the
  // H1 norm, while the selected algorithmic metric remains L2.
  inline ProblemSpec
  make_h1_regularised_scalar_diffusion_reaction_problem()
  {
    ProblemSpec specification = make_scalar_diffusion_reaction_problem();
    specification.id = "scalar_diffusion_reaction_h1_control_regularisation";
    specification.label =
      "Scalar diffusion-reaction with H1 control regularisation";
    specification.spaces.at(2) =
      {"control_space", "Continuous control", "domain", SpaceTopology::h1,
       SpaceRole::control};
    specification.spaces.at(4) =
      {"control_observation_space", "Continuous control observation", "domain",
       SpaceTopology::h1, SpaceRole::observation};
    specification.pairings.at(2) =
      {"control_pairing", "Continuous control coefficient pairing",
       "control_space", "control_space"};
    specification.pairings.at(4) =
      {"control_observation_pairing", "Continuous control observation pairing",
       "control_observation_space", "control_observation_space"};
    specification.observations.at(1) =
      {"control_observation", "Continuous control identity observation",
       ObservationKind::volume_restriction, "control", "domain",
       "control_observation_space", "control_observation_pairing"};
    specification.losses.at(1) =
      {"control_h1_regularisation", "Quadratic H1 control regularisation",
       LossKind::quadratic_h1_control_regularisation, "control_observation",
       "regularisation_weight", "control_observation_pairing"};
    specification.metrics.at(0) =
      {"control_l2_metric", "Continuous-control L2 metric", MetricKind::l2,
       "control", "control_pairing"};
    return specification;
  }

  // This changes only the Riesz map used for search directions. The inherited
  // H1 regularisation loss and every residual/objective action are unchanged.
  inline ProblemSpec
  make_h1_metric_scalar_diffusion_reaction_problem()
  {
    ProblemSpec specification =
      make_h1_regularised_scalar_diffusion_reaction_problem();
    specification.id = "scalar_diffusion_reaction_h1_control_metric";
    specification.label = "Scalar diffusion-reaction with H1 control metric";
    specification.metrics.at(0) =
      {"control_h1_metric", "Continuous-control H1 metric", MetricKind::h1,
       "control", "control_pairing"};
    specification.formulation.metric_id = "control_h1_metric";
    return specification;
  }

  // P3.1 uses the binary reduced DTO decision port for a physical diffusion
  // coefficient rather than a source control. Positivity is explicit through
  // the required cellwise box; a logarithmic parameterisation is deliberately
  // left to a later transformation realization.
  inline ProblemSpec
  make_coefficient_identification_problem()
  {
    ProblemSpec specification = make_scalar_diffusion_reaction_problem();
    specification.id = "scalar_diffusion_reaction_coefficient_identification";
    specification.label =
      "Scalar diffusion-reaction coefficient identification";
    specification.spaces.at(2) =
      {"parameter_space", "Cellwise diffusion parameter", "domain",
       SpaceTopology::l2, SpaceRole::parameter};
    specification.spaces.at(4) =
      {"parameter_observation_space", "Cellwise parameter observation",
       "domain", SpaceTopology::l2, SpaceRole::observation};
    specification.pairings.at(2) =
      {"parameter_pairing", "Parameter coefficient pairing", "parameter_space",
       "parameter_space"};
    specification.pairings.at(4) =
      {"parameter_observation_pairing", "Parameter observation pairing",
       "parameter_observation_space", "parameter_observation_space"};
    specification.variables.at(1) =
      {"diffusion_parameter", "Diffusion parameter", VariableRole::parameter,
       "parameter_space", ""};
    specification.data.erase(
      std::remove_if(specification.data.begin(),
                     specification.data.end(),
                     [](const DataSpec &datum) {
                       return datum.role == DataRole::diffusion;
                     }),
      specification.data.end());
    specification.data.push_back(
      {"parameter_lower_bound", "Strictly positive diffusion lower bound",
       DataKind::cellwise_bound, DataRole::lower_bound, "parameter_space"});
    specification.data.push_back(
      {"parameter_upper_bound", "Diffusion upper bound",
       DataKind::cellwise_bound, DataRole::upper_bound, "parameter_space"});
    specification.residual_terms.at(0) =
      {"parameter_diffusion_reaction", "Parameter diffusion and reaction",
       ResidualTermKind::parameter_diffusion_reaction, "state_equation",
       {"state", "diffusion_parameter"}, {"reaction"}, ""};
    specification.residual_terms.erase(specification.residual_terms.begin() + 2);
    specification.equations.at(0).residual_term_ids =
      {"parameter_diffusion_reaction", "volume_source"};
    specification.observations.at(1) =
      {"parameter_observation", "Full-domain parameter restriction",
       ObservationKind::volume_restriction, "diffusion_parameter", "domain",
       "parameter_observation_space", "parameter_observation_pairing"};
    specification.losses.at(1) =
      {"parameter_regularisation", "Quadratic parameter regularisation",
       LossKind::quadratic_parameter_regularisation, "parameter_observation",
       "regularisation_weight", "parameter_observation_pairing"};
    specification.metrics.at(0) =
      {"parameter_l2_metric", "Cellwise parameter L2 metric", MetricKind::l2,
       "diffusion_parameter", "parameter_pairing"};
    specification.constraints = {
      {"parameter_box", "Positive cellwise diffusion box",
       ConstraintKind::cellwise_box, "diffusion_parameter",
       "parameter_lower_bound", "parameter_upper_bound"}};
    specification.requirement_policies.push_back(
      {"parameter_positive_box_policy", "parameter_box",
       RequirementKind::discrete_cellwise_bounds,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "strictly positive FE_DGQ(0) lower bound with coefficientwise l2_cellwise_parameter clipping",
       "domain"});
    specification.formulation.control_variable_id = "diffusion_parameter";
    specification.formulation.metric_id = "parameter_l2_metric";
    specification.formulation.constraint_id = "parameter_box";
    return specification;
  }

  // The first boundary-control graph deliberately has a different control
  // space and residual term from volume control. It is a Neumann trace
  // pairing, not a Dirichlet lifting or a generic boundary-load switch.
  inline ProblemSpec
  make_neumann_boundary_control_problem(const bool with_facewise_box = false)
  {
    ProblemSpec specification;
    specification.id = "scalar_diffusion_reaction_neumann_boundary_control";
    specification.label = "Scalar diffusion-reaction Neumann boundary control";
    specification.regions = {
      {"domain", "Full volume domain", RegionKind::volume, true, {}, {}},
      {"dirichlet_boundary", "Homogeneous Dirichlet boundary", RegionKind::boundary,
       false, {0}, {}},
      {"control_boundary", "Neumann control boundary", RegionKind::boundary,
       false, {1}, {}},
      {"observation_boundary", "Boundary tracking region", RegionKind::boundary,
       false, {2}, {}}};
    specification.spaces = {
      {"state_space", "State", "domain", SpaceTopology::h1, SpaceRole::state},
      {"state_test_space", "State test", "domain", SpaceTopology::h1,
       SpaceRole::test},
      {"control_space", "Facewise Neumann control", "control_boundary",
       SpaceTopology::l2, SpaceRole::control},
      {"state_observation_space", "Boundary state trace", "observation_boundary",
       SpaceTopology::l2, SpaceRole::observation},
      {"control_observation_space", "Boundary control restriction",
       "control_boundary", SpaceTopology::l2, SpaceRole::observation}};
    specification.pairings = {
      {"state_pairing", "State coefficient pairing", "state_space", "state_space"},
      {"state_test_pairing", "State-test coefficient pairing",
       "state_test_space", "state_test_space"},
      {"control_pairing", "Facewise control coefficient pairing", "control_space",
       "control_space"},
      {"state_observation_pairing", "Boundary state-trace pairing",
       "state_observation_space", "state_observation_space"},
      {"control_observation_pairing", "Boundary control pairing",
       "control_observation_space", "control_observation_space"}};
    specification.variables = {
      {"state", "State", VariableRole::state, "state_space", ""},
      {"control", "Neumann control", VariableRole::control, "control_space", ""}};
    specification.data = {
      {"forcing", "Volume forcing", DataKind::function, DataRole::forcing,
       "state_test_space"},
      {"desired_state", "Desired boundary trace", DataKind::function,
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
       {"diffusion", "reaction"}, ""},
      {"volume_source", "Volume source", ResidualTermKind::volume_source,
       "state_equation", {}, {"forcing"}, ""},
      {"neumann_control", "Neumann control trace pairing",
       ResidualTermKind::neumann_control, "state_equation", {"control"}, {},
       "control_boundary"}};
    specification.equations = {
      {"state_equation", "State residual", "state_test_space",
       "state_test_pairing",
       {"diffusion_reaction", "volume_source", "neumann_control"}}};
    specification.observations = {
      {"state_boundary_trace", "Boundary state trace", ObservationKind::boundary_trace,
       "state", "observation_boundary", "state_observation_space",
       "state_observation_pairing"},
      {"control_boundary_restriction", "Boundary control restriction",
       ObservationKind::boundary_restriction, "control", "control_boundary",
       "control_observation_space", "control_observation_pairing"}};
    specification.losses = {
      {"state_tracking", "Quadratic boundary tracking", LossKind::quadratic_tracking,
       "state_boundary_trace", "desired_state", "state_observation_pairing"},
      {"control_regularisation", "Quadratic boundary control regularisation",
       LossKind::quadratic_control_regularisation, "control_boundary_restriction",
       "regularisation_weight", "control_observation_pairing"}};
    specification.metrics = {
      {"control_l2_metric", "Facewise L2 metric", MetricKind::l2, "control",
       "control_pairing"}};
    specification.requirement_policies = {
      {"state_fixed_dirichlet", "state", RequirementKind::fixed_dirichlet,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "homogeneous full-vector Dirichlet rows", "dirichlet_boundary"},
      {"neumann_control_trace_policy", "neumann_control",
       RequirementKind::boundary_trace,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "facewise-constant control with FEFaceValues trace pairing",
       "control_boundary"},
      {"state_boundary_trace_policy", "state_boundary_trace",
       RequirementKind::boundary_trace,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "FE_Q state trace at selected boundary face quadrature",
       "observation_boundary"},
      {"desired_state_quadrature_policy", "desired_state",
       RequirementKind::analytic_quadrature_evaluation,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "analytic Function evaluated at selected boundary face quadrature",
       "observation_boundary"}};
    specification.formulation = {"reduced_dto", FormulationKind::reduced_dto,
                                  "state", "control", "state_equation",
                                  "control_l2_metric", ""};

    if (with_facewise_box)
      {
        specification.data.push_back(
          {"lower_bound", "Boundary control lower bound", DataKind::facewise_bound,
           DataRole::lower_bound, "control_space"});
        specification.data.push_back(
          {"upper_bound", "Boundary control upper bound", DataKind::facewise_bound,
           DataRole::upper_bound, "control_space"});
        specification.constraints.push_back(
          {"control_box", "Facewise L2 box", ConstraintKind::facewise_box,
           "control", "lower_bound", "upper_bound"});
        specification.requirement_policies.push_back(
          {"control_box_policy", "control_box",
           RequirementKind::discrete_facewise_bounds,
           RequirementStatus::selected_discrete_realisation,
           RequirementScope::discrete_compilation,
           "one facewise-constant coefficient per marked boundary face; coefficientwise clipping in l2_facewise",
           "control_boundary"});
        specification.formulation.constraint_id = "control_box";
      }

    return specification;
  }

  // The pure-Neumann variant keeps the natural boundary residual of the
  // preceding graph but replaces the fixed-boundary uniqueness policy by the
  // selected discrete mean constraint.  The auxiliary multiplier belongs to
  // the compiled solve, not to the user-facing state/control graph.
  inline ProblemSpec
  make_pure_neumann_boundary_control_problem()
  {
    ProblemSpec specification = make_neumann_boundary_control_problem();
    specification.id = "scalar_diffusion_pure_neumann_boundary_control";
    specification.label =
      "Scalar diffusion pure-Neumann boundary control with mean constraint";
    specification.regions = {
      {"domain", "Full volume domain", RegionKind::volume, true, {}, {}},
      {"control_boundary", "Pure-Neumann control boundary", RegionKind::boundary,
       false, {0}, {}},
      {"observation_boundary", "Pure-Neumann boundary tracking region",
       RegionKind::boundary, false, {0}, {}}};
    specification.requirement_policies.erase(
      std::remove_if(specification.requirement_policies.begin(),
                     specification.requirement_policies.end(),
                     [](const RequirementPolicySpec &policy) {
                       return policy.kind == RequirementKind::fixed_dirichlet;
                     }),
      specification.requirement_policies.end());
    specification.requirement_policies.insert(
      specification.requirement_policies.begin(),
      {"state_mean_zero_gauge", "state", RequirementKind::mean_zero_multiplier,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "one mean-zero Lagrange multiplier in the state and adjoint solves",
       "domain"});
    return specification;
  }

  // This deliberately differs from the homogeneous reference graph above:
  // the state variable denotes independent coordinates and the declared map
  // reconstructs the physical field consumed by residuals and observations.
  inline ProblemSpec
  make_fixed_dirichlet_scalar_diffusion_reaction_problem(
    const bool with_cellwise_box = false)
  {
    ProblemSpec specification =
      make_scalar_diffusion_reaction_problem(with_cellwise_box);
    specification.id = "scalar_diffusion_reaction_fixed_dirichlet";
    specification.label = "Scalar diffusion-reaction with fixed Dirichlet lifting";
    specification.regions.at(1).label = "Fixed Dirichlet boundary";
    specification.data.push_back(
      {"fixed_dirichlet_data", "Fixed Dirichlet data", DataKind::function,
       DataRole::fixed_dirichlet_lifting, "state_space"});
    specification.transformations = {
      {"fixed_dirichlet_reconstruction", "Fixed Dirichlet reconstruction",
       TransformationKind::fixed_dirichlet_reconstruction, "state",
       "state_space", "fixed_dirichlet_data", ""}};
    specification.variables.at(0).physical_field_transform_id =
      "fixed_dirichlet_reconstruction";
    specification.requirement_policies.at(0).selected_policy =
      "P_h independent FE_Q coordinates plus nodal fixed Dirichlet lifting";
    return specification;
  }

  // P3.2 deliberately represents controlled essential data as a physical
  // state transformation, never as a boundary residual load. The first
  // registered lifting uses one continuous nodal trace coefficient for every
  // state DoF on the complete selected exterior boundary.
  inline ProblemSpec
  make_dirichlet_control_scalar_diffusion_reaction_problem()
  {
    ProblemSpec specification = make_scalar_diffusion_reaction_problem();
    specification.id = "scalar_diffusion_reaction_dirichlet_control";
    specification.label =
      "Scalar diffusion-reaction with Dirichlet control lifting";
    specification.regions.at(1) =
      {"control_boundary", "Complete controlled Dirichlet boundary",
       RegionKind::boundary, false, {0}, {}};
    specification.spaces.at(2) =
      {"control_space", "Nodal Dirichlet trace control", "control_boundary",
       SpaceTopology::h1, SpaceRole::control};
    specification.spaces.at(4) =
      {"control_observation_space", "Nodal Dirichlet trace observation",
       "control_boundary", SpaceTopology::h1, SpaceRole::observation};
    specification.pairings.at(2) =
      {"control_pairing", "Dirichlet trace control coefficient pairing",
       "control_space", "control_space"};
    specification.pairings.at(4) =
      {"control_observation_pairing",
       "Dirichlet trace observation coefficient pairing",
       "control_observation_space", "control_observation_space"};
    specification.variables.at(0).physical_field_transform_id =
      "dirichlet_control_lifting";
    specification.variables.at(1).label = "Dirichlet control";
    specification.residual_terms.erase(specification.residual_terms.begin() + 2);
    specification.equations.at(0).residual_term_ids =
      {"diffusion_reaction", "volume_source"};
    specification.observations.at(1) =
      {"control_boundary_restriction", "Dirichlet control restriction",
       ObservationKind::boundary_restriction, "control", "control_boundary",
       "control_observation_space", "control_observation_pairing"};
    specification.losses.at(1) =
      {"control_regularisation", "Quadratic Dirichlet control regularisation",
       LossKind::quadratic_control_regularisation,
       "control_boundary_restriction", "regularisation_weight",
       "control_observation_pairing"};
    specification.metrics.at(0) =
      {"control_l2_metric", "Dirichlet trace L2 metric", MetricKind::l2,
       "control", "control_pairing"};
    specification.transformations = {
      {"dirichlet_control_lifting", "Dirichlet control physical-state lifting",
       TransformationKind::dirichlet_control_lifting, "state", "state_space",
       "", "control"}};
    specification.requirement_policies.at(0) =
      {"state_controlled_dirichlet", "state",
       RequirementKind::controlled_dirichlet,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "complete-exterior-boundary nodal trace lifting with one shared coefficient per state boundary DoF",
       "control_boundary"};
    return specification;
  }
} // namespace nmopt::semantic::v1
