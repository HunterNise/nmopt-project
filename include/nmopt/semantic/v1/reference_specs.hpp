#pragma once

#include "nmopt/semantic/v1/types.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::semantic::v1
{
  namespace reference_detail
  {
    template <typename Component>
    Component &
    component_by_id(std::vector<Component> &components,
                    const std::string &     id,
                    const char *            component_name)
    {
      const auto component = std::find_if(
        components.begin(), components.end(), [&id](const Component &candidate) {
          return candidate.id == id;
        });
      const auto matches = std::count_if(
        components.begin(), components.end(), [&id](const Component &candidate) {
          return candidate.id == id;
        });
      if (matches != 1)
        throw std::logic_error("Reference specification requires exactly one " +
                               std::string(component_name) + " with id '" + id +
                               "'");
      return *component;
    }

    template <typename Component>
    void
    remove_component_by_id(std::vector<Component> &components,
                           const std::string &     id,
                           const char *            component_name)
    {
      (void)component_by_id(components, id, component_name);
      const auto component = std::find_if(
        components.begin(), components.end(), [&id](const Component &candidate) {
          return candidate.id == id;
        });
      components.erase(component);
    }
  } // namespace reference_detail

  // This is the current reference graph used to compare direct v0 assembly
  // with semantic compilation. It is a factory, not a PDE problem class.
  inline ProblemSpec
  make_scalar_diffusion_reaction_problem(const bool with_cellwise_box = false)
  {
    ProblemSpec specification;
    specification.id = "scalar_diffusion_reaction_volume_control";
    specification.label = "Scalar diffusion-reaction volume control";
    specification.regions = {
      {"domain", "Full volume domain", RegionKind::volume, true, {}, {}, {}},
      {"dirichlet_boundary", "Homogeneous Dirichlet boundary", RegionKind::boundary,
       false, {0}, {}, {}}};
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
       "state_observation_space", "state_observation_pairing", {}},
      {"control_observation", "Full-domain control restriction",
       ObservationKind::volume_restriction, "control", "domain",
       "control_observation_space", "control_observation_pairing", {}}};
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

  // P5.1's first registered target is a recombination of scalar residual
  // components. It deliberately keeps the existing volume-control,
  // observation, loss, metric, and optional cellwise-box declarations while
  // replacing the combined constant diffusion-reaction term by independently
  // bound coefficient Functions and adding one Robin boundary region.
  inline ProblemSpec
  make_general_scalar_elliptic_robin_problem(
    std::vector<unsigned int> fixed_dirichlet_boundary_ids,
    std::vector<unsigned int> robin_boundary_ids,
    const bool                with_cellwise_box = false)
  {
    ProblemSpec specification =
      make_scalar_diffusion_reaction_problem(with_cellwise_box);
    specification.id = "general_scalar_elliptic_robin_volume_control";
    specification.label =
      "General scalar elliptic volume control with Robin boundary";
    reference_detail::component_by_id(specification.regions,
                                      "dirichlet_boundary",
                                      "region")
      .boundary_ids = std::move(fixed_dirichlet_boundary_ids);
    specification.regions.push_back(
      {"robin_boundary", "Robin and transport-outflow boundary",
       RegionKind::boundary, false, std::move(robin_boundary_ids), {}, {}});

    reference_detail::component_by_id(specification.data,
                                      "diffusion",
                                      "data") =
      {"diffusion_tensor", "Tensor diffusion coefficient",
       DataKind::tensor_function, DataRole::diffusion, ""};
    reference_detail::component_by_id(specification.data,
                                      "reaction",
                                      "data") =
      {"reaction", "Reaction coefficient Function", DataKind::function,
       DataRole::reaction, ""};
    specification.data.push_back(
      {"conservative_transport", "Conservative transport coefficient",
       DataKind::vector_function, DataRole::conservative_transport, ""});
    specification.data.push_back(
      {"advective_transport", "Advective transport coefficient",
       DataKind::vector_function, DataRole::advective_transport, ""});
    specification.data.push_back(
      {"robin_coefficient", "Robin bilinear coefficient", DataKind::function,
       DataRole::robin_coefficient, ""});
    specification.data.push_back(
      {"robin_source", "Robin boundary source", DataKind::function,
       DataRole::robin_source, "state_test_space"});

    reference_detail::component_by_id(specification.residual_terms,
                                      "diffusion_reaction",
                                      "residual term") =
      {"tensor_diffusion", "Tensor diffusion", ResidualTermKind::tensor_diffusion,
       "state_equation", {"state"}, {"diffusion_tensor"}, ""};
    specification.residual_terms.push_back(
      {"conservative_transport", "Conservative transport",
       ResidualTermKind::conservative_transport, "state_equation", {"state"},
       {"conservative_transport"}, ""});
    specification.residual_terms.push_back(
      {"advective_transport", "Advective transport",
       ResidualTermKind::advective_transport, "state_equation", {"state"},
       {"advective_transport"}, ""});
    specification.residual_terms.push_back(
      {"reaction", "Reaction", ResidualTermKind::reaction, "state_equation",
       {"state"}, {"reaction"}, ""});
    specification.residual_terms.push_back(
      {"robin_bilinear", "Robin bilinear boundary term",
       ResidualTermKind::robin_bilinear, "state_equation", {"state"},
       {"robin_coefficient"}, "robin_boundary"});
    specification.residual_terms.push_back(
      {"robin_source", "Robin boundary source",
       ResidualTermKind::robin_source, "state_equation", {}, {"robin_source"},
       "robin_boundary"});
    reference_detail::component_by_id(specification.equations,
                                      "state_equation",
                                      "equation")
      .residual_term_ids = {"tensor_diffusion",
                            "conservative_transport",
                            "advective_transport",
                            "reaction",
                            "volume_source",
                            "volume_control",
                            "robin_bilinear",
                            "robin_source"};

    specification.requirement_policies.push_back(
      {"uniform_ellipticity_assumption", "diffusion_tensor",
       RequirementKind::uniform_ellipticity, RequirementStatus::user_assumed,
       RequirementScope::continuous_semantics,
       "bound tensor Function is uniformly elliptic", "domain"});
    specification.requirement_policies.push_back(
      {"coefficient_regularity_assumption", "state_equation",
       RequirementKind::coefficient_regularity,
       RequirementStatus::user_assumed,
       RequirementScope::continuous_semantics,
       "bound tensor, vector, scalar, and Robin coefficients have the declared bounded regularity",
       "domain"});
    specification.requirement_policies.push_back(
      {"coercivity_assumption", "state_equation", RequirementKind::coercivity,
       RequirementStatus::user_assumed,
       RequirementScope::continuous_semantics,
       "the composed scalar bilinear form is coercive on the fixed-Dirichlet state space",
       "domain"});
    specification.requirement_policies.push_back(
      {"scalar_boundary_partition", "state",
       RequirementKind::boundary_partition,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::both,
       "disjoint complete fixed-Dirichlet and Robin/transport-outflow regions; Neumann and separate transport-inflow regions are empty",
       ""});
    specification.requirement_policies.push_back(
      {"conormal_flux_convention", "robin_bilinear",
       RequirementKind::conormal_flux,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::both,
       "(A grad(y) - b y) dot n with outward unit normal", "robin_boundary"});
    specification.requirement_policies.push_back(
      {"robin_trace_policy", "robin_bilinear", RequirementKind::boundary_trace,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "FE_Q state trace with QGauss face quadrature", "robin_boundary"});
    specification.requirement_policies.push_back(
      {"transport_boundary_trace", "conservative_transport",
       RequirementKind::transport_boundary_trace,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::both,
       "Robin region is the selected natural transport outflow; remaining exterior faces are fixed Dirichlet",
       "robin_boundary"});
    return specification;
  }

  namespace reference_detail
  {
    inline void
    apply_subdomain_tracking_delta(ProblemSpec &       specification,
                                   const unsigned int observed_material_id)
    {
      specification.id = "scalar_diffusion_reaction_subdomain_tracking";
      specification.label = "Scalar diffusion-reaction with subdomain tracking";
      specification.regions.push_back(
        {"observation_subdomain", "Material subdomain observation region",
         RegionKind::volume, false, {}, {observed_material_id}, {}});
      component_by_id(specification.spaces,
                      "state_observation_space",
                      "space")
        .region_id = "observation_subdomain";
      auto &observation = component_by_id(
        specification.observations, "state_observation", "observation");
      observation.label = "Subdomain state restriction";
      observation.region_id = "observation_subdomain";
      component_by_id(specification.requirement_policies,
                      "desired_state_quadrature_policy",
                      "requirement policy")
        .region_id = "observation_subdomain";
    }
  } // namespace reference_detail

  inline ProblemSpec
  make_subdomain_tracking_scalar_diffusion_reaction_problem(
    const unsigned int observed_material_id,
    const bool         with_cellwise_box = false)
  {
    ProblemSpec specification =
      make_scalar_diffusion_reaction_problem(with_cellwise_box);
    reference_detail::apply_subdomain_tracking_delta(specification,
                                                     observed_material_id);
    return specification;
  }

  // P5.3's first registered target is the C5.10 finite point-sensor slice.
  // Coordinates are semantic data: the deal.II lowerer checks their mesh
  // dimension and owns the resulting immutable point-evaluation operator.
  inline ProblemSpec
  make_point_sensor_scalar_diffusion_reaction_problem(
    std::vector<std::vector<double>> sensor_coordinates)
  {
    ProblemSpec specification = make_scalar_diffusion_reaction_problem();
    specification.id = "scalar_diffusion_reaction_point_sensor";
    specification.label = "Scalar diffusion-reaction with point sensors";
    specification.regions.push_back(
      {"point_sensor_region", "Immutable point-sensor region",
       RegionKind::point_set, false, {}, {}, std::move(sensor_coordinates)});

    auto &sensor_region = reference_detail::component_by_id(
      specification.regions, "point_sensor_region", "region");
    auto &sensor_space = reference_detail::component_by_id(
      specification.spaces, "state_observation_space", "space");
    sensor_space = {"state_observation_space", "Finite point-sensor output",
                    "point_sensor_region", SpaceTopology::l2,
                    SpaceRole::observation, true,
                    sensor_region.point_coordinates.size()};
    reference_detail::component_by_id(specification.data,
                                      "desired_state",
                                      "data")
      .space_id = "state_observation_space";
    reference_detail::component_by_id(specification.observations,
                                      "state_observation",
                                      "observation") =
      {"state_observation", "Finite point-sensor state observation",
       ObservationKind::point_sensor, "state", "point_sensor_region",
       "state_observation_space", "state_observation_pairing", {}};
    reference_detail::component_by_id(specification.requirement_policies,
                                      "desired_state_quadrature_policy",
                                      "requirement policy") =
      {"desired_state_point_evaluation_policy", "desired_state",
       RequirementKind::analytic_quadrature_evaluation,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "analytic Function evaluated at the immutable sensor coordinates",
       "point_sensor_region"};
    specification.requirement_policies.push_back(
      {"point_sensor_transposition_policy", "state_equation",
       RequirementKind::transposition_formulation,
       RequirementStatus::provided, RequirementScope::continuous_semantics,
       "Y=H2(Omega) cap H1_0(Omega); T=-Delta+rI:Y->L2(Omega) is an isomorphism; the point residual is represented in Y* and its adjoint is very weak",
       "domain"});
    specification.requirement_policies.push_back(
      {"point_sensor_domain_regularity", "state_equation",
       RequirementKind::domain_regularity, RequirementStatus::user_assumed,
       RequirementScope::continuous_semantics,
       "the declared domain is convex or C2 so point evaluation on Y and the very-weak adjoint are meaningful",
       "domain"});
    specification.requirement_policies.push_back(
      {"point_sensor_evaluation_policy", "state_observation",
       RequirementKind::analytic_quadrature_evaluation,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "FE_Q shape functions are evaluated at immutable physical sensor coordinates; the transpose C_h^T is assembled as a finite-dimensional very-weak point load; no nearest-node or quadrature-point coincidence is used",
       "point_sensor_region"});
    return specification;
  }

  // P5.2's first observation target changes the state observation and its
  // pairing, not the residual or the control search metric. The selected
  // realization is the full H1_0 inner product assembled from mass and
  // stiffness contributions.
  inline ProblemSpec
  make_h1_state_tracking_scalar_diffusion_reaction_problem()
  {
    ProblemSpec specification = make_scalar_diffusion_reaction_problem();
    specification.id = "scalar_diffusion_reaction_h1_state_tracking";
    specification.label =
      "Scalar diffusion-reaction with H1 state tracking";
    reference_detail::component_by_id(specification.spaces,
                                      "state_observation_space",
                                      "space") =
      {"state_observation_space", "H1 state observation", "domain",
       SpaceTopology::h1, SpaceRole::observation};
    reference_detail::component_by_id(specification.pairings,
                                      "state_observation_pairing",
                                      "pairing") =
      {"state_observation_pairing", "H1_0 state-observation pairing",
       "state_observation_space", "state_observation_space"};
    reference_detail::component_by_id(specification.observations,
                                      "state_observation",
                                      "observation") =
      {"state_observation", "Full-domain H1 state restriction",
       ObservationKind::h1_state_restriction, "state", "domain",
       "state_observation_space", "state_observation_pairing", {}};
    reference_detail::component_by_id(specification.requirement_policies,
                                      "desired_state_quadrature_policy",
                                      "requirement policy")
      .selected_policy =
      "analytic Function value and gradient evaluated at selected volume quadrature";
    return specification;
  }

  namespace reference_detail
  {
    inline void
    apply_homogeneous_dirichlet_continuous_control_delta(
      ProblemSpec &specification)
    {
      component_by_id(specification.spaces, "control_space", "space") =
        {"control_space", "Homogeneous-Dirichlet continuous control", "domain",
         SpaceTopology::h1, SpaceRole::control};
      component_by_id(specification.pairings, "control_pairing", "pairing") =
        {"control_pairing", "Independent H1_0 control coefficient pairing",
         "control_space", "control_space"};
      specification.requirement_policies.push_back(
        {"control_homogeneous_dirichlet_policy", "control",
         RequirementKind::fixed_dirichlet,
         RequirementStatus::selected_discrete_realisation,
         RequirementScope::both,
         "P_h is continuous FE_Q with independent homogeneous-Dirichlet DoFs",
         "dirichlet_boundary"});
    }
  } // namespace reference_detail

  // Companion graph for P5.2's metric comparison. It changes the discrete
  // control search space but retains the L2 Riesz map.
  inline ProblemSpec
  make_l2_metric_h1_state_tracking_continuous_control_problem()
  {
    ProblemSpec specification =
      make_h1_state_tracking_scalar_diffusion_reaction_problem();
    specification.id =
      "scalar_diffusion_reaction_h1_state_tracking_continuous_control_l2_metric";
    specification.label =
      "Scalar diffusion-reaction with H1 state tracking and continuous-control L2 metric";
    reference_detail::apply_homogeneous_dirichlet_continuous_control_delta(
      specification);
    reference_detail::component_by_id(specification.metrics,
                                      "control_l2_metric",
                                      "metric")
      .label = "Continuous-control L2 metric";
    return specification;
  }

  // P5.2's negative metric is a separate search-geometry delta on the
  // companion energy-tracking graph. Its compiler realizes
  // G_h = M_h K_h^{-1} M_h and records the metric solve policy.
  inline ProblemSpec
  make_hminus1_metric_h1_state_tracking_scalar_diffusion_reaction_problem()
  {
    ProblemSpec specification =
      make_l2_metric_h1_state_tracking_continuous_control_problem();
    specification.id =
      "scalar_diffusion_reaction_h1_state_tracking_hminus1_metric";
    specification.label =
      "Scalar diffusion-reaction with H1 state tracking and H-1 control metric";
    reference_detail::component_by_id(specification.metrics,
                                      "control_l2_metric",
                                      "metric") =
      {"control_hminus1_metric", "Discrete H-1 control metric",
       MetricKind::hminus1, "control", "control_pairing"};
    specification.formulation.metric_id = "control_hminus1_metric";
    return specification;
  }

  // P2.3's first half changes the objective, not the search geometry.  The
  // control has continuous FE_Q coordinates and the declared loss is the
  // H1 norm, while the selected algorithmic metric remains L2.
  namespace reference_detail
  {
    inline void
    apply_h1_control_regularisation_delta(ProblemSpec &specification)
    {
      specification.id = "scalar_diffusion_reaction_h1_control_regularisation";
      specification.label =
        "Scalar diffusion-reaction with H1 control regularisation";
      component_by_id(specification.spaces, "control_space", "space") =
        {"control_space", "Continuous control", "domain", SpaceTopology::h1,
         SpaceRole::control};
      component_by_id(specification.spaces,
                      "control_observation_space",
                      "space") =
        {"control_observation_space", "Continuous control observation", "domain",
         SpaceTopology::h1, SpaceRole::observation};
      component_by_id(specification.pairings,
                      "control_pairing",
                      "pairing") =
        {"control_pairing", "Continuous control coefficient pairing",
         "control_space", "control_space"};
      component_by_id(specification.pairings,
                      "control_observation_pairing",
                      "pairing") =
        {"control_observation_pairing",
         "Continuous control observation pairing",
         "control_observation_space", "control_observation_space"};
      component_by_id(specification.observations,
                      "control_observation",
                      "observation") =
        {"control_observation", "Continuous control identity observation",
         ObservationKind::volume_restriction, "control", "domain",
         "control_observation_space", "control_observation_pairing", {}};
      component_by_id(specification.losses,
                      "control_regularisation",
                      "loss") =
        {"control_h1_regularisation", "Quadratic H1 control regularisation",
         LossKind::quadratic_h1_control_regularisation, "control_observation",
         "regularisation_weight", "control_observation_pairing"};
      component_by_id(specification.metrics, "control_l2_metric", "metric") =
        {"control_l2_metric", "Continuous-control L2 metric", MetricKind::l2,
         "control", "control_pairing"};
    }

    inline void
    apply_h1_metric_delta(ProblemSpec &specification)
    {
      specification.id = "scalar_diffusion_reaction_h1_control_metric";
      specification.label = "Scalar diffusion-reaction with H1 control metric";
      component_by_id(specification.metrics, "control_l2_metric", "metric") =
        {"control_h1_metric", "Continuous-control H1 metric", MetricKind::h1,
         "control", "control_pairing"};
      specification.formulation.metric_id = "control_h1_metric";
    }
  } // namespace reference_detail

  inline ProblemSpec
  make_h1_regularised_scalar_diffusion_reaction_problem()
  {
    ProblemSpec specification = make_scalar_diffusion_reaction_problem();
    reference_detail::apply_h1_control_regularisation_delta(specification);
    return specification;
  }

  // This changes only the Riesz map used for search directions. The inherited
  // H1 regularisation loss and every residual/objective action are unchanged.
  inline ProblemSpec
  make_h1_metric_scalar_diffusion_reaction_problem()
  {
    ProblemSpec specification =
      make_h1_regularised_scalar_diffusion_reaction_problem();
    reference_detail::apply_h1_metric_delta(specification);
    return specification;
  }

  // P3.1 uses the binary reduced DTO decision port for a physical diffusion
  // coefficient rather than a source control. Positivity is explicit through
  // the required cellwise box; a logarithmic parameterisation is deliberately
  // left to a later transformation realization.
  namespace reference_detail
  {
    inline void
    apply_coefficient_identification_delta(ProblemSpec &specification)
    {
      specification.id = "scalar_diffusion_reaction_coefficient_identification";
      specification.label =
        "Scalar diffusion-reaction coefficient identification";
      component_by_id(specification.spaces, "control_space", "space") =
        {"parameter_space", "Cellwise diffusion parameter", "domain",
         SpaceTopology::l2, SpaceRole::parameter};
      component_by_id(specification.spaces,
                      "control_observation_space",
                      "space") =
        {"parameter_observation_space", "Cellwise parameter observation",
         "domain", SpaceTopology::l2, SpaceRole::observation};
      component_by_id(specification.pairings,
                      "control_pairing",
                      "pairing") =
        {"parameter_pairing", "Parameter coefficient pairing", "parameter_space",
         "parameter_space"};
      component_by_id(specification.pairings,
                      "control_observation_pairing",
                      "pairing") =
        {"parameter_observation_pairing", "Parameter observation pairing",
         "parameter_observation_space", "parameter_observation_space"};
      component_by_id(specification.variables, "control", "variable") =
        {"diffusion_parameter", "Diffusion parameter", VariableRole::parameter,
         "parameter_space", ""};
      remove_component_by_id(specification.data, "diffusion", "data");
      specification.data.push_back(
        {"parameter_lower_bound", "Strictly positive diffusion lower bound",
         DataKind::cellwise_bound, DataRole::lower_bound, "parameter_space"});
      specification.data.push_back(
        {"parameter_upper_bound", "Diffusion upper bound",
         DataKind::cellwise_bound, DataRole::upper_bound, "parameter_space"});
      component_by_id(specification.residual_terms,
                      "diffusion_reaction",
                      "residual term") =
        {"parameter_diffusion_reaction", "Parameter diffusion and reaction",
         ResidualTermKind::parameter_diffusion_reaction, "state_equation",
         {"state", "diffusion_parameter"}, {"reaction"}, ""};
      remove_component_by_id(specification.residual_terms,
                             "volume_control",
                             "residual term");
      component_by_id(specification.equations,
                      "state_equation",
                      "equation")
        .residual_term_ids = {"parameter_diffusion_reaction", "volume_source"};
      component_by_id(specification.observations,
                      "control_observation",
                      "observation") =
        {"parameter_observation", "Full-domain parameter restriction",
         ObservationKind::volume_restriction, "diffusion_parameter", "domain",
         "parameter_observation_space", "parameter_observation_pairing", {}};
      component_by_id(specification.losses,
                      "control_regularisation",
                      "loss") =
        {"parameter_regularisation", "Quadratic parameter regularisation",
         LossKind::quadratic_parameter_regularisation, "parameter_observation",
         "regularisation_weight", "parameter_observation_pairing"};
      component_by_id(specification.metrics, "control_l2_metric", "metric") =
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
    }
  } // namespace reference_detail

  inline ProblemSpec
  make_coefficient_identification_problem()
  {
    ProblemSpec specification = make_scalar_diffusion_reaction_problem();
    reference_detail::apply_coefficient_identification_delta(specification);
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
      {"domain", "Full volume domain", RegionKind::volume, true, {}, {}, {}},
      {"dirichlet_boundary", "Homogeneous Dirichlet boundary", RegionKind::boundary,
       false, {0}, {}, {}},
      {"control_boundary", "Neumann control boundary", RegionKind::boundary,
       false, {1}, {}, {}},
      {"observation_boundary", "Boundary tracking region", RegionKind::boundary,
       false, {2}, {}, {}}};
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
       "state_observation_pairing", {}},
      {"control_boundary_restriction", "Boundary control restriction",
       ObservationKind::boundary_restriction, "control", "control_boundary",
       "control_observation_space", "control_observation_pairing", {}}};
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

  // C5.6 composes the natural boundary control with a volume restriction and
  // conservative transport.  It deliberately reuses the Neumann residual
  // and facewise control coordinate declarations: the observation and the
  // scalar state operator are the only changed components.
  inline ProblemSpec
  make_neumann_convection_subdomain_tracking_problem(
    const unsigned int observed_material_id,
    const bool         with_facewise_box = false)
  {
    ProblemSpec specification =
      make_neumann_boundary_control_problem(with_facewise_box);
    specification.id = "scalar_convection_neumann_subdomain_control";
    specification.label =
      "Scalar conservative transport with Neumann control and subdomain tracking";
    specification.regions.erase(
      std::remove_if(specification.regions.begin(), specification.regions.end(),
                     [](const RegionSpec &region) {
                       return region.id == "observation_boundary";
                     }),
      specification.regions.end());
    specification.regions.push_back(
      {"observation_subdomain", "Material subdomain observation region",
       RegionKind::volume, false, {}, {observed_material_id}, {}});

    reference_detail::component_by_id(specification.spaces,
                                      "state_observation_space",
                                      "space") =
      {"state_observation_space", "Subdomain state observation",
       "observation_subdomain", SpaceTopology::l2, SpaceRole::observation};
    reference_detail::component_by_id(specification.observations,
                                      "state_boundary_trace",
                                      "observation") =
      {"state_subdomain_restriction", "Subdomain state restriction",
       ObservationKind::volume_restriction, "state", "observation_subdomain",
       "state_observation_space", "state_observation_pairing", {}};
    reference_detail::component_by_id(specification.losses,
                                      "state_tracking",
                                      "loss")
      .source_observation_id = "state_subdomain_restriction";
    reference_detail::component_by_id(specification.data,
                                      "desired_state",
                                      "data")
      .label = "Desired subdomain state";
    reference_detail::component_by_id(specification.requirement_policies,
                                      "state_boundary_trace_policy",
                                      "requirement policy") =
      {"state_subdomain_restriction_policy", "state_subdomain_restriction",
       RequirementKind::analytic_quadrature_evaluation,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "FE_Q state restriction assembled on the selected material cells",
       "observation_subdomain"};
    reference_detail::component_by_id(specification.requirement_policies,
                                      "desired_state_quadrature_policy",
                                      "requirement policy") =
      {"desired_state_quadrature_policy", "desired_state",
       RequirementKind::analytic_quadrature_evaluation,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "analytic Function evaluated at selected volume quadrature on the observed material cells",
       "observation_subdomain"};

    specification.data.push_back(
      {"conservative_transport", "Conservative transport coefficient",
       DataKind::vector_function, DataRole::conservative_transport, ""});
    specification.residual_terms.insert(
      specification.residual_terms.begin() + 1,
      {"conservative_transport", "Conservative transport",
       ResidualTermKind::conservative_transport, "state_equation", {"state"},
       {"conservative_transport"}, ""});
    reference_detail::component_by_id(specification.equations,
                                      "state_equation",
                                      "equation")
      .residual_term_ids =
      {"diffusion_reaction", "conservative_transport", "volume_source",
       "neumann_control"};
    specification.requirement_policies.push_back(
      {"convection_coercivity_assumption", "state_equation",
       RequirementKind::coercivity, RequirementStatus::user_assumed,
       RequirementScope::continuous_semantics,
       "the selected conservative transport satisfies the C5.6 coercivity assumptions on the fixed-Dirichlet state space",
       "domain"});
    specification.requirement_policies.push_back(
      {"convection_regularity_assumption", "conservative_transport",
       RequirementKind::coefficient_regularity, RequirementStatus::user_assumed,
       RequirementScope::continuous_semantics,
       "transport components are Lipschitz on the volume domain",
       "domain"});
    specification.requirement_policies.push_back(
      {"convection_neumann_sign_assumption", "conservative_transport",
       RequirementKind::transport_boundary_trace,
       RequirementStatus::user_assumed,
       RequirementScope::continuous_semantics,
       "the conservative transport satisfies b dot n <= 0 on the Neumann-control boundary",
       "control_boundary"});
    specification.requirement_policies.push_back(
      {"neumann_convection_partition", "state",
       RequirementKind::boundary_partition,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::both,
       "disjoint complete fixed-Dirichlet and Neumann-control boundary regions; the Neumann datum is the conormal flux of the conservative transport form",
       ""});
    return specification;
  }

  namespace reference_detail
  {
    inline void
    apply_weighted_boundary_trace_delta(ProblemSpec &specification)
    {
      specification.id =
        "scalar_diffusion_reaction_weighted_boundary_trace_control";
      specification.label =
        "Scalar diffusion-reaction control with weighted boundary trace";

      auto &observation = component_by_id(specification.observations,
                                          "state_boundary_trace",
                                          "observation");
      observation.id = "weighted_state_boundary_trace";
      observation.label = "Weighted boundary state trace";
      observation.kind = ObservationKind::weighted_boundary_trace;
      observation.data_ids = {"boundary_weight"};

      component_by_id(specification.losses, "state_tracking", "loss")
        .source_observation_id = observation.id;
      component_by_id(specification.requirement_policies,
                      "state_boundary_trace_policy",
                      "requirement policy")
        .subject_id = observation.id;

      specification.data.push_back(
        {"boundary_weight", "Boundary observation weight", DataKind::function,
         DataRole::observation_weight, "state_observation_space"});
      specification.requirement_policies.push_back(
        {"boundary_weight_quadrature_policy", "boundary_weight",
         RequirementKind::analytic_quadrature_evaluation,
         RequirementStatus::selected_discrete_realisation,
         RequirementScope::discrete_compilation,
         "analytic Function evaluated at selected boundary face quadrature",
         "observation_boundary"});
      specification.requirement_policies.push_back(
        {"boundary_weight_boundedness", "boundary_weight",
         RequirementKind::coefficient_regularity,
         RequirementStatus::user_assumed,
         RequirementScope::continuous_semantics,
         "boundary weight belongs to L-infinity on the observation boundary",
         "observation_boundary"});
    }
  } // namespace reference_detail

  inline ProblemSpec
  make_weighted_boundary_trace_neumann_control_problem(
    const bool with_facewise_box = false)
  {
    ProblemSpec specification =
      make_neumann_boundary_control_problem(with_facewise_box);
    reference_detail::apply_weighted_boundary_trace_delta(specification);
    return specification;
  }

  // The pure-Neumann variant keeps the natural boundary residual of the
  // preceding graph but replaces the fixed-boundary uniqueness policy by the
  // selected discrete mean constraint.  The auxiliary multiplier belongs to
  // the compiled solve, not to the user-facing state/control graph.
  namespace reference_detail
  {
    inline void
    apply_pure_neumann_delta(ProblemSpec &specification)
    {
      specification.id = "scalar_diffusion_pure_neumann_boundary_control";
      specification.label =
        "Scalar diffusion pure-Neumann boundary control with mean constraint";
      specification.regions = {
        {"domain", "Full volume domain", RegionKind::volume, true, {}, {}, {}},
        {"control_boundary", "Pure-Neumann control boundary",
         RegionKind::boundary, false, {0}, {}, {}},
        {"observation_boundary", "Pure-Neumann boundary tracking region",
         RegionKind::boundary, false, {0}, {}, {}}};
      remove_component_by_id(specification.requirement_policies,
                             "state_fixed_dirichlet",
                             "requirement policy");
      specification.requirement_policies.insert(
        specification.requirement_policies.begin(),
        {"state_mean_zero_gauge", "state",
         RequirementKind::mean_zero_multiplier,
         RequirementStatus::selected_discrete_realisation,
         RequirementScope::discrete_compilation,
         "one mean-zero Lagrange multiplier in the state and adjoint solves",
         "domain"});
    }
  } // namespace reference_detail

  inline ProblemSpec
  make_pure_neumann_boundary_control_problem()
  {
    ProblemSpec specification = make_neumann_boundary_control_problem();
    reference_detail::apply_pure_neumann_delta(specification);
    return specification;
  }

  // This deliberately differs from the homogeneous reference graph above:
  // the state variable denotes independent coordinates and the declared map
  // reconstructs the physical field consumed by residuals and observations.
  namespace reference_detail
  {
    inline void
    apply_fixed_dirichlet_reconstruction_delta(ProblemSpec &specification)
    {
      specification.id = "scalar_diffusion_reaction_fixed_dirichlet";
      specification.label =
        "Scalar diffusion-reaction with fixed Dirichlet lifting";
      component_by_id(specification.regions,
                      "dirichlet_boundary",
                      "region")
        .label = "Fixed Dirichlet boundary";
      specification.data.push_back(
        {"fixed_dirichlet_data", "Fixed Dirichlet data", DataKind::function,
         DataRole::fixed_dirichlet_lifting, "state_space"});
      specification.transformations = {
        {"fixed_dirichlet_reconstruction", "Fixed Dirichlet reconstruction",
         TransformationKind::fixed_dirichlet_reconstruction, "state",
         "state_space", "fixed_dirichlet_data", ""}};
      component_by_id(specification.variables, "state", "variable")
        .physical_field_transform_id = "fixed_dirichlet_reconstruction";
      component_by_id(specification.requirement_policies,
                      "state_fixed_dirichlet",
                      "requirement policy")
        .selected_policy =
        "P_h independent FE_Q coordinates plus nodal fixed Dirichlet lifting";
    }

    inline void
    apply_dirichlet_control_delta(ProblemSpec &specification)
    {
      specification.id = "scalar_diffusion_reaction_dirichlet_control";
      specification.label =
        "Scalar diffusion-reaction with Dirichlet control lifting";
      component_by_id(specification.regions,
                      "dirichlet_boundary",
                      "region") =
        {"control_boundary", "Complete controlled Dirichlet boundary",
         RegionKind::boundary, false, {0}, {}, {}};
      component_by_id(specification.spaces, "control_space", "space") =
        {"control_space", "Nodal Dirichlet trace control", "control_boundary",
         SpaceTopology::h1, SpaceRole::control};
      component_by_id(specification.spaces,
                      "control_observation_space",
                      "space") =
        {"control_observation_space", "Nodal Dirichlet trace observation",
         "control_boundary", SpaceTopology::h1, SpaceRole::observation};
      component_by_id(specification.pairings,
                      "control_pairing",
                      "pairing") =
        {"control_pairing", "Dirichlet trace control coefficient pairing",
         "control_space", "control_space"};
      component_by_id(specification.pairings,
                      "control_observation_pairing",
                      "pairing") =
        {"control_observation_pairing",
         "Dirichlet trace observation coefficient pairing",
         "control_observation_space", "control_observation_space"};
      component_by_id(specification.variables, "state", "variable")
        .physical_field_transform_id = "dirichlet_control_lifting";
      component_by_id(specification.variables, "control", "variable").label =
        "Dirichlet control";
      remove_component_by_id(specification.residual_terms,
                             "volume_control",
                             "residual term");
      component_by_id(specification.equations,
                      "state_equation",
                      "equation")
        .residual_term_ids = {"diffusion_reaction", "volume_source"};
      component_by_id(specification.observations,
                      "control_observation",
                      "observation") =
        {"control_boundary_restriction", "Dirichlet control restriction",
         ObservationKind::boundary_restriction, "control", "control_boundary",
         "control_observation_space", "control_observation_pairing", {}};
      component_by_id(specification.losses,
                      "control_regularisation",
                      "loss") =
        {"control_regularisation", "Quadratic Dirichlet control regularisation",
         LossKind::quadratic_control_regularisation,
         "control_boundary_restriction", "regularisation_weight",
         "control_observation_pairing"};
      component_by_id(specification.metrics, "control_l2_metric", "metric") =
        {"control_l2_metric", "Dirichlet trace L2 metric", MetricKind::l2,
         "control", "control_pairing"};
      specification.transformations = {
        {"dirichlet_control_lifting",
         "Dirichlet control physical-state lifting",
         TransformationKind::dirichlet_control_lifting, "state", "state_space",
         "", "control"}};
      component_by_id(specification.requirement_policies,
                      "state_fixed_dirichlet",
                      "requirement policy") =
        {"state_controlled_dirichlet", "state",
         RequirementKind::controlled_dirichlet,
         RequirementStatus::selected_discrete_realisation,
         RequirementScope::discrete_compilation,
         "complete-exterior-boundary nodal trace lifting with one shared coefficient per state boundary DoF",
         "control_boundary"};
    }

    inline void
    apply_normalized_dirichlet_laplace_delta(ProblemSpec &specification)
    {
      remove_component_by_id(specification.data, "diffusion", "data");
      remove_component_by_id(specification.data, "reaction", "data");
      specification.spaces.push_back(
        {"forcing_space", "L2 volume forcing", "domain", SpaceTopology::l2,
         SpaceRole::data});
      component_by_id(specification.data, "forcing", "data").space_id =
        "forcing_space";
      component_by_id(specification.residual_terms,
                      "diffusion_reaction",
                      "residual term") =
        {"laplacian", "Normalized Laplace bilinear form",
         ResidualTermKind::laplacian, "state_equation", {"state"}, {}, ""};
      component_by_id(specification.equations,
                      "state_equation",
                      "equation")
        .residual_term_ids = {"laplacian", "volume_source"};
      component_by_id(specification.requirement_policies,
                      "state_controlled_dirichlet",
                      "requirement policy") =
        {"state_controlled_dirichlet", "state",
         RequirementKind::controlled_dirichlet,
         RequirementStatus::selected_discrete_realisation,
         RequirementScope::both,
         "complete-boundary H1 variational Dirichlet control lowered by the explicit conforming nodal trace lifting",
         "control_boundary"};
    }

    inline void
    apply_hhalf_dirichlet_control_delta(ProblemSpec &specification)
    {
      component_by_id(specification.spaces, "control_space", "space") =
        {"control_space", "H1/2 Dirichlet boundary control",
         "control_boundary", SpaceTopology::hhalf, SpaceRole::control};
      component_by_id(specification.spaces,
                      "control_observation_space",
                      "space") =
        {"control_observation_space", "H1/2 boundary-control observation",
         "control_boundary", SpaceTopology::hhalf, SpaceRole::observation};
      component_by_id(specification.pairings,
                      "control_pairing",
                      "pairing") =
        {"control_pairing", "H1/2 trace control coefficient pairing",
         "control_space", "control_space"};
      component_by_id(specification.pairings,
                      "control_observation_pairing",
                      "pairing") =
        {"control_observation_pairing",
         "H1/2 trace observation coefficient pairing",
         "control_observation_space", "control_observation_space"};
      component_by_id(specification.losses,
                      "control_regularisation",
                      "loss") =
        {"control_regularisation", "Quadratic H1/2 control regularisation",
         LossKind::quadratic_hhalf_control_regularisation,
         "control_boundary_restriction", "regularisation_weight",
         "control_observation_pairing"};
      component_by_id(specification.metrics,
                      "control_l2_metric",
                      "metric") =
        {"control_hhalf_metric", "Minimum-extension H1/2 trace metric",
         MetricKind::hhalf, "control", "control_pairing"};
      specification.formulation.metric_id = "control_hhalf_metric";
      specification.requirement_policies.push_back(
        {"control_hhalf_metric_realisation", "control_hhalf_metric",
         RequirementKind::fractional_trace_realisation,
         RequirementStatus::selected_discrete_realisation,
         RequirementScope::discrete_compilation,
         "Schur complement of the volume H1 mass-plus-stiffness matrix under minimum-energy trace extension",
         "control_boundary"});
      specification.requirement_policies.push_back(
        {"discrete_conormal_policy", "state_equation",
         RequirementKind::conormal_flux,
         RequirementStatus::selected_discrete_realisation,
         RequirementScope::discrete_compilation,
         "outward discrete conormal is the lifting pullback of the adjoint residual",
         "control_boundary"});
    }
  } // namespace reference_detail

  inline ProblemSpec
  make_fixed_dirichlet_scalar_diffusion_reaction_problem(
    const bool with_cellwise_box = false)
  {
    ProblemSpec specification =
      make_scalar_diffusion_reaction_problem(with_cellwise_box);
    reference_detail::apply_fixed_dirichlet_reconstruction_delta(specification);
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
    reference_detail::apply_dirichlet_control_delta(specification);
    return specification;
  }

  // Section 5.11.1, option 1: the state uses the ordinary H1 variational
  // formulation and both the control loss and search metric use the selected
  // minimum-extension H1/2 trace Riesz map.
  inline ProblemSpec
  make_hhalf_dirichlet_laplace_control_problem()
  {
    ProblemSpec specification =
      make_dirichlet_control_scalar_diffusion_reaction_problem();
    specification.id = "hhalf_dirichlet_laplace_control";
    specification.label =
      "H1/2 Dirichlet Laplace control with fractional regularisation";
    reference_detail::apply_normalized_dirichlet_laplace_delta(specification);
    reference_detail::apply_hhalf_dirichlet_control_delta(specification);
    return specification;
  }

  // Section 5.11.1, option 2: the control remains in H1/2 and its search
  // metric is fractional, while the objective combines H1 state tracking
  // with an L2 boundary-control loss.
  inline ProblemSpec
  make_h1_tracking_hhalf_dirichlet_laplace_control_problem()
  {
    ProblemSpec specification = make_hhalf_dirichlet_laplace_control_problem();
    specification.id = "h1_tracking_hhalf_dirichlet_laplace_control";
    specification.label =
      "H1-tracking H1/2 Dirichlet Laplace control with L2 regularisation";
    reference_detail::component_by_id(specification.spaces,
                                      "state_observation_space",
                                      "space") =
      {"state_observation_space", "H1 state observation", "domain",
       SpaceTopology::h1, SpaceRole::observation};
    reference_detail::component_by_id(specification.pairings,
                                      "state_observation_pairing",
                                      "pairing") =
      {"state_observation_pairing", "H1 state-observation pairing",
       "state_observation_space", "state_observation_space"};
    reference_detail::component_by_id(specification.observations,
                                      "state_observation",
                                      "observation") =
      {"state_observation", "Full-domain H1 state restriction",
       ObservationKind::h1_state_restriction, "state", "domain",
       "state_observation_space", "state_observation_pairing", {}};
    reference_detail::component_by_id(specification.spaces,
                                      "control_observation_space",
                                      "space") =
      {"control_observation_space", "L2 boundary-control observation",
       "control_boundary", SpaceTopology::l2, SpaceRole::observation};
    reference_detail::component_by_id(specification.pairings,
                                      "control_observation_pairing",
                                      "pairing") =
      {"control_observation_pairing", "L2 boundary-control pairing",
       "control_observation_space", "control_observation_space"};
    reference_detail::component_by_id(specification.losses,
                                      "control_regularisation",
                                      "loss") =
      {"control_regularisation", "Quadratic L2 control regularisation",
       LossKind::quadratic_control_regularisation,
       "control_boundary_restriction", "regularisation_weight",
       "control_observation_pairing"};
    reference_detail::component_by_id(specification.requirement_policies,
                                      "desired_state_quadrature_policy",
                                      "requirement policy")
      .selected_policy =
      "analytic Function value and gradient evaluated at selected volume quadrature";
    return specification;
  }

  // Section 5.11.3: both the control loss and search metric use the
  // mass-plus-tangential-stiffness H1 trace Riesz map.
  inline ProblemSpec
  make_h1_dirichlet_laplace_control_problem()
  {
    ProblemSpec specification = make_hhalf_dirichlet_laplace_control_problem();
    specification.id = "h1_dirichlet_laplace_control";
    specification.label =
      "Tangential H1 Dirichlet Laplace control";
    reference_detail::component_by_id(specification.spaces,
                                      "control_space",
                                      "space") =
      {"control_space", "H1 Dirichlet boundary control", "control_boundary",
       SpaceTopology::h1, SpaceRole::control};
    reference_detail::component_by_id(specification.spaces,
                                      "control_observation_space",
                                      "space") =
      {"control_observation_space", "H1 boundary-control observation",
       "control_boundary", SpaceTopology::h1, SpaceRole::observation};
    reference_detail::component_by_id(specification.pairings,
                                      "control_pairing",
                                      "pairing") =
      {"control_pairing", "Tangential H1 trace control pairing",
       "control_space", "control_space"};
    reference_detail::component_by_id(specification.pairings,
                                      "control_observation_pairing",
                                      "pairing") =
      {"control_observation_pairing", "Tangential H1 observation pairing",
       "control_observation_space", "control_observation_space"};
    reference_detail::component_by_id(specification.losses,
                                      "control_regularisation",
                                      "loss") =
      {"control_regularisation", "Quadratic tangential H1 regularisation",
       LossKind::quadratic_h1_control_regularisation,
       "control_boundary_restriction", "regularisation_weight",
       "control_observation_pairing"};
    reference_detail::component_by_id(specification.metrics,
                                      "control_hhalf_metric",
                                      "metric") =
      {"control_h1_metric", "Tangential H1 trace metric", MetricKind::h1,
       "control", "control_pairing"};
    specification.formulation.metric_id = "control_h1_metric";
    reference_detail::remove_component_by_id(
      specification.requirement_policies,
      "control_hhalf_metric_realisation",
      "requirement policy");
    specification.requirement_policies.push_back(
      {"control_h1_metric_realisation", "control_h1_metric",
       RequirementKind::tangential_gradient_realisation,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "boundary mass plus tangential stiffness assembled on the complete controlled trace",
       "control_boundary"});
    return specification;
  }

  // Chapter 5.11.2 is a different continuous residual, not a relabeling of
  // the H1 diffusion-reaction graph. The state and control are L2 fields, the
  // test space is H2 cap H1_0, and the boundary datum enters through the
  // normal derivative of the test function. Its registered deal.II lowerer
  // may use the equivalent variational problem only because it separately
  // selects a conforming trace-control subspace.
  inline ProblemSpec
  make_l2_dirichlet_laplace_control_problem()
  {
    ProblemSpec specification =
      make_dirichlet_control_scalar_diffusion_reaction_problem();
    specification.id = "l2_dirichlet_laplace_control";
    specification.label = "L2 Dirichlet Laplace control by transposition";

    reference_detail::component_by_id(specification.spaces,
                                      "state_space",
                                      "space") =
      {"state_space", "L2 very-weak state", "domain", SpaceTopology::l2,
       SpaceRole::state};
    reference_detail::component_by_id(specification.spaces,
                                      "state_test_space",
                                      "space") =
      {"state_test_space", "H2 cap H1_0 transposition test", "domain",
       SpaceTopology::h2, SpaceRole::test};
    reference_detail::component_by_id(specification.spaces,
                                      "control_space",
                                      "space") =
      {"control_space", "L2 Dirichlet boundary control",
       "control_boundary", SpaceTopology::l2, SpaceRole::control};
    reference_detail::component_by_id(specification.spaces,
                                      "control_observation_space",
                                      "space") =
      {"control_observation_space", "L2 boundary-control observation",
       "control_boundary", SpaceTopology::l2, SpaceRole::observation};
    specification.spaces.push_back(
      {"forcing_space", "L2 volume forcing", "domain", SpaceTopology::l2,
       SpaceRole::data});

    reference_detail::component_by_id(specification.data, "forcing", "data")
      .space_id = "forcing_space";
    reference_detail::remove_component_by_id(specification.data,
                                              "diffusion",
                                              "data");
    reference_detail::remove_component_by_id(specification.data,
                                              "reaction",
                                              "data");

    reference_detail::component_by_id(specification.variables,
                                      "state",
                                      "variable")
      .physical_field_transform_id.clear();
    specification.transformations.clear();
    specification.residual_terms = {
      {"transposition_state_action", "Very-weak Laplace state action",
       ResidualTermKind::transposition_laplacian, "state_equation", {"state"},
       {}, ""},
      {"volume_source", "Volume source", ResidualTermKind::volume_source,
       "state_equation", {}, {"forcing"}, ""},
      {"dirichlet_transposition_control",
       "Dirichlet datum paired with outward normal test derivative",
       ResidualTermKind::dirichlet_transposition_control, "state_equation",
       {"control"}, {}, "control_boundary"}};
    reference_detail::component_by_id(specification.equations,
                                      "state_equation",
                                      "equation")
      .residual_term_ids = {"transposition_state_action",
                            "volume_source",
                            "dirichlet_transposition_control"};

    specification.requirement_policies = {
      {"state_controlled_dirichlet", "state",
       RequirementKind::controlled_dirichlet,
       RequirementStatus::selected_discrete_realisation, RequirementScope::both,
       "the L2 boundary datum enters the declared transposition residual; the registered discrete realization uses the complete-boundary conforming trace subspace",
       "control_boundary"},
      {"transposition_formulation", "state_equation",
       RequirementKind::transposition_formulation, RequirementStatus::provided,
       RequirementScope::continuous_semantics,
       "E(y,u;f)(psi)=(y,-Delta psi)_Omega-(f,psi)_Omega+(u,partial_n psi)_Gamma on H2 cap H1_0 tests",
       "domain"},
      {"transposition_domain_regularity", "state_equation",
       RequirementKind::domain_regularity, RequirementStatus::user_assumed,
       RequirementScope::continuous_semantics,
       "the Dirichlet Laplacian maps H2 cap H1_0 isomorphically to L2 on the declared convex or C2 domain; C2 is required before claiming the Proposition 5.16 optimal regularity bootstrap",
       "domain"},
      {"conforming_trace_subspace", "control",
       RequirementKind::conforming_trace_subspace,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "U_h=trace(V_h) is contained in H1/2(Gamma) and is lowered by the equivalent lifted H1 Galerkin state problem",
       "control_boundary"},
      {"discrete_conormal_policy", "state_equation",
       RequirementKind::conormal_flux,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "outward discrete conormal is the lifting pullback of the adjoint residual, never a pointwise FE_Q normal derivative",
       "control_boundary"},
      {"desired_state_quadrature_policy", "desired_state",
       RequirementKind::analytic_quadrature_evaluation,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "analytic Function evaluated at selected volume quadrature", "domain"}};
    return specification;
  }

  // P5.4's first target keeps the nodal L2 trace metric but makes the
  // fixed/controlled boundary partition and its interface ownership
  // explicit. At shared corner DoFs, fixed data owns the value, so the
  // control is the relative-interior trace with zero endpoint extension.
  inline ProblemSpec
  make_partial_dirichlet_control_scalar_diffusion_reaction_problem()
  {
    ProblemSpec specification =
      make_dirichlet_control_scalar_diffusion_reaction_problem();
    specification.id = "scalar_diffusion_reaction_partial_dirichlet_control";
    specification.label =
      "Scalar diffusion-reaction with partial Dirichlet control and fixed lifting";
    reference_detail::component_by_id(specification.regions,
                                      "control_boundary",
                                      "region")
      .boundary_ids = {1};
    specification.regions.push_back(
      {"fixed_dirichlet_boundary", "Fixed nonzero Dirichlet boundary",
       RegionKind::boundary, false, {0}, {}, {}});
    specification.data.push_back(
      {"fixed_dirichlet_data", "Fixed Dirichlet data", DataKind::function,
       DataRole::fixed_dirichlet_lifting, "state_space"});
    reference_detail::component_by_id(specification.transformations,
                                      "dirichlet_control_lifting",
                                      "transformation")
      .fixed_data_id = "fixed_dirichlet_data";
    specification.requirement_policies.push_back(
      {"state_fixed_dirichlet", "dirichlet_control_lifting", RequirementKind::fixed_dirichlet,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "nodal fixed-data lifting on the declared fixed boundary",
       "fixed_dirichlet_boundary"});
    specification.requirement_policies.push_back(
      {"partial_dirichlet_boundary_partition", "state",
       RequirementKind::boundary_partition,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::both,
       "complete disjoint fixed and controlled Dirichlet boundary partition",
       ""});
    specification.requirement_policies.push_back(
      {"partial_dirichlet_interface_policy", "dirichlet_control_lifting",
       RequirementKind::controlled_dirichlet,
       RequirementStatus::selected_discrete_realisation,
       RequirementScope::discrete_compilation,
       "fixed-data precedence at every fixed/controlled corner or interface DoF; the controlled nodal trace has zero endpoint extension",
       "control_boundary"});
    return specification;
  }
} // namespace nmopt::semantic::v1
