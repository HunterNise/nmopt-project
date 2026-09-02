#pragma once

#include "nmopt/application/chapter6.hpp"
#include "nmopt/application/dealii/chapter6_b2.hpp"
#include "parameter_binding.hpp"

#include <stdexcept>
#include <string>

namespace nmopt::application::runner::binding
{
  inline nmopt::application::chapter6::B2ObservationRegion
  b2_observation_region_from_combination(
    const ParameterCombination &combination)
  {
    const auto region = combination.values.find("observation-region");
    if (region != combination.values.end())
      {
        using namespace nmopt::application::chapter6;
        if (region->second == "wings")
          return B2ObservationRegion::wings;
        if (region->second == "full")
          return B2ObservationRegion::full;
      }
    throw std::invalid_argument(
      "B2 combinations need an observation-region axis with value wings or full");
  }

  inline void
  bind_b1_scenario(
    nmopt::application::chapter6::B1Scenario &scenario,
    const ParameterFile &file,
    const ParameterCombination &combination,
    const std::string &method_id,
    const double beta,
    const std::string &scenario_id)
  {
    using namespace nmopt::application;
    require_parameter(file,
                      "Benchmark/id",
                      "chapter-6.b1.distributed-laplace");
    require_parameter(file,
                      "Benchmark/recipe",
                      chapter6::b1_recipe_id);
    require_parameter(file, "Problem/observation", "full-domain");
    require_parameter(file,
                      "Output/selected fields",
                      "state, control, adjoint, negative-adjoint, target, forcing");
    if (combination.values.size() < 2)
      throw std::invalid_argument("B1 needs method and regularisation matrix axes");

    const auto control_representation =
      file.value("Problem/control representation");
    if (control_representation == "cellwise-volume")
      scenario.problem.recipe.discretisation =
        chapter5::DistributedControlDiscretisation::cellwise_constant;
    else if (control_representation ==
             "continuous-volume-homogeneous-dirichlet")
      scenario.problem.recipe.discretisation = chapter5::
        DistributedControlDiscretisation::homogeneous_dirichlet_continuous;
    else
      throw std::invalid_argument(
        "B1 has an unknown control representation '" +
        control_representation + "'");
    scenario.problem.recipe.with_cellwise_box =
      parameter_bool(file, "Problem/cellwise box constraint");
    scenario.problem.data.diffusion = parameter_double(file, "Runtime/diffusion");
    scenario.problem.data.reaction = parameter_double(file, "Runtime/reaction");
    scenario.problem.data.regularisation_weight = beta;
    scenario.problem.desired_state =
      parameter_scalar_function_definition(file, "Functions/desired state");
    scenario.problem.data.desired_state_provenance =
      scenario.problem.desired_state.provenance;
    scenario.problem.regularisation_sweep = {beta};
    scenario.problem.forcing = parameter_scalar_function_definition(
      file, "Functions/forcing");
    scenario.problem.data.forcing_provenance =
      scenario.problem.forcing.provenance;
    scenario.solver.method = parse_method(method_id);
    apply_common_parameter_options(scenario, file, scenario_id);
    apply_solver_options(scenario.solver, file, method_id);
    if (scenario.solver.method == chapter6::ReducedMethod::bfgs)
      throw std::invalid_argument("B1 does not support full BFGS");
  }

  inline void
  bind_b2_scenario(
    nmopt::application::chapter6::B2Scenario &scenario,
    const ParameterFile &file,
    const ParameterCombination &combination,
    const std::string &scenario_id)
  {
    using namespace nmopt::application;
    require_parameter(file, "Benchmark/id", "chapter-6.b2.graetz-flow");
    require_parameter(file, "Benchmark/recipe", chapter6::b2_recipe_id);
    scenario.problem.transport_boundary_form =
      b2_transport_boundary_form(file);
    scenario.problem.recipe.control_discretisation =
      b2_neumann_control_discretisation(file);
    scenario.compile.volume_observation =
      b2_volume_observation_options(file);
    require_parameter(file, "Boundary/normal orientation", "outward");
    require_parameter(file, "Boundary/trace evaluation", "fe-q-state-trace");
    require_parameter(file, "Boundary/face quadrature", "qgauss-face");
    require_parameter(file, "Boundary/fixed region", "dirichlet-boundary");
    require_parameter(file, "Boundary/control region", "control-boundary");
    require_parameter(file, "Boundary/outflow region", "outflow-boundary");
    require_parameter(file, "Boundary/fixed boundary id", "0");
    require_parameter(file, "Boundary/control boundary id", "1");
    require_parameter(file, "Boundary/outflow boundary id", "2");
    require_parameter(file, "Boundary/upstream transition x", "1.0");
    require_parameter(file, "Boundary/outflow x", "4.0");
    require_parameter(file, "Compile/stabilization", "galerkin");
    const auto method_id = file.value("Solver/method");
    const auto method = parse_method(method_id);
    if (method != chapter6::ReducedMethod::bfgs)
      throw std::invalid_argument("B2 selects the BFGS method");
    scenario.solver.method = method;
    require_parameter(
      file,
      "Output/selected fields",
      "state, state-uncontrolled, control, adjoint, negative-adjoint, target, forcing, observation-region");

    scenario.problem.recipe.observed_material_id =
      parameter_unsigned(file, "Observation/material id");
    scenario.problem.recipe.with_facewise_box =
      parameter_bool(file, "Problem/facewise box constraint");
    scenario.problem.data.diffusion = parameter_double(file, "Runtime/diffusion");
    scenario.problem.data.reaction = parameter_double(file, "Runtime/reaction");
    scenario.problem.data.regularisation_weight =
      combination.values.count("regularisation")
        ? parse_number_text(combination_value(combination, "regularisation"),
                            "Matrix/regularisation")
        : parameter_double(file, "Runtime/regularisation");
    scenario.problem.fixed_dirichlet_data =
      parameter_scalar_function_definition(file,
                                           "Functions/fixed Dirichlet data");
    scenario.problem.data.fixed_dirichlet_data_provenance =
      scenario.problem.fixed_dirichlet_data.provenance;
    scenario.problem.conservative_transport =
      parameter_rank_one_vector_function_definition(
        file, "Functions/conservative transport");
    scenario.problem.data.conservative_transport_provenance =
      scenario.problem.conservative_transport.provenance;
    scenario.problem.observation_region =
      b2_observation_region_from_combination(combination);
    const auto target_axis = combination.values.find("target-profile");
    if (target_axis == combination.values.end())
      throw std::invalid_argument(
        "B2 combinations need a target-profile axis");
    scenario.problem.target_profile = target_axis->second;
    scenario.problem.target_catalog =
      b2_target_catalog(file, scenario.problem.target_profile);

    const auto forcing_axis = combination.values.find("forcing");
    if (forcing_axis == combination.values.end())
      scenario.problem.forcing = parameter_scalar_function_definition(
        file, "Functions/forcing");
    else
      scenario.problem.forcing =
        parameter_scalar_function_definition_from_catalog(
          file, "Functions/forcing definitions/", forcing_axis->second);
    scenario.problem.data.forcing_provenance =
      scenario.problem.forcing.provenance;
    scenario.problem.data.desired_state_provenance =
      chapter6::b2_target_definition(scenario.problem.target_catalog).provenance;

    apply_common_parameter_options(scenario, file, scenario_id);
    apply_solver_options(scenario.solver, file, method_id);
  }
} // namespace nmopt::application::runner::binding
