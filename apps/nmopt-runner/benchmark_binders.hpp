#pragma once

#include "nmopt/application/chapter6.hpp"
#include "parameter_binding.hpp"

#include <stdexcept>
#include <string>

namespace nmopt::application::runner::binding
{
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
    require_parameter(file, "Functions/desired state", "b1-polynomial");
    require_parameter(file, "Functions/desired state/kind", "polynomial");
    require_parameter(file,
                      "Functions/desired state/expression",
                      "10*x0*(1-x0)*x1*(1-x1)");
    require_parameter(file,
                      "Mesh/geometry",
                      "unit-hypercube");
    require_parameter(file, "Solver/method", "from-matrix");
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
    scenario.problem.data.desired_state_provenance =
      file.value("Functions/desired state/provenance");
    scenario.problem.regularisation_sweep = {beta};
    scenario.problem.forcing =
      parameter_scalar_function_definition(file, "Functions/forcing");
    scenario.problem.data.forcing_provenance =
      scenario.problem.forcing.provenance;
    scenario.solver.method = parse_method(method_id);
    apply_common_parameter_options(scenario, file, scenario_id);
    apply_solver_options(scenario.solver, file, method_id);
    if (scenario.solver.method == chapter6::ReducedMethod::bfgs)
      throw std::invalid_argument("B1 does not support full BFGS");
  }
} // namespace nmopt::application::runner::binding
