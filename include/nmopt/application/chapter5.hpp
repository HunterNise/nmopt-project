#pragma once

#include "nmopt/application/catalog.hpp"
#include "nmopt/application/recipe.hpp"
#include "nmopt/semantic/v1/reference_specs.hpp"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::application::chapter5
{
  inline constexpr const char *scalar_distributed_recipe_id =
    "chapter-5.scalar-diffusion-reaction-volume";
  inline constexpr const char *general_scalar_recipe_id =
    "chapter-5.scalar-elliptic-robin-volume";
  inline constexpr const char *subdomain_recipe_id =
    "chapter-5.scalar-subdomain-tracking";
  inline constexpr const char *neumann_recipe_id =
    "chapter-5.scalar-neumann-boundary";
  inline constexpr const char *neumann_convection_recipe_id =
    "chapter-5.scalar-neumann-convection-subdomain";

  // These records describe semantic choices only. Runtime Function objects,
  // meshes, compiler policies, and solver policies belong to the separate
  // scenario and backend option records.
  enum class DistributedControlDiscretisation
  {
    cellwise_constant,
    homogeneous_dirichlet_continuous
  };

  inline const char *
  distributed_control_discretisation_name(
    const DistributedControlDiscretisation discretisation)
  {
    switch (discretisation)
      {
        case DistributedControlDiscretisation::cellwise_constant:
          return "cellwise-volume";
        case DistributedControlDiscretisation::homogeneous_dirichlet_continuous:
          return "continuous-volume-homogeneous-dirichlet";
      }
    throw std::invalid_argument(
      "unknown distributed-control discretisation");
  }

  struct ScalarDistributedControlParameters
  {
    bool with_cellwise_box = false;
    DistributedControlDiscretisation discretisation =
      DistributedControlDiscretisation::cellwise_constant;
  };

  struct GeneralScalarParameters
  {
    std::vector<unsigned int> fixed_dirichlet_boundary_ids = {0};
    std::vector<unsigned int> robin_boundary_ids = {1};
    bool                      with_cellwise_box = false;
  };

  struct SubdomainTrackingParameters
  {
    unsigned int observed_material_id = 1;
    bool         with_cellwise_box = false;
  };

  struct NeumannBoundaryParameters
  {
    bool with_facewise_box = false;
  };

  struct NeumannConvectionParameters
  {
    unsigned int observed_material_id = 1;
    bool         with_facewise_box = false;
  };

  using ScalarDistributedRecipe =
    ProblemRecipeT<ScalarDistributedControlParameters>;
  using GeneralScalarRecipe = ProblemRecipeT<GeneralScalarParameters>;
  using SubdomainTrackingRecipe = ProblemRecipeT<SubdomainTrackingParameters>;
  using NeumannBoundaryRecipe = ProblemRecipeT<NeumannBoundaryParameters>;
  using NeumannConvectionRecipe =
    ProblemRecipeT<NeumannConvectionParameters>;

  inline ScalarDistributedRecipe
  make_scalar_distributed_recipe()
  {
    return ScalarDistributedRecipe{
      {scalar_distributed_recipe_id,
       "Scalar diffusion-reaction volume control",
       "Chapter 5 baseline distributed-control recipe",
       "chapter-5",
       {"scalar", "volume-control", "l2-cellwise", "l2-continuous"}},
      [](const ScalarDistributedControlParameters &parameters) {
        switch (parameters.discretisation)
          {
            case DistributedControlDiscretisation::cellwise_constant:
              return semantic::v1::make_scalar_diffusion_reaction_problem(
                parameters.with_cellwise_box);
            case DistributedControlDiscretisation::
              homogeneous_dirichlet_continuous:
              if (parameters.with_cellwise_box)
                throw std::invalid_argument(
                  "continuous distributed control does not support the cellwise box");
              return semantic::v1::
                make_l2_state_tracking_continuous_control_problem();
          }
        throw std::invalid_argument(
          "distributed-control recipe has an unknown discretisation");
      }};
  }

  inline GeneralScalarRecipe
  make_general_scalar_recipe()
  {
    return GeneralScalarRecipe{
      {general_scalar_recipe_id,
       "General scalar elliptic Robin volume control",
       "Chapter 5 scalar residual-term and Robin composition",
       "chapter-5",
       {"tensor-diffusion", "transport", "robin", "scalar"}},
      [](const GeneralScalarParameters &parameters) {
        return semantic::v1::make_general_scalar_elliptic_robin_problem(
          parameters.fixed_dirichlet_boundary_ids,
          parameters.robin_boundary_ids,
          parameters.with_cellwise_box);
      }};
  }

  inline SubdomainTrackingRecipe
  make_subdomain_tracking_recipe()
  {
    return SubdomainTrackingRecipe{
      {subdomain_recipe_id,
       "Scalar material-subdomain tracking",
       "Chapter 5 volume observation on a declared material region",
       "chapter-5",
       {"scalar", "material-observation", "l2-cellwise"}},
      [](const SubdomainTrackingParameters &parameters) {
        return semantic::v1::make_subdomain_tracking_scalar_diffusion_reaction_problem(
          parameters.observed_material_id,
          parameters.with_cellwise_box);
      }};
  }

  inline NeumannBoundaryRecipe
  make_neumann_boundary_recipe()
  {
    return NeumannBoundaryRecipe{
      {neumann_recipe_id,
       "Scalar Neumann boundary control",
       "Chapter 5 facewise boundary-control recipe",
       "chapter-5",
       {"scalar", "neumann-control", "l2-facewise"}},
      [](const NeumannBoundaryParameters &parameters) {
        return semantic::v1::make_neumann_boundary_control_problem(
          parameters.with_facewise_box);
      }};
  }

  inline NeumannConvectionRecipe
  make_neumann_convection_recipe()
  {
    return NeumannConvectionRecipe{
      {neumann_convection_recipe_id,
       "Scalar conservative-transport Neumann control",
       "Chapter 5 C5.6 boundary control with material-subdomain tracking",
       "chapter-5",
       {"scalar", "neumann-control", "conservative-transport",
        "material-observation"}},
      [](const NeumannConvectionParameters &parameters) {
        return semantic::v1::make_neumann_convection_subdomain_tracking_problem(
          parameters.observed_material_id,
          parameters.with_facewise_box);
      }};
  }

  inline ApplicationCatalog
  make_catalog()
  {
    ApplicationCatalog catalog;
    catalog.add(make_scalar_distributed_recipe().metadata());
    catalog.add(make_general_scalar_recipe().metadata());
    catalog.add(make_subdomain_tracking_recipe().metadata());
    catalog.add(make_neumann_boundary_recipe().metadata());
    catalog.add(make_neumann_convection_recipe().metadata());
    return catalog;
  }
} // namespace nmopt::application::chapter5
