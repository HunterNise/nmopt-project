#pragma once

#include "nmopt/contract/linalg.hpp"

#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/grid/tria.h>

#include <memory>
#include <string>

namespace nmopt::compiler::v1::detail
{
  template <int dim>
  std::unique_ptr<dealii::FiniteElement<dim>>
  make_scalar_lagrange_element(
    const dealii::Triangulation<dim> &triangulation,
    const unsigned int                degree,
    const char *                      target_name)
  {
    if (triangulation.all_reference_cells_are_hyper_cube())
      return std::make_unique<dealii::FE_Q<dim>>(degree);
    contract::require(
      triangulation.all_reference_cells_are_simplex(),
      std::string(target_name) +
        " supports only one simplex or hypercube reference-cell family");
    return std::make_unique<dealii::FE_SimplexP<dim>>(degree);
  }

  template <int dim>
  std::unique_ptr<dealii::Quadrature<dim>>
  make_gauss_volume_quadrature(
    const dealii::Triangulation<dim> &triangulation,
    const unsigned int                order,
    const char *                      target_name)
  {
    if (triangulation.all_reference_cells_are_hyper_cube())
      return std::make_unique<dealii::QGauss<dim>>(order);
    contract::require(
      triangulation.all_reference_cells_are_simplex(),
      std::string(target_name) +
        " supports only one simplex or hypercube reference-cell family");
    return std::make_unique<dealii::QGaussSimplex<dim>>(order);
  }

  template <int dim>
  std::unique_ptr<dealii::Quadrature<dim - 1>>
  make_gauss_face_quadrature(
    const dealii::Triangulation<dim> &triangulation,
    const unsigned int                order,
    const char *                      target_name)
  {
    if (triangulation.all_reference_cells_are_hyper_cube())
      return std::make_unique<dealii::QGauss<dim - 1>>(order);
    contract::require(
      triangulation.all_reference_cells_are_simplex(),
      std::string(target_name) +
        " supports only one simplex or hypercube reference-cell family");
    return std::make_unique<dealii::QGaussSimplex<dim - 1>>(order);
  }
} // namespace nmopt::compiler::v1::detail
