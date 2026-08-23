#pragma once

#include "nmopt/compiler/v1/dealii_reference_cell.hpp"

#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <set>
#include <utility>
#include <vector>

namespace nmopt::compiler::v1::detail
{
  // Discrete material-subdomain realization of
  //
  //   1/2 ||y - y_d||^2_{L2(omega_o)}.
  //
  // The stored load and constant include the exact affine correction induced
  // by fixed physical state values. The free-coordinate objective therefore
  // remains 1/2 y^T M y - q^T y + 1/2 c.
  template <int dim>
  class VolumeObservationAssembly final
  {
  public:
    using Vector = dealii::Vector<double>;

    VolumeObservationAssembly(
      const dealii::DoFHandler<dim> &             state_dof_handler,
      const dealii::FiniteElement<dim> &           state_fe,
      const std::vector<bool> &                    constrained_state_dofs,
      const Vector &                               fixed_state_values,
      std::set<dealii::types::material_id>         observation_material_ids,
      const dealii::Function<dim> &                desired_state,
      const unsigned int                           quadrature_order)
      : observation_material_ids_(std::move(observation_material_ids))
    {
      contract::require(!observation_material_ids_.empty(),
                        "Volume observation needs a material selection");
      contract::require(quadrature_order > 0,
                        "Volume observation needs a positive quadrature order");
      contract::require(
        constrained_state_dofs.size() == state_dof_handler.n_dofs() &&
          fixed_state_values.size() == state_dof_handler.n_dofs(),
        "Volume observation received incompatible fixed-state coordinates");

      dealii::DynamicSparsityPattern state_dsp(state_dof_handler.n_dofs(),
                                                state_dof_handler.n_dofs());
      dealii::DoFTools::make_sparsity_pattern(state_dof_handler, state_dsp);
      state_sparsity_.copy_from(state_dsp);
      state_tracking_matrix_.reinit(state_sparsity_);
      desired_state_load_.reinit(state_dof_handler.n_dofs());
      assemble(state_dof_handler,
               state_fe,
               constrained_state_dofs,
               fixed_state_values,
               desired_state,
               quadrature_order);
    }

    const dealii::SparseMatrix<double> &
    state_tracking_matrix() const
    {
      return state_tracking_matrix_;
    }

    const Vector &
    desired_state_load() const
    {
      return desired_state_load_;
    }

    double
    desired_state_norm() const
    {
      return desired_state_norm_;
    }

  private:
    void
    assemble(const dealii::DoFHandler<dim> &state_dof_handler,
             const dealii::FiniteElement<dim> &state_fe,
             const std::vector<bool> &constrained_state_dofs,
             const Vector &fixed_state_values,
             const dealii::Function<dim> &desired_state,
             const unsigned int quadrature_order)
    {
      const auto volume_quadrature = make_gauss_volume_quadrature(
        state_dof_handler.get_triangulation(),
        quadrature_order,
        "The volume-observation realization");
      dealii::FEValues<dim> state_values(
        state_fe,
        *volume_quadrature,
        dealii::update_values | dealii::update_quadrature_points |
          dealii::update_JxW_values);
      dealii::FullMatrix<double> local_state_tracking(state_fe.dofs_per_cell,
                                                       state_fe.dofs_per_cell);
      dealii::Vector<double> local_desired_state(state_fe.dofs_per_cell);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe.dofs_per_cell);

      for (auto cell = state_dof_handler.begin_active();
           cell != state_dof_handler.end();
           ++cell)
        {
          if (observation_material_ids_.count(cell->material_id()) == 0)
            continue;
          state_values.reinit(cell);
          local_state_tracking = 0.0;
          local_desired_state = 0.0;
          double fixed_tracking_value = 0.0;
          double fixed_desired_value = 0.0;
          for (unsigned int q = 0; q < volume_quadrature->size(); ++q)
            {
              const double weight = state_values.JxW(q);
              const double desired_value =
                desired_state.value(state_values.quadrature_point(q));
              desired_state_norm_ += desired_value * desired_value * weight;
              for (unsigned int i = 0; i < state_fe.dofs_per_cell; ++i)
                {
                  const double phi_i = state_values.shape_value(i, q);
                  local_desired_state(i) += desired_value * phi_i * weight;
                  for (unsigned int j = 0; j < state_fe.dofs_per_cell; ++j)
                    local_state_tracking(i, j) +=
                      phi_i * state_values.shape_value(j, q) * weight;
                }
            }

          cell->get_dof_indices(state_indices);
          for (unsigned int i = 0; i < state_fe.dofs_per_cell; ++i)
            {
              const auto global_i = state_indices[i];
              if (constrained_state_dofs.at(global_i))
                {
                  fixed_desired_value +=
                    local_desired_state(i) * fixed_state_values[global_i];
                  for (unsigned int j = 0; j < state_fe.dofs_per_cell; ++j)
                    {
                      const auto global_j = state_indices[j];
                      if (constrained_state_dofs.at(global_j))
                        fixed_tracking_value +=
                          fixed_state_values[global_i] *
                          local_state_tracking(i, j) *
                          fixed_state_values[global_j];
                    }
                  continue;
                }

              desired_state_load_[global_i] += local_desired_state(i);
              for (unsigned int j = 0; j < state_fe.dofs_per_cell; ++j)
                {
                  const auto global_j = state_indices[j];
                  if (constrained_state_dofs.at(global_j))
                    desired_state_load_[global_i] -=
                      local_state_tracking(i, j) * fixed_state_values[global_j];
                  else
                    state_tracking_matrix_.add(global_i,
                                               global_j,
                                               local_state_tracking(i, j));
                }
            }
          desired_state_norm_ +=
            fixed_tracking_value - 2.0 * fixed_desired_value;
        }
    }

    const std::set<dealii::types::material_id> observation_material_ids_;
    dealii::SparsityPattern state_sparsity_;
    dealii::SparseMatrix<double> state_tracking_matrix_;
    Vector desired_state_load_;
    double desired_state_norm_ = 0.0;
  };
} // namespace nmopt::compiler::v1::detail
