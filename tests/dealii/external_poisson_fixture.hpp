#pragma once

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <locale>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace external_poisson_application
{
  class PoissonControlApplication final
  {
  public:
    using Vector = dealii::Vector<double>;
    using Matrix = dealii::SparseMatrix<double>;

    struct ObjectiveDerivative
    {
      Vector state;
      Vector control;
    };

    PoissonControlApplication(const unsigned int refinement = 3,
                              const double       forcing = 1.0,
                              const double       desired_state = 0.25,
                              const double       regularisation = 0.2)
      : fe_(1)
      , dof_handler_(triangulation_)
      , forcing_(forcing)
      , desired_state_(desired_state)
      , regularisation_(regularisation)
    {
      if (!(regularisation_ > 0.0))
        throw std::invalid_argument(
          "external Poisson regularisation must be positive");

      dealii::GridGenerator::hyper_cube(triangulation_, 0.0, 1.0);
      triangulation_.refine_global(refinement);
      dof_handler_.distribute_dofs(fe_);
      build_constraints();
      initialise_storage();
      assemble();
    }

    std::size_t
    state_dimension() const
    {
      return static_cast<std::size_t>(dof_handler_.n_dofs());
    }

    std::size_t
    control_dimension() const
    {
      return state_dimension();
    }

    std::shared_ptr<const Matrix>
    control_mass_matrix() const
    {
      return control_mass_;
    }

    Vector
    residual(const Vector &state, const Vector &control) const
    {
      require_vector_size(state, "state");
      require_vector_size(control, "control");

      Vector value(state_dimension());
      system_matrix_.vmult(value, state);
      value.add(-1.0, forcing_load_);

      Vector control_contribution(state_dimension());
      control_coupling_.vmult(control_contribution, control);
      value.add(-1.0, control_contribution);
      return value;
    }

    Vector
    residual_jvp(const Vector &state_tangent,
                 const Vector &control_tangent) const
    {
      require_vector_size(state_tangent, "state tangent");
      require_vector_size(control_tangent, "control tangent");

      Vector value(state_dimension());
      system_matrix_.vmult(value, state_tangent);

      Vector control_contribution(state_dimension());
      control_coupling_.vmult(control_contribution, control_tangent);
      value.add(-1.0, control_contribution);
      return value;
    }

    Vector
    residual_vjp(const Vector &test_seed) const
    {
      require_vector_size(test_seed, "test seed");

      Vector state(state_dimension());
      system_matrix_.Tvmult(state, test_seed);
      Vector control(state_dimension());
      control_coupling_.Tvmult(control, test_seed);
      control *= -1.0;

      Vector result(2 * state_dimension());
      for (std::size_t index = 0; index < state_dimension(); ++index)
        {
          result[index] = state[index];
          result[state_dimension() + index] = control[index];
        }
      return result;
    }

    double
    objective(const Vector &state, const Vector &control) const
    {
      require_vector_size(state, "state");
      require_vector_size(control, "control");

      Vector state_mass_times_state(state_dimension());
      state_mass_.vmult(state_mass_times_state, state);
      const double state_value =
        0.5 * (state * state_mass_times_state) -
        (desired_state_load_ * state) + 0.5 * desired_state_norm_;

      Vector control_mass_times_control(state_dimension());
      control_mass_->vmult(control_mass_times_control, control);
      const double control_value =
        0.5 * regularisation_ * (control * control_mass_times_control);
      return state_value + control_value;
    }

    ObjectiveDerivative
    objective_derivative(const Vector &state, const Vector &control) const
    {
      require_vector_size(state, "state");
      require_vector_size(control, "control");

      ObjectiveDerivative result{Vector(state_dimension()),
                                 Vector(state_dimension())};
      state_mass_.vmult(result.state, state);
      result.state.add(-1.0, desired_state_load_);
      control_mass_->vmult(result.control, control);
      result.control *= regularisation_;
      return result;
    }

    Vector
    solve_state(const Vector &control) const
    {
      require_vector_size(control, "control");

      Vector right_hand_side = forcing_load_;
      Vector control_contribution(state_dimension());
      control_coupling_.vmult(control_contribution, control);
      right_hand_side.add(1.0, control_contribution);

      Vector state(state_dimension());
      dealii::SparseDirectUMFPACK solver;
      solver.initialize(system_matrix_);
      solver.vmult(state, right_hand_side);
      return state;
    }

    Vector
    solve_adjoint(const Vector &state_objective_derivative) const
    {
      require_vector_size(state_objective_derivative,
                          "state objective derivative");

      Vector adjoint(state_dimension());
      dealii::SparseDirectUMFPACK solver;
      solver.initialize(system_matrix_);
      solver.vmult(adjoint, state_objective_derivative);
      return adjoint;
    }

    double
    residual_norm(const Vector &state, const Vector &control) const
    {
      return residual(state, control).l2_norm();
    }

    void
    write_native_output(const std::filesystem::path &directory,
                        const Vector &                 state) const
    {
      require_vector_size(state, "output state");
      std::filesystem::create_directories(directory);

      dealii::DataOut<1> data_out;
      data_out.attach_dof_handler(dof_handler_);
      data_out.add_data_vector(state, "state");
      data_out.build_patches();

      std::ofstream output(directory / "fields-volume.vtu");
      if (!output)
        throw std::runtime_error(
          "could not open external Poisson native output");
      output.imbue(std::locale::classic());
      data_out.write_vtu(output);
      if (!output)
        throw std::runtime_error(
          "could not write external Poisson native output");
    }

  private:
    void
    require_vector_size(const Vector &vector, const char *name) const
    {
      if (static_cast<std::size_t>(vector.size()) != state_dimension())
        throw std::invalid_argument(std::string("external Poisson ") + name +
                                    " has the wrong dimension");
    }

    void
    build_constraints()
    {
      constraints_.clear();
      dealii::DoFTools::make_hanging_node_constraints(dof_handler_,
                                                       constraints_);
      dealii::Functions::ZeroFunction<1> zero;
      dealii::VectorTools::interpolate_boundary_values(
        dof_handler_, 0, zero, constraints_);
      dealii::VectorTools::interpolate_boundary_values(
        dof_handler_, 1, zero, constraints_);
      constraints_.close();

      constrained_.assign(dof_handler_.n_dofs(), false);
      for (dealii::types::global_dof_index index = 0;
           index < dof_handler_.n_dofs();
           ++index)
        constrained_[index] = constraints_.is_constrained(index);
    }

    void
    initialise_storage()
    {
      const auto size = dof_handler_.n_dofs();
      dealii::DynamicSparsityPattern sparsity(size, size);
      dealii::DoFTools::make_sparsity_pattern(dof_handler_, sparsity);

      state_sparsity_.copy_from(sparsity);
      system_matrix_.reinit(state_sparsity_);
      state_mass_.reinit(state_sparsity_);
      control_coupling_.reinit(state_sparsity_);

      control_mass_ = std::make_shared<Matrix>();
      control_mass_->reinit(state_sparsity_);
      forcing_load_.reinit(size);
      desired_state_load_.reinit(size);
    }

    void
    assemble()
    {
      const dealii::QGauss<1> quadrature(fe_.degree + 1);
      dealii::FEValues<1> values(
        fe_,
        quadrature,
        dealii::update_values | dealii::update_gradients |
          dealii::update_JxW_values);
      dealii::FullMatrix<double> local_system(fe_.dofs_per_cell,
                                              fe_.dofs_per_cell);
      dealii::FullMatrix<double> local_mass(fe_.dofs_per_cell,
                                            fe_.dofs_per_cell);
      dealii::Vector<double> local_forcing(fe_.dofs_per_cell);
      dealii::Vector<double> local_desired_state(fe_.dofs_per_cell);
      std::vector<dealii::types::global_dof_index> indices(fe_.dofs_per_cell);

      for (const auto &cell : dof_handler_.active_cell_iterators())
        {
          values.reinit(cell);
          local_system = 0.0;
          local_mass = 0.0;
          local_forcing = 0.0;
          local_desired_state = 0.0;

          for (unsigned int q = 0; q < quadrature.size(); ++q)
            for (unsigned int i = 0; i < fe_.dofs_per_cell; ++i)
              {
                const double phi_i = values.shape_value(i, q);
                local_forcing(i) += forcing_ * phi_i * values.JxW(q);
                local_desired_state(i) +=
                  desired_state_ * phi_i * values.JxW(q);
                for (unsigned int j = 0; j < fe_.dofs_per_cell; ++j)
                  {
                    const double phi_j = values.shape_value(j, q);
                    local_system(i, j) +=
                      (values.shape_grad(i, q) * values.shape_grad(j, q) +
                       0.25 * phi_i * phi_j) *
                      values.JxW(q);
                    local_mass(i, j) += phi_i * phi_j * values.JxW(q);
                  }
              }

          cell->get_dof_indices(indices);
          for (unsigned int i = 0; i < fe_.dofs_per_cell; ++i)
            {
              const auto global_i = indices[i];
              if (!constrained_[global_i])
                {
                  forcing_load_[global_i] += local_forcing(i);
                  desired_state_load_[global_i] += local_desired_state(i);
                  for (unsigned int j = 0; j < fe_.dofs_per_cell; ++j)
                    {
                      const auto global_j = indices[j];
                      control_coupling_.add(global_i,
                                            global_j,
                                            local_mass(i, j));
                      if (!constrained_[global_j])
                        {
                          system_matrix_.add(global_i,
                                             global_j,
                                             local_system(i, j));
                          state_mass_.add(global_i,
                                           global_j,
                                           local_mass(i, j));
                        }
                    }
                }

              for (unsigned int j = 0; j < fe_.dofs_per_cell; ++j)
                control_mass_->add(global_i,
                                   indices[j],
                                   local_mass(i, j));
            }
      }

      desired_state_norm_ = desired_state_ * desired_state_;

      for (dealii::types::global_dof_index index = 0;
           index < dof_handler_.n_dofs();
           ++index)
        if (constrained_[index])
          system_matrix_.set(index, index, 1.0);
    }

    dealii::Triangulation<1>              triangulation_;
    dealii::FE_Q<1>                       fe_;
    dealii::DoFHandler<1>                 dof_handler_;
    dealii::AffineConstraints<double>     constraints_;
    std::vector<bool>                     constrained_;
    dealii::SparsityPattern              state_sparsity_;
    Matrix                                system_matrix_;
    Matrix                                state_mass_;
    Matrix                                control_coupling_;
    std::shared_ptr<Matrix>               control_mass_;
    Vector                                forcing_load_;
    Vector                                desired_state_load_;
    double                                desired_state_norm_ = 0.0;
    double                                forcing_;
    double                                desired_state_;
    double                                regularisation_;
  };
} // namespace external_poisson_application
