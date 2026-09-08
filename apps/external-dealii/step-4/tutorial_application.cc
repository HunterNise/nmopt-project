#include "tutorial_application.hpp"

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#define NMOPT_EXTERNAL_TUTORIAL_STEP_4_NO_MAIN
#include "step-4.cc"
#undef NMOPT_EXTERNAL_TUTORIAL_STEP_4_NO_MAIN

namespace
{
  using Index = dealii::types::global_dof_index;

  void
  require_size(const TutorialApplication::Vector &vector,
               const std::size_t                  expected,
               const char *const                  name)
  {
    if (static_cast<std::size_t>(vector.size()) != expected)
      throw std::invalid_argument(std::string("tutorial ") + name +
                                  " has the wrong dimension");
  }
}

struct TutorialApplication::Impl
{
  using Matrix = TutorialApplication::Matrix;
  using Vector = TutorialApplication::Vector;

  Step4<2> tutorial;
  dealii::SparsityPattern metric_sparsity;
  std::shared_ptr<Matrix> metric;
  double regularisation = 0.2;

  Impl()
  {
    tutorial.prepare_for_external_use();
    assemble_metric();
  }

  void
  assemble_metric()
  {
    const auto &dof_handler = tutorial.dof_handler_view();
    const auto  size = dof_handler.n_dofs();

    dealii::DynamicSparsityPattern dsp(size, size);
    dealii::DoFTools::make_sparsity_pattern(dof_handler, dsp);
    metric_sparsity.copy_from(dsp);

    metric = std::make_shared<Matrix>();
    metric->reinit(metric_sparsity);

    const auto &finite_element = tutorial.finite_element_view();
    const dealii::QGauss<2> quadrature(finite_element.degree + 1);
    dealii::FEValues<2> values(
      finite_element, quadrature, dealii::update_values |
                                    dealii::update_JxW_values);
    dealii::FullMatrix<double> local_mass(finite_element.dofs_per_cell,
                                          finite_element.dofs_per_cell);
    std::vector<Index> indices(finite_element.dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        values.reinit(cell);
        local_mass = 0.0;
        for (const auto q : values.quadrature_point_indices())
          for (const auto i : values.dof_indices())
            for (const auto j : values.dof_indices())
              local_mass(i, j) += values.shape_value(i, q) *
                                  values.shape_value(j, q) * values.JxW(q);

        cell->get_dof_indices(indices);
        for (const auto i : values.dof_indices())
          for (const auto j : values.dof_indices())
            metric->add(indices[i], indices[j], local_mass(i, j));
      }
  }

  std::size_t
  dimension() const
  {
    return static_cast<std::size_t>(tutorial.dof_handler_view().n_dofs());
  }
};

TutorialApplication::TutorialApplication()
  : impl_(std::make_unique<Impl>())
{}

TutorialApplication::~TutorialApplication() = default;

TutorialApplication::TutorialApplication(TutorialApplication &&) noexcept =
  default;

TutorialApplication &
TutorialApplication::operator=(TutorialApplication &&) noexcept = default;

std::size_t
TutorialApplication::state_dimension() const
{
  return impl_->dimension();
}

std::size_t
TutorialApplication::control_dimension() const
{
  return impl_->dimension();
}

TutorialApplication::Vector
TutorialApplication::residual(const Vector &state, const Vector &control) const
{
  require_size(state, state_dimension(), "state");
  require_size(control, control_dimension(), "control");

  Vector value(state_dimension());
  impl_->tutorial.system_matrix_view().vmult(value, state);
  value.add(-1.0, impl_->tutorial.system_rhs_view());
  value.add(-1.0, control);
  return value;
}

TutorialApplication::Vector
TutorialApplication::residual_jvp(const Vector &state_tangent,
                                  const Vector &control_tangent) const
{
  require_size(state_tangent, state_dimension(), "state tangent");
  require_size(control_tangent, control_dimension(), "control tangent");

  Vector value(state_dimension());
  impl_->tutorial.system_matrix_view().vmult(value, state_tangent);
  value.add(-1.0, control_tangent);
  return value;
}

TutorialApplication::ResidualDerivative
TutorialApplication::residual_vjp(const Vector &test_seed) const
{
  require_size(test_seed, state_dimension(), "test seed");

  ResidualDerivative result{Vector(state_dimension()),
                            Vector(control_dimension())};
  impl_->tutorial.system_matrix_view().Tvmult(result.state, test_seed);
  result.control = test_seed;
  result.control *= -1.0;
  return result;
}

double
TutorialApplication::objective(const Vector &state, const Vector &control) const
{
  require_size(state, state_dimension(), "state");
  require_size(control, control_dimension(), "control");

  Vector state_metric(state_dimension());
  Vector control_metric(control_dimension());
  impl_->metric->vmult(state_metric, state);
  impl_->metric->vmult(control_metric, control);
  return 0.5 * (state * state_metric) +
         0.5 * impl_->regularisation * (control * control_metric);
}

TutorialApplication::ObjectiveDerivative
TutorialApplication::objective_derivative(const Vector &state,
                                           const Vector &control) const
{
  require_size(state, state_dimension(), "state");
  require_size(control, control_dimension(), "control");

  ObjectiveDerivative result{Vector(state_dimension()),
                             Vector(control_dimension())};
  impl_->metric->vmult(result.state, state);
  impl_->metric->vmult(result.control, control);
  result.control *= impl_->regularisation;
  return result;
}

TutorialApplication::Vector
TutorialApplication::solve_state(const Vector &control) const
{
  require_size(control, control_dimension(), "control");

  Vector right_hand_side = impl_->tutorial.system_rhs_view();
  right_hand_side.add(1.0, control);

  Vector state(state_dimension());
  dealii::SparseDirectUMFPACK solver;
  solver.initialize(impl_->tutorial.system_matrix_view());
  solver.vmult(state, right_hand_side);
  return state;
}

TutorialApplication::Vector
TutorialApplication::solve_adjoint(
  const Vector &state_objective_derivative) const
{
  require_size(state_objective_derivative,
               state_dimension(),
               "state objective derivative");

  Vector adjoint(state_dimension());
  dealii::SparseDirectUMFPACK solver;
  solver.initialize(impl_->tutorial.system_matrix_view());
  solver.vmult(adjoint, state_objective_derivative);
  return adjoint;
}

std::shared_ptr<const TutorialApplication::Matrix>
TutorialApplication::control_metric_matrix() const
{
  return impl_->metric;
}

void
TutorialApplication::write_native_output(
  const std::filesystem::path &directory, const Vector &state) const
{
  require_size(state, state_dimension(), "output state");
  std::filesystem::create_directories(directory);

  dealii::DataOut<2> data_out;
  data_out.attach_dof_handler(impl_->tutorial.dof_handler_view());
  data_out.add_data_vector(state, "state");
  data_out.build_patches();

  std::ofstream output(directory / "fields-volume.vtu");
  if (!output)
    throw std::runtime_error("could not open tutorial native output");
  data_out.write_vtu(output);
  if (!output)
    throw std::runtime_error("could not write tutorial native output");
}
