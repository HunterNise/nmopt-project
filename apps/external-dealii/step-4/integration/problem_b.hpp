#pragma once

#include "problem_b_mass.hpp"
#include "../diagnostics/instrumentation.hpp"

#include <deal.II/lac/vector.h>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

namespace external_dealii_step4
{
  template <int dim, typename Application>
  class ProblemB final
  {
  public:
    using Coordinates   = ProblemBCoordinates<dim>;
    using Mass          = ProblemBMass<dim>;
    using Vector        = dealii::Vector<double>;
    using SolveEvidence = decltype(std::declval<const Application &>().solve(
      std::declval<const Vector &>(), std::declval<Vector &>()));

    struct ResidualDerivative
    {
      Vector state;
      Vector control;
    };

    struct ObjectiveDerivative
    {
      Vector state;
      Vector control;
    };

    struct SolveResult
    {
      Vector         solution;
      Vector         full_solution;
      SolveEvidence evidence;
    };

    explicit ProblemB(Application &application,
                      Instrumentation *instrumentation = nullptr)
      : application_(application)
      , instrumentation_(instrumentation)
      , coordinates_(application.dof_handler_view(),
                     application.boundary_values_view())
      , mass_(application.dof_handler_view(), coordinates_, instrumentation_)
      , free_system_rhs_(coordinates_.restrict(application.system_rhs_view()))
    {
      if (application.system_matrix_view().m() != coordinates_.full_dimension() ||
          application.system_matrix_view().n() != coordinates_.full_dimension() ||
          application.system_rhs_view().size() != coordinates_.full_dimension())
        throw std::invalid_argument(
          "Problem B received incompatible Step-4 system dimensions");
    }

    ProblemB(const ProblemB &) = delete;
    ProblemB &operator=(const ProblemB &) = delete;
    ProblemB(ProblemB &&) = delete;
    ProblemB &operator=(ProblemB &&) = delete;

    std::size_t
    state_dimension() const
    {
      return coordinates_.free_dimension();
    }

    std::size_t
    control_dimension() const
    {
      return coordinates_.full_dimension();
    }

    double
    alpha() const
    {
      return 1.0;
    }

    const Coordinates &
    coordinates() const
    {
      return coordinates_;
    }

    const Mass &
    mass() const
    {
      return mass_;
    }

    const Vector &
    free_system_rhs() const
    {
      return free_system_rhs_;
    }

    Vector
    residual(const Vector &state, const Vector &control) const
    {
      require_size(state, state_dimension(), "state");
      require_size(control, control_dimension(), "control");

      Vector value =
        apply_state_operator(state, ProblemBMatrixPurpose::residual);
      value.add(-1.0, free_system_rhs_);
      value.add(-1.0,
                mass_.coupling_apply(control, ProblemBMatrixPurpose::residual));
      return value;
    }

    Vector
    residual_jvp(const Vector &state_tangent,
                 const Vector &control_tangent) const
    {
      require_size(state_tangent, state_dimension(), "state tangent");
      require_size(control_tangent, control_dimension(), "control tangent");

      Vector value = apply_state_operator(
        state_tangent, ProblemBMatrixPurpose::residual_jvp);
      value.add(-1.0,
                mass_.coupling_apply(control_tangent,
                                     ProblemBMatrixPurpose::residual_jvp));
      return value;
    }

    ResidualDerivative
    residual_vjp(const Vector &test_seed) const
    {
      require_size(test_seed, state_dimension(), "test seed");

      Vector state = apply_transpose_state_operator(
        test_seed, ProblemBMatrixPurpose::residual_vjp);
      return {std::move(state),
              control_vjp(test_seed, ProblemBMatrixPurpose::residual_vjp)};
    }

    Vector
    control_vjp(const Vector &test_seed) const
    {
      return control_vjp(test_seed, ProblemBMatrixPurpose::control_vjp);
    }

    double
    objective(const Vector &state, const Vector &control) const
    {
      require_size(state, state_dimension(), "state");
      require_size(control, control_dimension(), "control");

      const Vector physical_state = coordinates_.reconstruct(state);
      const Vector state_mass =
        mass_.mass_apply(physical_state, ProblemBMatrixPurpose::objective);
      const Vector control_mass =
        mass_.mass_apply(control, ProblemBMatrixPurpose::objective);
      return 0.5 * (physical_state * state_mass) +
             0.5 * (control * control_mass);
    }

    ObjectiveDerivative
    objective_derivative(const Vector &state, const Vector &control) const
    {
      require_size(state, state_dimension(), "state");
      require_size(control, control_dimension(), "control");

      const Vector physical_state = coordinates_.reconstruct(state);
      const Vector state_mass = mass_.mass_apply(
        physical_state, ProblemBMatrixPurpose::objective_derivative);
      const Vector control_mass =
        mass_.mass_apply(control, ProblemBMatrixPurpose::objective_derivative);
      return {coordinates_.restrict(state_mass), control_mass};
    }

    SolveResult
    solve_state(const Vector &control) const
    {
      require_size(control, control_dimension(), "control");

      Vector rhs = application_.system_rhs_view();
      rhs.add(1.0,
              coordinates_.embed_free(
                mass_.coupling_apply(control, ProblemBMatrixPurpose::state_solve)));
      return solve_full_system(rhs);
    }

    SolveResult
    solve_adjoint(const Vector &state_objective_derivative) const
    {
      require_size(state_objective_derivative,
                   state_dimension(),
                   "state objective derivative");
      return solve_full_system(coordinates_.embed_free(
        state_objective_derivative));
    }

  private:
    Vector
    apply_state_operator(const Vector &                state,
                         const ProblemBMatrixPurpose purpose) const
    {
      Vector full_state = coordinates_.embed_free(state);
      Vector full_value(coordinates_.full_dimension());
      record(ProblemBMatrixAction::stiffness_apply, purpose);
      application_.system_matrix_view().vmult(full_value, full_state);
      return coordinates_.restrict(full_value);
    }

    Vector
    apply_transpose_state_operator(const Vector &                state,
                                   const ProblemBMatrixPurpose purpose) const
    {
      Vector full_state = coordinates_.embed_free(state);
      Vector full_value(coordinates_.full_dimension());
      record(ProblemBMatrixAction::stiffness_transpose_apply, purpose);
      application_.system_matrix_view().Tvmult(full_value, full_state);
      return coordinates_.restrict(full_value);
    }

    Vector
    control_vjp(const Vector &                test_seed,
                const ProblemBMatrixPurpose purpose) const
    {
      require_size(test_seed, state_dimension(), "control test seed");

      Vector control = mass_.coupling_transpose_apply(test_seed, purpose);
      control *= -1.0;
      return control;
    }

    void
    record(const ProblemBMatrixAction  action,
           const ProblemBMatrixPurpose purpose) const
    {
      if (instrumentation_ != nullptr)
        instrumentation_->record_problem_b_matrix_action(action, purpose);
    }

    SolveResult
    solve_full_system(const Vector &rhs) const
    {
      if (rhs.size() != coordinates_.full_dimension())
        throw std::invalid_argument(
          "Problem B full-system RHS has the wrong dimension");

      Vector full_solution(rhs.size());
      full_solution = 0.0;
      auto evidence = application_.solve(rhs, full_solution);
      return {coordinates_.restrict(full_solution),
              std::move(full_solution),
              std::move(evidence)};
    }

    static void
    require_size(const Vector & vector,
                 const std::size_t expected,
                 const char *const name)
    {
      if (static_cast<std::size_t>(vector.size()) != expected)
        throw std::invalid_argument(std::string("Problem B ") + name +
                                    " has the wrong dimension");
    }

    Application &application_;
    Instrumentation *instrumentation_;
    Coordinates  coordinates_;
    Mass         mass_;
    Vector       free_system_rhs_;
  };
} // namespace external_dealii_step4
