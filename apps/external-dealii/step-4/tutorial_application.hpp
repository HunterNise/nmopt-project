#pragma once

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <cstddef>
#include <filesystem>
#include <memory>

class TutorialApplication final
{
public:
  using Matrix = dealii::SparseMatrix<double>;
  using Vector = dealii::Vector<double>;

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

  TutorialApplication();
  ~TutorialApplication();

  TutorialApplication(TutorialApplication &&) noexcept;
  TutorialApplication &operator=(TutorialApplication &&) noexcept;

  TutorialApplication(const TutorialApplication &) = delete;
  TutorialApplication &operator=(const TutorialApplication &) = delete;

  std::size_t state_dimension() const;
  std::size_t control_dimension() const;

  Vector residual(const Vector &state, const Vector &control) const;
  Vector residual_jvp(const Vector &state_tangent,
                      const Vector &control_tangent) const;
  ResidualDerivative residual_vjp(const Vector &test_seed) const;

  double objective(const Vector &state, const Vector &control) const;
  ObjectiveDerivative objective_derivative(const Vector &state,
                                           const Vector &control) const;

  Vector solve_state(const Vector &control) const;
  Vector solve_adjoint(const Vector &state_objective_derivative) const;

  std::shared_ptr<const Matrix> control_metric_matrix() const;

  void write_native_output(const std::filesystem::path &directory,
                           const Vector &                 state) const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
