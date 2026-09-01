#pragma once

#include "capability_registry.hpp"
#include "nmopt/application/scalar_function.hpp"

#include <stdexcept>
#include <string_view>

namespace nmopt::application::runner
{
  // Future benchmark binders can carry scalar lower and upper definitions
  // without introducing a B3/B4-specific function representation. This is a
  // data contract only; sampling a definition into a backend constraint is a
  // separate, not-yet-registered capability.
  struct BoxBoundDataSelection
  {
    ScalarFunctionDefinition lower;
    ScalarFunctionDefinition upper;
  };

  inline void
  validate_box_bound_data_selection(const BoxBoundDataSelection &selection)
  {
    validate_scalar_function_definition(selection.lower, "box lower bound");
    validate_scalar_function_definition(selection.upper, "box upper bound");
  }

  // This record names the implemented public quadratic-KKT contract without
  // claiming that the runner has a B5 adapter. Product, KKT method, and
  // preconditioner remain independently typed selections.
  struct QuadraticKktSolverSelection
  {
    chapter6::ProductSelection              product =
      chapter6::ProductSelection::quadratic_kkt;
    contract::QuadraticKKTSolverPolicy      solver;
    PreconditionerSelection                 preconditioner =
      PreconditionerSelection::identity_baseline;
  };

  inline QuadraticKktSolverSelection
  resolve_quadratic_kkt_solver_selection(
    const CapabilityRegistry<chapter6::ProductSelection> &product_registry,
    const CapabilityRegistry<contract::QuadraticKKTSolverMethod> &method_registry,
    const CapabilityRegistry<PreconditionerSelection> &preconditioner_registry,
    const std::string_view product_id,
    const std::string_view method_id,
    const std::string_view preconditioner_id)
  {
    QuadraticKktSolverSelection selection;
    selection.product = product_registry.resolve(product_id, "product");
    if (selection.product != chapter6::ProductSelection::quadratic_kkt)
      throw std::invalid_argument(
        "quadratic-KKT solver selection needs the quadratic-kkt product");
    selection.solver.method = method_registry.resolve(method_id, "KKT method");
    selection.preconditioner =
      preconditioner_registry.resolve(preconditioner_id, "preconditioner");
    if (!contract::valid(selection.solver))
      throw std::invalid_argument("quadratic-KKT solver policy is invalid");
    return selection;
  }

  inline QuadraticKktSolverSelection
  resolve_quadratic_kkt_solver_selection(
    const std::string_view product_id,
    const std::string_view method_id,
    const std::string_view preconditioner_id)
  {
    return resolve_quadratic_kkt_solver_selection(
      product_capability_registry(),
      quadratic_kkt_method_capability_registry(),
      preconditioner_capability_registry(),
      product_id,
      method_id,
      preconditioner_id);
  }
} // namespace nmopt::application::runner
