#pragma once

#include "nmopt/contract/layout.hpp"

namespace nmopt::contract
{
  template <typename Backend>
  class ExecutableModelT
  {
  public:
    virtual ~ExecutableModelT() = default;

    virtual const LayoutPtr &
    variable_layout() const = 0;

    virtual const LayoutPtr &
    test_layout() const = 0;

    virtual CovectorBlockT<Backend>
    residual(const PrimalBlockT<Backend> &variables) const = 0;

    virtual CovectorBlockT<Backend>
    residual_jvp(const PrimalBlockT<Backend> &variables,
                 const PrimalBlockT<Backend> &variable_tangent) const = 0;

    // The seed is primal in the residual test space. The result is a
    // covector in the model variable layout.
    virtual CovectorBlockT<Backend>
    residual_vjp(const PrimalBlockT<Backend> &variables,
                 const PrimalBlockT<Backend> &test_seed) const = 0;

    virtual double
    objective(const PrimalBlockT<Backend> &variables) const = 0;

    virtual CovectorBlockT<Backend>
    objective_derivative(const PrimalBlockT<Backend> &variables) const = 0;
  };

  using ExecutableModel = ExecutableModelT<DenseBackend>;
} // namespace nmopt::contract
