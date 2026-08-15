#pragma once

#include "nmopt/contract/layout.hpp"

namespace nmopt::contract
{
  // A reduced Hessian is an explicit second-order capability. First-order
  // residual JVP/VJP and objective-derivative ports do not satisfy this
  // interface by implication.
  template <typename Backend>
  class ReducedHessianT
  {
  public:
    virtual ~ReducedHessianT() = default;

    virtual const LayoutPtr &
    layout() const = 0;

    // Apply H(u) to a primal control direction and return the reduced
    // covector H(u)w. The provider owns any tangent-state and
    // incremental-adjoint work needed to produce the action.
    virtual CovectorBlockT<Backend>
    apply(const PrimalBlockT<Backend> &control,
          const PrimalBlockT<Backend> &direction) const = 0;
  };

  using ReducedHessian = ReducedHessianT<DenseBackend>;
} // namespace nmopt::contract
