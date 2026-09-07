#pragma once

#include "nmopt/contract/executable_model.hpp"

#include <functional>
#include <string>
#include <utility>

namespace nmopt::contract
{
  template <typename Backend>
  class CallbackExecutableModelT final : public ExecutableModelT<Backend>
  {
  public:
    using Primal = PrimalBlockT<Backend>;
    using Covector = CovectorBlockT<Backend>;

    using ResidualAction = std::function<Covector(const Primal &)>;
    using LinearizedAction =
      std::function<Covector(const Primal &, const Primal &)>;
    using TransposeAction =
      std::function<Covector(const Primal &, const Primal &)>;
    using ObjectiveAction = std::function<double(const Primal &)>;
    using ObjectiveDerivativeAction = std::function<Covector(const Primal &)>;

    CallbackExecutableModelT(LayoutPtr                    variable_layout,
                             LayoutPtr                    test_layout,
                             ResidualAction                residual,
                             LinearizedAction              residual_jvp,
                             TransposeAction               residual_vjp,
                             ObjectiveAction               objective,
                             ObjectiveDerivativeAction     objective_derivative)
      : variable_layout_(std::move(variable_layout))
      , test_layout_(std::move(test_layout))
      , residual_(std::move(residual))
      , residual_jvp_(std::move(residual_jvp))
      , residual_vjp_(std::move(residual_vjp))
      , objective_(std::move(objective))
      , objective_derivative_(std::move(objective_derivative))
    {
      require(static_cast<bool>(variable_layout_),
              "Callback executable model needs a variable layout");
      require(static_cast<bool>(test_layout_),
              "Callback executable model needs a test layout");
      require(static_cast<bool>(residual_),
              "Callback executable model needs a residual action");
      require(static_cast<bool>(residual_jvp_),
              "Callback executable model needs a residual JVP action");
      require(static_cast<bool>(residual_vjp_),
              "Callback executable model needs a residual VJP action");
      require(static_cast<bool>(objective_),
              "Callback executable model needs an objective action");
      require(static_cast<bool>(objective_derivative_),
              "Callback executable model needs an objective derivative action");
    }

    const LayoutPtr &
    variable_layout() const override
    {
      return variable_layout_;
    }

    const LayoutPtr &
    test_layout() const override
    {
      return test_layout_;
    }

    Covector
    residual(const Primal &variables) const override
    {
      require_variable(variables, "Callback executable model residual");
      return checked_test_covector(residual_(variables),
                                   "Callback executable model residual");
    }

    Covector
    residual_jvp(const Primal &variables,
                 const Primal &variable_tangent) const override
    {
      require_variable(variables,
                       "Callback executable model residual JVP point");
      require_variable(variable_tangent,
                       "Callback executable model residual JVP tangent");
      return checked_test_covector(
        residual_jvp_(variables, variable_tangent),
        "Callback executable model residual JVP");
    }

    Covector
    residual_vjp(const Primal &variables,
                 const Primal &test_seed) const override
    {
      require_variable(variables,
                       "Callback executable model residual VJP point");
      require(test_seed.layout()->compatible_with(*test_layout_),
              "Callback executable model residual VJP seed has an incompatible test layout");
      return checked_variable_covector(
        residual_vjp_(variables, test_seed),
        "Callback executable model residual VJP");
    }

    double
    objective(const Primal &variables) const override
    {
      require_variable(variables, "Callback executable model objective");
      return objective_(variables);
    }

    Covector
    objective_derivative(const Primal &variables) const override
    {
      require_variable(variables,
                       "Callback executable model objective derivative");
      return checked_variable_covector(
        objective_derivative_(variables),
        "Callback executable model objective derivative");
    }

  private:
    void
    require_variable(const Primal &value, const char *operation) const
    {
      require(value.layout()->compatible_with(*variable_layout_),
              std::string(operation) +
                " has an incompatible variable layout");
    }

    Covector
    checked_test_covector(Covector value, const char *operation) const
    {
      require(value.layout()->compatible_with(*test_layout_),
              std::string(operation) +
                " returned an incompatible test layout");
      return value;
    }

    Covector
    checked_variable_covector(Covector value, const char *operation) const
    {
      require(value.layout()->compatible_with(*variable_layout_),
              std::string(operation) +
                " returned an incompatible variable layout");
      return value;
    }

    LayoutPtr                variable_layout_;
    LayoutPtr                test_layout_;
    ResidualAction           residual_;
    LinearizedAction         residual_jvp_;
    TransposeAction          residual_vjp_;
    ObjectiveAction          objective_;
    ObjectiveDerivativeAction objective_derivative_;
  };

  using CallbackExecutableModel = CallbackExecutableModelT<DenseBackend>;
} // namespace nmopt::contract
