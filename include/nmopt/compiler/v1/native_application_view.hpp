#pragma once

#include "nmopt/contract/layout.hpp"

#include <cstddef>
#include <filesystem>
#include <functional>
#include <utility>

namespace nmopt::compiler::v1
{
  template <typename Backend>
  struct NativeApplicationDimensionsT
  {
    std::size_t physical_state = 0;
    std::size_t independent_state = 0;
    std::size_t physical_control = 0;
    std::size_t independent_control = 0;
    std::size_t realized_observation = 0;
  };

  template <typename Backend>
  struct NativeObjectiveComponentsT
  {
    double state_tracking = 0.0;
    double control_regularisation = 0.0;
  };

  // Compiler-path application access retained beside the erased solver view.
  // This is deliberately a small callback/value seam, not a public numerical
  // realization hierarchy or an output contract on ExecutableModelT. The
  // callback may borrow the compile-time data bindings, so callers must keep
  // those bindings alive while invoking the view; the retained model itself
  // is owned by the view.
  template <typename Backend>
  class NativeApplicationViewT final
  {
  public:
    using Primal = contract::PrimalBlockT<Backend>;
    using Dimensions = NativeApplicationDimensionsT<Backend>;
    using ObjectiveComponents = NativeObjectiveComponentsT<Backend>;
    using OutputAction = std::function<void(const std::filesystem::path &,
                                            const Primal &,
                                            const Primal &,
                                            const Primal &,
                                            const Primal *)>;
    using ObjectiveComponentsAction =
      std::function<ObjectiveComponents(const Primal &)>;

    NativeApplicationViewT(
      Dimensions               dimensions,
      OutputAction              output,
      ObjectiveComponentsAction objective_components = {})
      : dimensions_(dimensions)
      , output_(std::move(output))
      , objective_components_(std::move(objective_components))
    {
      contract::require(dimensions_.physical_state > 0,
                        "Native application view needs a physical state dimension");
      contract::require(dimensions_.independent_state > 0,
                        "Native application view needs an independent state dimension");
      contract::require(dimensions_.physical_control > 0,
                        "Native application view needs a physical control dimension");
      contract::require(dimensions_.independent_control > 0,
                        "Native application view needs an independent control dimension");
      contract::require(static_cast<bool>(output_),
                        "Native application view needs an output action");
    }

    const Dimensions &
    dimensions() const noexcept
    {
      return dimensions_;
    }

    bool
    has_objective_components() const noexcept
    {
      return static_cast<bool>(objective_components_);
    }

    ObjectiveComponents
    objective_components(const Primal &full_point) const
    {
      contract::require(static_cast<bool>(objective_components_),
                        "Native application view has no objective components");
      return objective_components_(full_point);
    }

    void
    write_native_output(const std::filesystem::path &directory,
                        const Primal &                 state,
                        const Primal &                 control,
                        const Primal &                 adjoint,
                        const Primal *                 uncontrolled_state = nullptr) const
    {
      output_(directory, state, control, adjoint, uncontrolled_state);
    }

  private:
    Dimensions                 dimensions_;
    OutputAction               output_;
    ObjectiveComponentsAction  objective_components_;
  };

  template <typename Backend>
  using NativeApplicationView = NativeApplicationViewT<Backend>;
} // namespace nmopt::compiler::v1
