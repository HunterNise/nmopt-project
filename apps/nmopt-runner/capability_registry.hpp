#pragma once

#include "nmopt/application/chapter6.hpp"

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace nmopt::application::runner
{
  // A capability ID is the file-facing name of one typed selection. The
  // registry deliberately owns no construction or execution logic; it only
  // turns a declared name into the value consumed by a typed resolver.
  template <typename Selection>
  struct CapabilityDefinition
  {
    std::string_view id;
    Selection        selection;
  };

  template <typename Selection>
  class CapabilityRegistry final
  {
  public:
    using definition_type = CapabilityDefinition<Selection>;

    CapabilityRegistry() = default;

    CapabilityRegistry(
      const std::initializer_list<definition_type> definitions)
      : definitions_(definitions)
    {
      for (std::size_t index = 0; index < definitions_.size(); ++index)
        {
          if (definitions_[index].id.empty())
            throw std::logic_error("capability IDs must not be empty");
          for (std::size_t previous = 0; previous < index; ++previous)
            if (definitions_[previous].id == definitions_[index].id)
              throw std::logic_error("duplicate capability ID '" +
                                     std::string(definitions_[index].id) +
                                     "'");
        }
    }

    const definition_type *
    find(const std::string_view id) const
    {
      for (const auto &definition : definitions_)
        if (definition.id == id)
          return &definition;
      return nullptr;
    }

    Selection
    resolve(const std::string_view id,
            const std::string_view capability_name) const
    {
      const auto *definition = find(id);
      if (definition == nullptr)
        throw std::invalid_argument("unknown " + std::string(capability_name) +
                                    " '" + std::string(id) + "'");
      return definition->selection;
    }

    const std::vector<definition_type> &
    definitions() const
    {
      return definitions_;
    }

  private:
    std::vector<definition_type> definitions_;
  };

  enum class PreconditionerSelection
  {
    identity_baseline
  };

  inline const CapabilityRegistry<chapter6::ProductSelection> &
  product_capability_registry()
  {
    static const CapabilityRegistry<chapter6::ProductSelection> registry{
      {"reduced-dto", chapter6::ProductSelection::reduced_dto},
      {"quadratic-kkt", chapter6::ProductSelection::quadratic_kkt},
      {"pdas", chapter6::ProductSelection::pdas}};
    return registry;
  }

  inline const CapabilityRegistry<chapter6::ExecutionSelection> &
  execution_capability_registry()
  {
    static const CapabilityRegistry<chapter6::ExecutionSelection> registry{
      {"assembled", chapter6::ExecutionSelection::assembled},
      {"matrix-free", chapter6::ExecutionSelection::matrix_free}};
    return registry;
  }

  inline const CapabilityRegistry<chapter6::ReducedMethod> &
  reduced_method_capability_registry()
  {
    static const CapabilityRegistry<chapter6::ReducedMethod> registry{
      {"steepest-descent", chapter6::ReducedMethod::steepest_descent},
      {"bfgs", chapter6::ReducedMethod::bfgs},
      {"l-bfgs", chapter6::ReducedMethod::limited_memory_bfgs}};
    return registry;
  }

  inline const CapabilityRegistry<PreconditionerSelection> &
  preconditioner_capability_registry()
  {
    static const CapabilityRegistry<PreconditionerSelection> registry{
      {"identity", PreconditionerSelection::identity_baseline}};
    return registry;
  }
} // namespace nmopt::application::runner
