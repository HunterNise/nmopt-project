#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace nmopt::application
{
  enum class ScalarFunctionKind
  {
    zero,
    constant,
    expression
  };

  inline const char *
  scalar_function_kind_name(const ScalarFunctionKind kind)
  {
    switch (kind)
      {
        case ScalarFunctionKind::zero:
          return "zero";
        case ScalarFunctionKind::constant:
          return "constant";
        case ScalarFunctionKind::expression:
          return "expression";
      }
    throw std::invalid_argument("unknown scalar function kind");
  }

  inline ScalarFunctionKind
  scalar_function_kind_from_name(const std::string_view name)
  {
    if (name == "zero")
      return ScalarFunctionKind::zero;
    if (name == "constant")
      return ScalarFunctionKind::constant;
    if (name == "expression")
      return ScalarFunctionKind::expression;
    throw std::invalid_argument("unknown scalar function kind '" +
                                std::string(name) + "'");
  }

  struct ScalarFunctionDefinition
  {
    std::string        id;
    ScalarFunctionKind kind = ScalarFunctionKind::zero;
    double             value = 0.0;
    std::string        expression;
    std::string        provenance;
  };

  struct ScalarFunctionCatalog
  {
    std::vector<ScalarFunctionDefinition> definitions;
    std::string                           selected_id;
  };

  inline bool
  operator==(const ScalarFunctionDefinition &lhs,
             const ScalarFunctionDefinition &rhs)
  {
    return lhs.id == rhs.id && lhs.kind == rhs.kind &&
           lhs.value == rhs.value && lhs.expression == rhs.expression &&
           lhs.provenance == rhs.provenance;
  }

  inline bool
  operator!=(const ScalarFunctionDefinition &lhs,
             const ScalarFunctionDefinition &rhs)
  {
    return !(lhs == rhs);
  }

  inline const ScalarFunctionDefinition &
  selected_scalar_function_definition(
    const ScalarFunctionCatalog &catalog,
    const std::string_view        description = "scalar function catalog")
  {
    if (catalog.selected_id.empty())
      throw std::invalid_argument(std::string(description) +
                                  " needs a selected definition ID");
    const auto selected = std::find_if(
      catalog.definitions.begin(),
      catalog.definitions.end(),
      [&](const auto &definition) {
        return definition.id == catalog.selected_id;
      });
    if (selected == catalog.definitions.end())
      throw std::invalid_argument(std::string(description) +
                                  " selected definition ID is not registered");
    return *selected;
  }

  inline void
  validate_scalar_function_definition(
    const ScalarFunctionDefinition &definition,
    const std::string_view           description = "scalar function")
  {
    const auto prefix = std::string(description);
    if (definition.id.empty())
      throw std::invalid_argument(prefix + " needs a stable ID");
    if (definition.provenance.empty())
      throw std::invalid_argument(prefix + " needs provenance");

    switch (definition.kind)
      {
        case ScalarFunctionKind::zero:
          if (definition.value != 0.0 || !definition.expression.empty())
            throw std::invalid_argument(
              prefix + " zero kind cannot carry a value or expression");
          break;
        case ScalarFunctionKind::constant:
          if (!std::isfinite(definition.value))
            throw std::invalid_argument(
              prefix + " constant value must be finite");
          if (!definition.expression.empty())
            throw std::invalid_argument(
              prefix + " constant kind cannot carry an expression");
          break;
        case ScalarFunctionKind::expression:
          if (!std::isfinite(definition.value) || definition.value != 0.0)
            throw std::invalid_argument(
              prefix + " expression kind cannot carry a value");
          if (definition.expression.empty())
            throw std::invalid_argument(
              prefix + " expression kind needs an expression");
          if (definition.expression.find(';') != std::string::npos)
            throw std::invalid_argument(
              prefix + " must contain one scalar expression");
          break;
        default:
          throw std::invalid_argument(prefix + " has an unknown kind");
      }
  }

  inline void
  validate_scalar_function_catalog(
    const ScalarFunctionCatalog &catalog,
    const std::string_view        description = "scalar function catalog")
  {
    if (catalog.definitions.empty())
      throw std::invalid_argument(std::string(description) +
                                  " needs at least one definition");
    for (std::size_t index = 0; index < catalog.definitions.size(); ++index)
      {
        validate_scalar_function_definition(catalog.definitions[index],
                                            description);
        for (std::size_t previous = 0; previous < index; ++previous)
          if (catalog.definitions[previous].id == catalog.definitions[index].id)
            throw std::invalid_argument(std::string(description) +
                                        " contains duplicate definition ID '" +
                                        catalog.definitions[index].id + "'");
      }
    (void)selected_scalar_function_definition(catalog, description);
  }
} // namespace nmopt::application
