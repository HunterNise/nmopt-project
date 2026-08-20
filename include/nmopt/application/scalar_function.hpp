#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>

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
} // namespace nmopt::application
