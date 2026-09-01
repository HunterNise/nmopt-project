#pragma once

#include "nmopt/application/scalar_function.hpp"

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/numbers.h>

#include <cctype>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace nmopt::application::dealii_support
{
  namespace detail
  {
    template <int dim>
    std::string
    scalar_coordinate_variable_names()
    {
      static_assert(dim > 0, "scalar functions need a positive dimension");
      std::ostringstream names;
      for (unsigned int coordinate = 0;
           coordinate < static_cast<unsigned int>(dim);
           ++coordinate)
        {
          if (coordinate != 0)
            names << ',';
          names << 'x' << coordinate;
        }
      return names.str();
    }

    inline bool
    scalar_expression_contains_identifier(const std::string_view expression,
                                          const std::string_view identifier)
    {
      const auto is_identifier_character = [](const char character) {
        const auto value = static_cast<unsigned char>(character);
        return std::isalnum(value) != 0 || character == '_';
      };
      for (auto position = expression.find(identifier);
           position != std::string_view::npos;
           position = expression.find(identifier, position + 1))
        {
          const bool starts_identifier =
            position == 0 || !is_identifier_character(expression[position - 1]);
          const auto end = position + identifier.size();
          const bool ends_identifier =
            end == expression.size() ||
            !is_identifier_character(expression[end]);
          if (starts_identifier && ends_identifier)
            return true;
        }
      return false;
    }
  } // namespace detail

  template <int dim>
  std::unique_ptr<::dealii::Function<dim>>
  make_scalar_function(const ScalarFunctionDefinition &definition,
                       const std::string_view           description =
                         "scalar function")
  {
    validate_scalar_function_definition(definition, description);
    switch (definition.kind)
      {
        case ScalarFunctionKind::zero:
          return std::make_unique<::dealii::Functions::ZeroFunction<dim>>();
        case ScalarFunctionKind::constant:
          return std::make_unique<
            ::dealii::Functions::ConstantFunction<dim>>(definition.value);
        case ScalarFunctionKind::expression:
          {
            if (detail::scalar_expression_contains_identifier(
                  definition.expression, "rand") ||
                detail::scalar_expression_contains_identifier(
                  definition.expression, "rand_seed"))
              throw std::invalid_argument(std::string(description) +
                                          " expressions cannot use random functions");
            auto parser = std::make_unique<::dealii::FunctionParser<dim>>();
            parser->initialize(
              detail::scalar_coordinate_variable_names<dim>(),
              definition.expression,
              {{"e", ::dealii::numbers::E},
               {"pi", ::dealii::numbers::PI}});
            ::dealii::Point<dim> validation_point;
            for (unsigned int coordinate = 0;
                 coordinate < static_cast<unsigned int>(dim);
                 ++coordinate)
              validation_point[coordinate] = 0.371 + 0.113 * coordinate;
            try
              {
                (void)parser->value(validation_point);
              }
            catch (const std::exception &exception)
              {
                throw std::invalid_argument(
                  std::string(description) + " expression is invalid: " +
                  exception.what());
              }
            return parser;
          }
      }
    throw std::invalid_argument(std::string(description) +
                                " has an unknown scalar function kind");
  }
} // namespace nmopt::application::dealii_support
