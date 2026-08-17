#pragma once

#include "nmopt/contract/linalg.hpp"

#include <string>
#include <utility>

namespace nmopt::test_support
{
  template <typename Callable>
  void
  require_contract_error(Callable &&        callable,
                         const std::string &expected_message,
                         const std::string &description)
  {
    try
      {
        std::forward<Callable>(callable)();
      }
    catch (const contract::ContractError &exception)
      {
        contract::require(exception.what() == expected_message,
                          description + ": expected ContractError \"" +
                            expected_message + "\", got \"" +
                            exception.what() + "\"");
        return;
      }
    throw contract::ContractError(description +
                                  ": expected ContractError was not thrown");
  }
} // namespace nmopt::test_support
