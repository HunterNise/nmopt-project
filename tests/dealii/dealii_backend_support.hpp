#pragma once

#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/compiler/v1/dealii_neumann_control_realisation.hpp"
#include "nmopt/dealii/cellwise_box_constraint.hpp"
#include "nmopt/dealii/facewise_box_constraint.hpp"
#include "nmopt/dealii/hminus1_metric.hpp"
#include "nmopt/dealii/quadratic_form.hpp"
#include "nmopt/dealii/serial_spd_solver.hpp"
#include "../support/contract_errors.hpp"
#include "../support/scenario_dispatch.hpp"

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  using namespace nmopt;
  using Backend = dealii_backend::SerialBackend;
  using Primal = contract::PrimalBlockT<Backend>;
  using Covector = contract::CovectorBlockT<Backend>;

  void
  require_close(const double              actual,
                const double              expected,
                const double              tolerance,
                const std::string &description)
  {
    if (std::abs(actual - expected) > tolerance)
      throw contract::ContractError(description + ": expected " +
                                    std::to_string(expected) + ", got " +
                                    std::to_string(actual));
  }
} // namespace
