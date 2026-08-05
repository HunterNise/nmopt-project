#pragma once

#include "nmopt/contract/metric_constraint.hpp"
#include "nmopt/dealii/serial_backend.hpp"

#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>

namespace nmopt::dealii_backend
{
  // Coefficientwise bounds for a cellwise-constant serial decision space. The
  // owning lowerer is responsible for constructing this only for its declared
  // FE_DGQ(0) control or parameter layout.
  class CellwiseBoxConstraint final : public contract::ConstraintT<SerialBackend>
  {
  public:
    using Vector = dealii::Vector<double>;
    using Primal = contract::PrimalBlockT<SerialBackend>;

    CellwiseBoxConstraint(contract::LayoutPtr layout,
                          Vector              lower,
                          Vector              upper,
                          std::string l2_metric_id = "l2_cellwise")
      : layout_(std::move(layout))
      , lower_(std::move(lower))
      , upper_(std::move(upper))
      , l2_metric_id_(std::move(l2_metric_id))
    {
      contract::require(static_cast<bool>(layout_),
                        "Cellwise box constraint needs a layout");
      contract::require(layout_->n_blocks() == 1,
                        "Cellwise box constraint supports exactly one block");
      contract::require(lower_.size() == layout_->dimension(0) &&
                          upper_.size() == layout_->dimension(0),
                        "Cellwise box bounds do not match their layout");
      contract::require(!l2_metric_id_.empty(),
                        "Cellwise box metric identifier must not be empty");
      for (std::size_t index = 0; index < lower_.size(); ++index)
        contract::require(lower_[index] <= upper_[index],
                          "Cellwise box lower bound exceeds upper bound");
    }

    CellwiseBoxConstraint(contract::LayoutPtr layout,
                          const double        lower,
                          const double        upper,
                          std::string l2_metric_id = "l2_cellwise")
      : CellwiseBoxConstraint(layout,
                              constant_bound(layout, lower),
                              constant_bound(layout, upper),
                              std::move(l2_metric_id))
    {}

    const contract::LayoutPtr &
    layout() const override
    {
      return layout_;
    }

    bool
    is_feasible(const Primal &primal) const override
    {
      require_primal(primal, "Cellwise box primal has an incompatible layout");
      for (std::size_t index = 0; index < lower_.size(); ++index)
        if (primal.block(0)[index] < lower_[index] ||
            primal.block(0)[index] > upper_[index])
          return false;
      return true;
    }

    bool
    supports_projection_in(const contract::MetricT<SerialBackend> &metric) const override
    {
      return layout_->compatible_with(*metric.layout()) &&
             metric.id() == l2_metric_id_;
    }

    Primal
    project_in(const Primal &primal,
               const contract::MetricT<SerialBackend> &metric) const override
    {
      contract::require(supports_projection_in(metric),
                        "Cellwise box projection is unavailable for this metric");
      require_primal(primal, "Cellwise box primal has an incompatible layout");

      Vector projected = primal.block(0);
      for (std::size_t index = 0; index < projected.size(); ++index)
        projected[index] = std::max(lower_[index],
                                    std::min(upper_[index], projected[index]));
      return Primal(layout_, {std::move(projected)});
    }

  private:
    static Vector
    constant_bound(const contract::LayoutPtr &layout, const double value)
    {
      contract::require(static_cast<bool>(layout),
                        "Cellwise box constraint needs a layout");
      contract::require(layout->n_blocks() == 1,
                        "Cellwise box constraint supports exactly one block");

      Vector bound(layout->dimension(0));
      for (std::size_t index = 0; index < bound.size(); ++index)
        bound[index] = value;
      return bound;
    }

    void
    require_primal(const Primal &primal, const char *message) const
    {
      contract::require(primal.layout()->compatible_with(*layout_), message);
    }

    contract::LayoutPtr layout_;
    Vector              lower_;
    Vector              upper_;
    std::string         l2_metric_id_;
  };
} // namespace nmopt::dealii_backend
