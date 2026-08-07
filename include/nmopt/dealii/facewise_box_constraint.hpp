#pragma once

#include "nmopt/dealii/mass_metric.hpp"
#include "nmopt/dealii/serial_backend.hpp"

#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>

namespace nmopt::dealii_backend
{
  // Coefficientwise bounds for a facewise-constant boundary control. The
  // owning lowerer creates one coefficient for each selected boundary face;
  // unlike CellwiseBoxConstraint, this class has no volume-cell semantics.
  class FacewiseBoxConstraint final : public contract::ConstraintT<SerialBackend>
  {
  public:
    using Vector = dealii::Vector<double>;
    using Primal = contract::PrimalBlockT<SerialBackend>;

    FacewiseBoxConstraint(contract::LayoutPtr layout,
                          Vector              lower,
                          Vector              upper,
                          const MassMetric &  projection_metric)
      : layout_(std::move(layout))
      , lower_(std::move(lower))
      , upper_(std::move(upper))
      , projection_metric_(projection_metric.realisation_witness())
    {
      contract::require(static_cast<bool>(layout_),
                        "Facewise box constraint needs a layout");
      contract::require(layout_->n_blocks() == 1,
                        "Facewise box constraint supports exactly one block");
      contract::require(lower_.size() == layout_->dimension(0) &&
                          upper_.size() == layout_->dimension(0),
                        "Facewise box bounds do not match their layout");
      contract::require(layout_->compatible_with(*projection_metric.layout()),
                        "Facewise box projection metric has an incompatible layout");
      contract::require(
        projection_metric.supports_coefficientwise_box_projection(),
        "Facewise box projection needs a positive diagonal metric realization");
      for (std::size_t index = 0; index < lower_.size(); ++index)
        contract::require(lower_[index] <= upper_[index],
                          "Facewise box lower bound exceeds upper bound");
    }

    FacewiseBoxConstraint(contract::LayoutPtr layout,
                          const double        lower,
                          const double        upper,
                          const MassMetric &  projection_metric)
      : FacewiseBoxConstraint(layout,
                              constant_bound(layout, lower),
                              constant_bound(layout, upper),
                              projection_metric)
    {}

    const contract::LayoutPtr &
    layout() const override
    {
      return layout_;
    }

    bool
    is_feasible(const Primal &primal) const override
    {
      require_primal(primal, "Facewise box primal has an incompatible layout");
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
             projection_metric_.matches(metric.realisation_witness());
    }

    Primal
    project_in(const Primal &primal,
               const contract::MetricT<SerialBackend> &metric) const override
    {
      contract::require(supports_projection_in(metric),
                        "Facewise box projection is unavailable for this metric");
      require_primal(primal, "Facewise box primal has an incompatible layout");

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
                        "Facewise box constraint needs a layout");
      contract::require(layout->n_blocks() == 1,
                        "Facewise box constraint supports exactly one block");

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
    contract::MetricRealisationWitness projection_metric_;
  };
} // namespace nmopt::dealii_backend
