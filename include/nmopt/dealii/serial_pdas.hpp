#pragma once

#include "nmopt/contract/complementarity.hpp"
#include "nmopt/dealii/mass_metric.hpp"
#include "nmopt/dealii/serial_backend.hpp"

#include <deal.II/lac/vector.h>

#include <memory>
#include <utility>

namespace nmopt::dealii_backend
{
  /**
   * Owning serial deal.II realization of a cellwise L2 box complementarity.
   *
   * The typed contract stores conversion callbacks, so the metric that owns
   * those callbacks must outlive the complementarity object. This bundle owns
   * the selected diagonal mass metric and exposes the resulting backend-
   * parameterized contract without reducing the dual to a display ID.
   */
  class SerialCellwiseBoxComplementarity final
  {
  public:
    using Vector = dealii::Vector<double>;
    using Primal = contract::PrimalBlockT<SerialBackend>;
    using Complementarity = contract::BoxComplementarityT<SerialBackend>;

    SerialCellwiseBoxComplementarity(
      contract::LayoutPtr          layout,
      Vector                       lower,
      Vector                       upper,
      std::shared_ptr<const MassMetric> metric)
      : metric_(std::move(metric))
      , complementarity_(make_complementarity(std::move(layout),
                                               std::move(lower),
                                               std::move(upper)))
    {}

    SerialCellwiseBoxComplementarity(
      contract::LayoutPtr          layout,
      const double                 lower,
      const double                 upper,
      std::shared_ptr<const MassMetric> metric)
      : SerialCellwiseBoxComplementarity(
          layout,
          constant_bound(layout, lower),
          constant_bound(layout, upper),
          std::move(metric))
    {}

    const MassMetric &
    metric() const
    {
      return *metric_;
    }

    const Complementarity &
    contract() const
    {
      return complementarity_;
    }

  private:
    static Vector
    constant_bound(const contract::LayoutPtr &layout, const double value)
    {
      contract::require(static_cast<bool>(layout),
                        "Serial cellwise box needs a layout");
      contract::require(layout->n_blocks() == 1,
                        "Serial cellwise box supports exactly one block");
      Vector bound(SerialBackend::checked_native_size(layout->dimension(0)));
      bound = value;
      return bound;
    }

    Complementarity
    make_complementarity(contract::LayoutPtr layout,
                         Vector            lower,
                         Vector            upper)
    {
      contract::require(static_cast<bool>(layout),
                        "Serial cellwise box needs a layout");
      contract::require(static_cast<bool>(metric_),
                        "Serial cellwise box needs an owned mass metric");
      contract::require(metric_->layout()->compatible_with(*layout),
                        "Serial cellwise box metric has an incompatible layout");
      contract::require(metric_->supports_coefficientwise_box_projection(),
                        "Serial cellwise box needs a positive diagonal L2 metric");
      const Primal lower_block(layout, {std::move(lower)});
      const Primal upper_block(layout, {std::move(upper)});
      return Complementarity(
        contract::BoxBoundsT<SerialBackend>(
          layout, lower_block, upper_block),
        contract::make_metric_multiplier_representation(
          std::shared_ptr<const contract::MetricT<SerialBackend>>(metric_)));
    }

    std::shared_ptr<const MassMetric> metric_;
    Complementarity                  complementarity_;
  };
} // namespace nmopt::dealii_backend
