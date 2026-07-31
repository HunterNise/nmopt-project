#pragma once

#include "nmopt/contract/layout.hpp"

#include <algorithm>
#include <string>
#include <utility>

namespace nmopt::contract
{
  template <typename Backend>
  class MetricT
  {
  public:
    virtual ~MetricT() = default;

    virtual const std::string &
    id() const = 0;

    virtual const LayoutPtr &
    layout() const = 0;

    virtual CovectorBlockT<Backend>
    apply(const PrimalBlockT<Backend> &primal) const = 0;

    virtual PrimalBlockT<Backend>
    inverse_apply(const CovectorBlockT<Backend> &covector) const = 0;
  };

  using Metric = MetricT<DenseBackend>;

  class DiagonalMetric final : public Metric
  {
  public:
    DiagonalMetric(std::string              id,
                   LayoutPtr                layout,
                   std::vector<DenseVector> diagonal)
      : id_(std::move(id))
      , layout_(std::move(layout))
      , diagonal_(std::move(diagonal))
    {
      require(!id_.empty(), "Metric identifier must not be empty");
      const PrimalBlock probe(layout_, diagonal_);
      for (std::size_t block = 0; block < probe.n_blocks(); ++block)
        for (std::size_t entry = 0; entry < probe.block(block).size(); ++entry)
          require(probe.block(block)[entry] > 0.0,
                  "Metric diagonal must be strictly positive");
    }

    const std::string &
    id() const override
    {
      return id_;
    }

    const LayoutPtr &
    layout() const override
    {
      return layout_;
    }

    CovectorBlock
    apply(const PrimalBlock &primal) const override
    {
      require_compatible(primal, PrimalBlock(layout_, diagonal_),
                         "Metric primal has an incompatible layout");
      CovectorBlock result = CovectorBlock::zeros(layout_);
      for (std::size_t block = 0; block < primal.n_blocks(); ++block)
        for (std::size_t entry = 0; entry < primal.block(block).size(); ++entry)
          result.block(block)[entry] =
            diagonal_.at(block)[entry] * primal.block(block)[entry];
      return result;
    }

    PrimalBlock
    inverse_apply(const CovectorBlock &covector) const override
    {
      require_compatible(covector, PrimalBlock(layout_, diagonal_),
                         "Metric covector has an incompatible layout");
      PrimalBlock result = PrimalBlock::zeros(layout_);
      for (std::size_t block = 0; block < covector.n_blocks(); ++block)
        for (std::size_t entry = 0; entry < covector.block(block).size(); ++entry)
          result.block(block)[entry] =
            covector.block(block)[entry] / diagonal_.at(block)[entry];
      return result;
    }

  private:
    std::string              id_;
    LayoutPtr                layout_;
    std::vector<DenseVector> diagonal_;
  };

  template <typename Backend>
  class ConstraintT
  {
  public:
    virtual ~ConstraintT() = default;

    virtual const LayoutPtr &
    layout() const = 0;

    virtual bool
    is_feasible(const PrimalBlockT<Backend> &primal) const = 0;

    virtual bool
    supports_projection_in(const MetricT<Backend> &metric) const = 0;

    virtual PrimalBlockT<Backend>
    project_in(const PrimalBlockT<Backend> &primal,
               const MetricT<Backend>      &metric) const = 0;
  };

  using Constraint = ConstraintT<DenseBackend>;

  // This is the first-default box policy: every coefficient represents one
  // cellwise-constant control value, and the metric is the declared L2 metric.
  class CellwiseBoxConstraint final : public Constraint
  {
  public:
    CellwiseBoxConstraint(LayoutPtr                layout,
                          std::vector<DenseVector> lower,
                          std::vector<DenseVector> upper,
                          std::string              l2_metric_id = "l2_cellwise")
      : layout_(std::move(layout))
      , lower_(std::move(lower))
      , upper_(std::move(upper))
      , l2_metric_id_(std::move(l2_metric_id))
    {
      const PrimalBlock lower_probe(layout_, lower_);
      const PrimalBlock upper_probe(layout_, upper_);
      for (std::size_t block = 0; block < lower_probe.n_blocks(); ++block)
        for (std::size_t entry = 0; entry < lower_probe.block(block).size();
             ++entry)
          require(lower_probe.block(block)[entry] <= upper_probe.block(block)[entry],
                  "Cellwise box lower bound exceeds upper bound");
    }

    const LayoutPtr &
    layout() const override
    {
      return layout_;
    }

    bool
    is_feasible(const PrimalBlock &primal) const override
    {
      require_compatible(primal, PrimalBlock(layout_, lower_),
                         "Constraint primal has an incompatible layout");
      for (std::size_t block = 0; block < primal.n_blocks(); ++block)
        for (std::size_t entry = 0; entry < primal.block(block).size(); ++entry)
          if (primal.block(block)[entry] < lower_.at(block)[entry] ||
              primal.block(block)[entry] > upper_.at(block)[entry])
            return false;
      return true;
    }

    bool
    supports_projection_in(const Metric &metric) const override
    {
      return layout_->compatible_with(*metric.layout()) &&
             metric.id() == l2_metric_id_;
    }

    PrimalBlock
    project_in(const PrimalBlock &primal, const Metric &metric) const override
    {
      require(supports_projection_in(metric),
              "Cellwise box projection is unavailable for this metric");
      require_compatible(primal, PrimalBlock(layout_, lower_),
                         "Constraint primal has an incompatible layout");

      PrimalBlock projected = primal;
      for (std::size_t block = 0; block < primal.n_blocks(); ++block)
        for (std::size_t entry = 0; entry < primal.block(block).size(); ++entry)
          projected.block(block)[entry] =
            std::max(lower_.at(block)[entry],
                     std::min(upper_.at(block)[entry],
                              primal.block(block)[entry]));
      return projected;
    }

  private:
    LayoutPtr                layout_;
    std::vector<DenseVector> lower_;
    std::vector<DenseVector> upper_;
    std::string              l2_metric_id_;
  };
} // namespace nmopt::contract
