#pragma once

#include "nmopt/contract/metric_constraint.hpp"

#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::contract
{
  enum class BoxActivity
  {
    lower,
    inactive,
    upper
  };

  template <typename Values>
  void
  require_complementarity_layout(const Values &value,
                                 const LayoutPtr &expected,
                                 const std::string &operation)
  {
    require(value.layout()->compatible_with(*expected),
            operation + " has an incompatible layout");
  }

  template <typename Backend>
  struct BoxMultiplierRepresentationT
  {
    using Primal = PrimalBlockT<Backend>;
    using Covector = CovectorBlockT<Backend>;
    using ToPrimalAction = std::function<Primal(const Covector &)>;
    using ToDualAction = std::function<Covector(const Primal &)>;

    LayoutPtr       primal_layout;
    LayoutPtr       dual_layout;
    std::string     description;
    ToPrimalAction  to_primal;
    ToDualAction    to_dual;
    std::shared_ptr<const MetricT<Backend>> metric_owner;
    MetricRealisationWitness metric_witness;
  };

  using BoxMultiplierRepresentation =
    BoxMultiplierRepresentationT<DenseBackend>;

  template <typename Backend>
  BoxMultiplierRepresentationT<Backend>
  make_metric_multiplier_representation(
    std::shared_ptr<const MetricT<Backend>> metric)
  {
    using Representation = BoxMultiplierRepresentationT<Backend>;
    using Primal = typename Representation::Primal;
    using Covector = typename Representation::Covector;
    require(static_cast<bool>(metric),
            "Metric multiplier representation needs an owned metric");
    const auto metric_witness = metric->realisation_witness();

    return Representation{
      metric->layout(),
      metric->layout(),
      "metric Riesz map '" + metric->id() + "'",
      [metric](const Covector &covector) {
        return metric->inverse_apply(covector);
      },
      [metric](const Primal &primal) { return metric->apply(primal); },
      std::move(metric),
      metric_witness};
  }

  template <typename Backend>
  class BoxBoundsT final
  {
  public:
    using Primal = PrimalBlockT<Backend>;

    BoxBoundsT(LayoutPtr layout, Primal lower, Primal upper)
      : layout_(std::move(layout))
      , lower_(std::move(lower))
      , upper_(std::move(upper))
    {
      require(static_cast<bool>(layout_), "Box bounds need a layout");
      require_complementarity_layout(lower_, layout_, "Box lower bounds");
      require_complementarity_layout(upper_, layout_, "Box upper bounds");
      for (std::size_t block = 0; block < layout_->n_blocks(); ++block)
        for (std::size_t index = 0; index < layout_->dimension(block); ++index)
          {
            const double lower_value =
              Backend::value(lower_.block(block), index);
            const double upper_value =
              Backend::value(upper_.block(block), index);
            require(std::isfinite(lower_value) && std::isfinite(upper_value),
                    "Box bounds must be finite");
            require(lower_value <= upper_value,
                    "Box lower bound exceeds upper bound");
          }
    }

    const LayoutPtr &
    layout() const
    {
      return layout_;
    }

    const Primal &
    lower() const
    {
      return lower_;
    }

    const Primal &
    upper() const
    {
      return upper_;
    }

    bool
    is_feasible(const Primal &primal) const
    {
      require_complementarity_layout(primal,
                                     layout_,
                                     "Box feasibility primal");
      for (std::size_t block = 0; block < layout_->n_blocks(); ++block)
        for (std::size_t index = 0; index < layout_->dimension(block); ++index)
          {
            const double value = Backend::value(primal.block(block), index);
            if (value < Backend::value(lower_.block(block), index) ||
                value > Backend::value(upper_.block(block), index))
              return false;
          }
      return true;
    }

  private:
    LayoutPtr layout_;
    Primal    lower_;
    Primal    upper_;
  };

  using BoxBounds = BoxBoundsT<DenseBackend>;

  template <typename Backend>
  class ActiveSetSelectionT final
  {
  public:
    using Vector = typename Backend::Vector;

    ActiveSetSelectionT(LayoutPtr                    layout,
                        std::vector<BoxActivity>     activities)
      : layout_(std::move(layout))
      , activities_(std::move(activities))
    {
      require(static_cast<bool>(layout_),
              "Active-set selection needs a control layout");
      require(layout_->n_blocks() == 1,
              "Active-set selection supports exactly one control block");
      require(activities_.size() == layout_->dimension(0),
              "Active-set selection size does not match its layout");

      for (std::size_t index = 0; index < activities_.size(); ++index)
        {
          if (activities_[index] == BoxActivity::inactive)
            free_indices_.push_back(index);
          else
            active_indices_.push_back(index);
        }
    }

    const LayoutPtr &
    layout() const
    {
      return layout_;
    }

    const std::vector<BoxActivity> &
    activities() const
    {
      return activities_;
    }

    const std::vector<std::size_t> &
    free_indices() const
    {
      return free_indices_;
    }

    const std::vector<std::size_t> &
    active_indices() const
    {
      return active_indices_;
    }

    std::size_t
    free_size() const
    {
      return free_indices_.size();
    }

    std::size_t
    active_size() const
    {
      return active_indices_.size();
    }

    bool
    operator==(const ActiveSetSelectionT &other) const
    {
      return layout_->compatible_with(*other.layout_) &&
             activities_ == other.activities_;
    }

    Vector
    restrict_free(const Vector &full) const
    {
      require_vector_size(full, "Active-set free restriction input");
      Vector result = Backend::zeros(free_indices_.size());
      for (std::size_t selected = 0; selected < free_indices_.size(); ++selected)
        Backend::set_value(result,
                           selected,
                           Backend::value(full, free_indices_[selected]));
      return result;
    }

    Vector
    restrict_active(const Vector &full) const
    {
      require_vector_size(full, "Active-set active restriction input");
      Vector result = Backend::zeros(active_indices_.size());
      for (std::size_t selected = 0; selected < active_indices_.size(); ++selected)
        Backend::set_value(result,
                           selected,
                           Backend::value(full, active_indices_[selected]));
      return result;
    }

    Vector
    prolong(const Vector &free, const Vector &active) const
    {
      require(free.size() == free_indices_.size(),
              "Active-set prolongation free vector has the wrong size");
      require(active.size() == active_indices_.size(),
              "Active-set prolongation active vector has the wrong size");

      Vector result = Backend::zeros(layout_->dimension(0));
      for (std::size_t selected = 0; selected < free_indices_.size(); ++selected)
        Backend::set_value(result,
                           free_indices_[selected],
                           Backend::value(free, selected));
      for (std::size_t selected = 0; selected < active_indices_.size(); ++selected)
        Backend::set_value(result,
                           active_indices_[selected],
                           Backend::value(active, selected));
      return result;
    }

  private:
    void
    require_vector_size(const Vector &vector, const char *operation) const
    {
      require(Backend::size(vector) == layout_->dimension(0),
              std::string(operation) + " has the wrong size");
    }

    LayoutPtr               layout_;
    std::vector<BoxActivity> activities_;
    std::vector<std::size_t> free_indices_;
    std::vector<std::size_t> active_indices_;
  };

  using ActiveSetSelection = ActiveSetSelectionT<DenseBackend>;

  template <typename Backend>
  class BoxComplementarityT final
  {
  public:
    using Primal = PrimalBlockT<Backend>;
    using Covector = CovectorBlockT<Backend>;
    using Selection = ActiveSetSelectionT<Backend>;
    using Bounds = BoxBoundsT<Backend>;
    using Representation = BoxMultiplierRepresentationT<Backend>;

    BoxComplementarityT(Bounds bounds, Representation representation)
      : bounds_(std::move(bounds))
      , representation_(std::move(representation))
    {
      require(static_cast<bool>(representation_.primal_layout),
              "Box multiplier representation needs a primal layout");
      require(static_cast<bool>(representation_.dual_layout),
              "Box multiplier representation needs a dual layout");
      require(static_cast<bool>(representation_.metric_owner),
              "Box multiplier representation needs an owned metric");
      require(bounds_.layout()->compatible_with(
                *representation_.primal_layout),
              "Box multiplier primal representation has an incompatible layout");
      require(bounds_.layout()->compatible_with(
                *representation_.dual_layout),
              "Box multiplier dual representation has an incompatible layout");
      require(!representation_.description.empty(),
              "Box multiplier representation needs a description");
      require(static_cast<bool>(representation_.to_primal),
              "Box multiplier representation needs a dual-to-primal action");
      require(static_cast<bool>(representation_.to_dual),
              "Box multiplier representation needs a primal-to-dual action");
      require(representation_.metric_owner->layout()->compatible_with(
                *representation_.primal_layout),
              "Box multiplier representation metric has an incompatible layout");
      require(representation_.metric_witness.matches(
                representation_.metric_owner->realisation_witness()),
              "Box multiplier representation metric witness does not match its owner");
    }

    const LayoutPtr &
    layout() const
    {
      return bounds_.layout();
    }

    const Bounds &
    bounds() const
    {
      return bounds_;
    }

    const Representation &
    multiplier_representation() const
    {
      return representation_;
    }

    Primal
    multiplier_to_primal(const Covector &multiplier) const
    {
      require_complementarity_layout(multiplier,
                                     representation_.dual_layout,
                                     "Box multiplier dual-to-primal input");
      Primal result = representation_.to_primal(multiplier);
      require_complementarity_layout(result,
                                     bounds_.layout(),
                                     "Box multiplier dual-to-primal result");
      return result;
    }

    Covector
    primal_to_multiplier(const Primal &primal) const
    {
      require_complementarity_layout(primal,
                                     bounds_.layout(),
                                     "Box multiplier primal-to-dual input");
      Covector result = representation_.to_dual(primal);
      require_complementarity_layout(result,
                                     representation_.dual_layout,
                                     "Box multiplier primal-to-dual result");
      return result;
    }

    Selection
    classify(const Primal &primal,
             const Covector &multiplier,
             const double c) const
    {
      require_complementarity_layout(primal,
                                     bounds_.layout(),
                                     "Box classification primal");
      require(std::isfinite(c) && c > 0.0,
              "Box classification parameter must be positive and finite");
      const Primal representative = multiplier_to_primal(multiplier);
      std::vector<BoxActivity> activities;
      activities.reserve(bounds_.layout()->dimension(0));
      for (std::size_t index = 0; index < bounds_.layout()->dimension(0); ++index)
        {
          const double value = Backend::value(primal.block(0), index);
          const double multiplier_value =
            Backend::value(representative.block(0), index);
          const double lower =
            Backend::value(bounds_.lower().block(0), index);
          const double upper =
            Backend::value(bounds_.upper().block(0), index);
          const double upper_expression =
            multiplier_value + c * (value - upper);
          const double lower_expression =
            multiplier_value + c * (value - lower);
          if (upper_expression > 0.0)
            activities.push_back(BoxActivity::upper);
          else if (lower_expression < 0.0)
            activities.push_back(BoxActivity::lower);
          else
            activities.push_back(BoxActivity::inactive);
        }
      return Selection(bounds_.layout(), std::move(activities));
    }

  private:
    Bounds          bounds_;
    Representation  representation_;
  };

  using BoxComplementarity = BoxComplementarityT<DenseBackend>;
} // namespace nmopt::contract
