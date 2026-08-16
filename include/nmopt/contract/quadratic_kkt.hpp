#pragma once

#include "nmopt/contract/layout.hpp"

#include <algorithm>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::contract
{
  enum class QuadraticKKTSymmetry
  {
    symmetric_indefinite,
    nonsymmetric
  };

  struct QuadraticKKTAssumptions
  {
    bool        rank_condition_declared = false;
    bool        kernel_positivity_declared = false;
    std::string rank_policy;
    std::string kernel_policy;
    bool        d_transpose_consistency_declared = false;
    bool        kkt_transpose_consistency_declared = false;
    std::string transpose_consistency_policy;
  };

  struct QuadraticKKTBlockPairing
  {
    std::string              id;
    std::vector<std::size_t> domain_blocks;
    std::vector<std::size_t> range_blocks;
    std::vector<std::string> pairing_ids;
  };

  template <typename Backend>
  struct QuadraticKKTLayoutT
  {
    LayoutPtr primal;
    LayoutPtr multiplier;
    LayoutPtr adjoint;
    LayoutPtr stationarity;
    LayoutPtr equality;
    QuadraticKKTBlockPairing primal_stationarity_pairing;
    QuadraticKKTBlockPairing multiplier_equality_pairing;

    QuadraticKKTLayoutT(LayoutPtr primal_layout,
                        LayoutPtr multiplier_layout,
                        LayoutPtr adjoint_layout,
                        LayoutPtr stationarity_layout,
                        LayoutPtr equality_layout,
                        QuadraticKKTBlockPairing primal_stationarity = {},
                        QuadraticKKTBlockPairing multiplier_equality = {})
      : primal(std::move(primal_layout))
      , multiplier(std::move(multiplier_layout))
      , adjoint(std::move(adjoint_layout))
      , stationarity(std::move(stationarity_layout))
      , equality(std::move(equality_layout))
      , primal_stationarity_pairing(std::move(primal_stationarity))
      , multiplier_equality_pairing(std::move(multiplier_equality))
    {
      require(static_cast<bool>(primal) && static_cast<bool>(multiplier) &&
                static_cast<bool>(adjoint) && static_cast<bool>(stationarity) &&
                static_cast<bool>(equality),
              "Quadratic KKT layouts must all be present");
      validate_pairing(primal_stationarity_pairing,
                       *primal,
                       *stationarity,
                       "primal/stationarity");
      validate_pairing(multiplier_equality_pairing,
                       *multiplier,
                       *equality,
                       "multiplier/equality");
    }

    bool
    has_complete_pairings() const
    {
      return !primal_stationarity_pairing.id.empty() &&
             !multiplier_equality_pairing.id.empty();
    }

  private:
    static void
    validate_pairing(const QuadraticKKTBlockPairing &pairing,
                     const BlockLayout &             domain,
                     const BlockLayout &             range,
                     const char *                    name)
    {
      const std::string prefix = "Quadratic KKT " + std::string(name) +
                                 " pairing ";
      require(!pairing.id.empty(), prefix + "needs an identifier");
      require(pairing.domain_blocks.size() == domain.n_blocks() &&
                pairing.range_blocks.size() == range.n_blocks() &&
                pairing.pairing_ids.size() == domain.n_blocks(),
              prefix + "needs a complete block map");

      std::vector<bool> domain_seen(domain.n_blocks(), false);
      std::vector<bool> range_seen(range.n_blocks(), false);
      for (std::size_t index = 0; index < pairing.domain_blocks.size(); ++index)
        {
          const std::size_t domain_block = pairing.domain_blocks[index];
          const std::size_t range_block = pairing.range_blocks[index];
          require(domain_block < domain.n_blocks() &&
                    range_block < range.n_blocks(),
                  prefix + "contains an out-of-range block");
          require(!domain_seen[domain_block] && !range_seen[range_block],
                  prefix + "must be one-to-one");
          require(domain.dimension(domain_block) == range.dimension(range_block),
                  prefix + "has incompatible block dimensions");
          require(!pairing.pairing_ids[index].empty(),
                  prefix + "needs an identifier for every block pair");
          require(std::find(pairing.pairing_ids.begin(),
                            pairing.pairing_ids.begin() + index,
                            pairing.pairing_ids[index]) ==
                    pairing.pairing_ids.begin() + index,
                  prefix + "pairing identifiers must be unique");
          domain_seen[domain_block] = true;
          range_seen[range_block] = true;
        }

      require(std::all_of(domain_seen.begin(), domain_seen.end(),
                          [](const bool value) { return value; }) &&
                std::all_of(range_seen.begin(), range_seen.end(),
                            [](const bool value) { return value; }),
              prefix + "must cover every domain and range block");
    }
  };

  using QuadraticKKTLayout = QuadraticKKTLayoutT<DenseBackend>;

  template <typename Backend>
  struct QuadraticKKTMultiplierConversionT
  {
    using Primal = PrimalBlockT<Backend>;
    using ConversionAction = std::function<Primal(const Primal &)>;

    std::string     description;
    ConversionAction to_adjoint;
    ConversionAction to_multiplier;
  };

  using QuadraticKKTMultiplierConversion =
    QuadraticKKTMultiplierConversionT<DenseBackend>;

  template <typename Backend>
  struct QuadraticKKTPointT
  {
    PrimalBlockT<Backend> primal;
    PrimalBlockT<Backend> multiplier;
  };

  template <typename Backend>
  struct QuadraticKKTSeedT
  {
    PrimalBlockT<Backend> stationarity;
    PrimalBlockT<Backend> equality;
  };

  template <typename Backend>
  struct QuadraticKKTResidualT
  {
    CovectorBlockT<Backend> stationarity;
    CovectorBlockT<Backend> equality;
  };

  template <typename Backend>
  struct QuadraticKKTTransposeResultT
  {
    CovectorBlockT<Backend> primal;
    CovectorBlockT<Backend> multiplier;
  };

  template <typename Backend>
  class EqualityConstrainedQuadraticKKTProductT final
  {
  public:
    using Layout = QuadraticKKTLayoutT<Backend>;
    using Primal = PrimalBlockT<Backend>;
    using Covector = CovectorBlockT<Backend>;
    using Point = QuadraticKKTPointT<Backend>;
    using Seed = QuadraticKKTSeedT<Backend>;
    using Residual = QuadraticKKTResidualT<Backend>;
    using TransposeResult = QuadraticKKTTransposeResultT<Backend>;
    using MultiplierConversion = QuadraticKKTMultiplierConversionT<Backend>;

    using QuadraticAction = std::function<Covector(const Primal &)>;
    using EqualityAction = std::function<Covector(const Primal &)>;
    using MultiplierAction = std::function<Covector(const Primal &)>;
    using TransposeAction = std::function<TransposeResult(const Seed &)>;

    EqualityConstrainedQuadraticKKTProductT(
      Layout                         layout,
      QuadraticAction                quadratic_action,
      EqualityAction                equality_action,
      MultiplierAction              multiplier_action,
      TransposeAction               transpose_action,
      Covector                      stationarity_rhs,
      Covector                      equality_rhs,
      MultiplierConversion           multiplier_conversion,
      QuadraticKKTAssumptions        assumptions,
      const QuadraticKKTSymmetry     symmetry)
      : layout_(std::move(layout))
      , quadratic_action_(std::move(quadratic_action))
      , equality_action_(std::move(equality_action))
      , multiplier_action_(std::move(multiplier_action))
      , transpose_action_(std::move(transpose_action))
      , stationarity_rhs_(std::move(stationarity_rhs))
      , equality_rhs_(std::move(equality_rhs))
      , multiplier_conversion_(std::move(multiplier_conversion))
      , assumptions_(std::move(assumptions))
      , symmetry_(symmetry)
    {
      require(static_cast<bool>(quadratic_action_),
              "Quadratic KKT product needs a Q action");
      require(static_cast<bool>(equality_action_),
              "Quadratic KKT product needs a D action");
      require(static_cast<bool>(multiplier_action_),
              "Quadratic KKT product needs a D transpose action");
      require(static_cast<bool>(transpose_action_),
              "Quadratic KKT product needs a KKT transpose action");
      require_layout(stationarity_rhs_, layout_.stationarity,
                     "Quadratic KKT stationarity right-hand side");
      require_layout(equality_rhs_, layout_.equality,
                     "Quadratic KKT equality right-hand side");
      require(!multiplier_conversion_.description.empty(),
              "Quadratic KKT multiplier conversion needs a description");
      require(static_cast<bool>(multiplier_conversion_.to_adjoint),
              "Quadratic KKT product needs multiplier-to-adjoint conversion");
      require(static_cast<bool>(multiplier_conversion_.to_multiplier),
              "Quadratic KKT product needs adjoint-to-multiplier conversion");
      require(assumptions_.rank_condition_declared,
              "Quadratic KKT product needs a declared rank condition");
      require(!assumptions_.rank_policy.empty(),
              "Quadratic KKT rank condition needs a policy");
      require(assumptions_.kernel_positivity_declared,
              "Quadratic KKT product needs a declared kernel-positivity condition");
      require(!assumptions_.kernel_policy.empty(),
              "Quadratic KKT kernel positivity needs a policy");
      if (symmetry_ == QuadraticKKTSymmetry::symmetric_indefinite)
        {
          require(layout_.has_complete_pairings(),
                  "Symmetric quadratic KKT product needs validated domain-range pairings");
          require(assumptions_.d_transpose_consistency_declared,
                  "Symmetric quadratic KKT product needs declared D-transpose consistency");
          require(assumptions_.kkt_transpose_consistency_declared,
                  "Symmetric quadratic KKT product needs declared KKT-transpose consistency");
          require(!assumptions_.transpose_consistency_policy.empty(),
                  "Symmetric quadratic KKT transpose consistency needs a policy");
        }
    }

    const Layout &
    layout() const
    {
      return layout_;
    }

    const QuadraticKKTAssumptions &
    assumptions() const
    {
      return assumptions_;
    }

    QuadraticKKTSymmetry
    symmetry() const
    {
      return symmetry_;
    }

    bool
    supports_minres() const
    {
      return symmetry_ == QuadraticKKTSymmetry::symmetric_indefinite &&
             layout_.has_complete_pairings() &&
             assumptions_.d_transpose_consistency_declared &&
             assumptions_.kkt_transpose_consistency_declared;
    }

    Covector
    apply_q(const Primal &primal) const
    {
      require_layout(primal, layout_.primal, "Quadratic KKT Q input");
      return checked_covector(quadratic_action_(primal),
                              layout_.stationarity,
                              "Quadratic KKT Q action");
    }

    Covector
    apply_d(const Primal &primal) const
    {
      require_layout(primal, layout_.primal, "Quadratic KKT D input");
      return checked_covector(equality_action_(primal),
                              layout_.equality,
                              "Quadratic KKT D action");
    }

    Covector
    apply_d_transpose(const Primal &multiplier) const
    {
      require_layout(multiplier,
                     layout_.multiplier,
                     "Quadratic KKT D transpose input");
      return checked_covector(multiplier_action_(multiplier),
                              layout_.stationarity,
                              "Quadratic KKT D transpose action");
    }

    Residual
    apply_kkt(const Point &point) const
    {
      require_layout(point.primal, layout_.primal, "Quadratic KKT primal input");
      require_layout(point.multiplier,
                     layout_.multiplier,
                     "Quadratic KKT multiplier input");

      Covector stationarity = apply_q(point.primal);
      add_covector(stationarity, apply_d_transpose(point.multiplier), 1.0,
                   "Quadratic KKT stationarity action");
      return {std::move(stationarity), apply_d(point.primal)};
    }

    Residual
    residual(const Point &point) const
    {
      Residual value = apply_kkt(point);
      add_covector(value.stationarity,
                   stationarity_rhs_,
                   -1.0,
                   "Quadratic KKT stationarity residual");
      add_covector(value.equality,
                   equality_rhs_,
                   -1.0,
                   "Quadratic KKT equality residual");
      return value;
    }

    TransposeResult
    apply_kkt_transpose(const Seed &seed) const
    {
      require_layout(seed.stationarity,
                     layout_.stationarity,
                     "Quadratic KKT transpose stationarity seed");
      require_layout(seed.equality,
                     layout_.equality,
                     "Quadratic KKT transpose equality seed");
      TransposeResult result = transpose_action_(seed);
      require_layout(result.primal,
                     layout_.primal,
                     "Quadratic KKT transpose primal result");
      require_layout(result.multiplier,
                     layout_.multiplier,
                     "Quadratic KKT transpose multiplier result");
      return result;
    }

    Primal
    multiplier_to_adjoint(const Primal &multiplier) const
    {
      require_layout(multiplier,
                     layout_.multiplier,
                     "Quadratic KKT multiplier conversion input");
      Primal result = multiplier_conversion_.to_adjoint(multiplier);
      require_layout(result,
                     layout_.adjoint,
                     "Quadratic KKT multiplier-to-adjoint result");
      return result;
    }

    Primal
    adjoint_to_multiplier(const Primal &adjoint) const
    {
      require_layout(adjoint,
                     layout_.adjoint,
                     "Quadratic KKT adjoint conversion input");
      Primal result = multiplier_conversion_.to_multiplier(adjoint);
      require_layout(result,
                     layout_.multiplier,
                     "Quadratic KKT adjoint-to-multiplier result");
      return result;
    }

  private:
    template <typename Values>
    static void
    require_layout(const Values &value,
                   const LayoutPtr &expected,
                   const std::string &operation)
    {
      require(value.layout()->compatible_with(*expected),
              operation + " has an incompatible layout");
    }

    static Covector
    checked_covector(Covector value,
                     const LayoutPtr &expected,
                     const std::string &operation)
    {
      require_layout(value, expected, operation + " returned");
      return value;
    }

    static void
    add_covector(Covector &       target,
                 const Covector & source,
                 const double     factor,
                 const std::string &operation)
    {
      require_layout(source, target.layout(), operation + " source");
      for (std::size_t block = 0; block < target.n_blocks(); ++block)
        target.add_scaled_block(block, factor, source.block(block));
    }

    Layout                       layout_;
    QuadraticAction              quadratic_action_;
    EqualityAction               equality_action_;
    MultiplierAction             multiplier_action_;
    TransposeAction              transpose_action_;
    Covector                     stationarity_rhs_;
    Covector                     equality_rhs_;
    MultiplierConversion         multiplier_conversion_;
    QuadraticKKTAssumptions      assumptions_;
    QuadraticKKTSymmetry         symmetry_;
  };

  using EqualityConstrainedQuadraticKKTProduct =
    EqualityConstrainedQuadraticKKTProductT<DenseBackend>;
} // namespace nmopt::contract
