#pragma once

#include "nmopt/contract/linalg.hpp"

#include <memory>
#include <string>
#include <unordered_set>
#include <utility>

namespace nmopt::contract
{
  struct SpaceId
  {
    std::string value;

    bool
    operator==(const SpaceId &other) const
    {
      return value == other.value;
    }
  };

  class BlockLayout
  {
  public:
    BlockLayout(std::string              label,
                std::vector<SpaceId>     spaces,
                std::vector<std::size_t> dimensions)
      : label_(std::move(label))
      , spaces_(std::move(spaces))
      , dimensions_(std::move(dimensions))
    {
      require(!label_.empty(), "BlockLayout label must not be empty");
      require(!spaces_.empty(), "BlockLayout must contain at least one block");
      require(spaces_.size() == dimensions_.size(),
              "BlockLayout spaces and dimensions do not match");

      std::unordered_set<std::string> identifiers;
      for (std::size_t block = 0; block < spaces_.size(); ++block)
        {
          require(!spaces_[block].value.empty(),
                  "BlockLayout space identifier must not be empty");
          require(dimensions_[block] > 0,
                  "BlockLayout block dimension must be positive");
          require(identifiers.insert(spaces_[block].value).second,
                  "BlockLayout space identifiers must be unique");
        }
    }

    const std::string &
    label() const
    {
      return label_;
    }

    std::size_t
    n_blocks() const
    {
      return spaces_.size();
    }

    const SpaceId &
    space(const std::size_t block) const
    {
      return spaces_.at(block);
    }

    std::size_t
    dimension(const std::size_t block) const
    {
      return dimensions_.at(block);
    }

    bool
    compatible_with(const BlockLayout &other) const
    {
      return spaces_ == other.spaces_ && dimensions_ == other.dimensions_;
    }

    std::shared_ptr<const BlockLayout>
    single_block(const std::size_t block, std::string label) const
    {
      return std::make_shared<const BlockLayout>(
        std::move(label), std::vector<SpaceId>{space(block)},
        std::vector<std::size_t>{dimension(block)});
    }

  private:
    std::string              label_;
    std::vector<SpaceId>     spaces_;
    std::vector<std::size_t> dimensions_;
  };

  using LayoutPtr = std::shared_ptr<const BlockLayout>;

  template <typename Backend>
  class BlockValuesT
  {
  public:
    const LayoutPtr &
    layout() const
    {
      return layout_;
    }

    std::size_t
    n_blocks() const
    {
      return blocks_.size();
    }

    const typename Backend::Vector &
    block(const std::size_t index) const
    {
      return blocks_.at(index);
    }

    void
    add_scaled_block(const std::size_t               index,
                     const double                    factor,
                     const typename Backend::Vector &source)
    {
      require(Backend::size(source) == layout_->dimension(index),
              "BlockValues update vector dimension does not match layout");
      Backend::add_scaled(blocks_.at(index), factor, source);
    }

    void
    scale_block(const std::size_t index, const double factor)
    {
      Backend::scale(blocks_.at(index), factor);
    }

  protected:
    BlockValuesT(LayoutPtr                              layout,
                 std::vector<typename Backend::Vector> blocks)
      : layout_(std::move(layout))
      , blocks_(std::move(blocks))
    {
      require(static_cast<bool>(layout_), "BlockValues needs a layout");
      require(blocks_.size() == layout_->n_blocks(),
              "BlockValues block count does not match layout");
      for (std::size_t block = 0; block < blocks_.size(); ++block)
        require(Backend::size(blocks_[block]) == layout_->dimension(block),
                "BlockValues vector dimension does not match layout");
    }

    static std::vector<typename Backend::Vector>
    zero_blocks(const LayoutPtr &layout)
    {
      require(static_cast<bool>(layout), "BlockValues needs a layout");
      std::vector<typename Backend::Vector> blocks;
      blocks.reserve(layout->n_blocks());
      for (std::size_t block = 0; block < layout->n_blocks(); ++block)
        blocks.emplace_back(Backend::zeros(layout->dimension(block)));
      return blocks;
    }

  private:
    LayoutPtr                layout_;
    std::vector<typename Backend::Vector> blocks_;
  };

  template <typename Backend>
  class PrimalBlockT final : public BlockValuesT<Backend>
  {
  public:
    PrimalBlockT(LayoutPtr layout, std::vector<typename Backend::Vector> blocks)
      : BlockValuesT<Backend>(std::move(layout), std::move(blocks))
    {}

    static PrimalBlockT
    zeros(const LayoutPtr &layout)
    {
      return PrimalBlockT(layout, BlockValuesT<Backend>::zero_blocks(layout));
    }
  };

  template <typename Backend>
  class CovectorBlockT final : public BlockValuesT<Backend>
  {
  public:
    CovectorBlockT(LayoutPtr layout, std::vector<typename Backend::Vector> blocks)
      : BlockValuesT<Backend>(std::move(layout), std::move(blocks))
    {}

    static CovectorBlockT
    zeros(const LayoutPtr &layout)
    {
      return CovectorBlockT(layout,
                            BlockValuesT<Backend>::zero_blocks(layout));
    }
  };

  using BlockValues = BlockValuesT<DenseBackend>;
  using PrimalBlock = PrimalBlockT<DenseBackend>;
  using CovectorBlock = CovectorBlockT<DenseBackend>;

  template <typename Backend>
  inline void
  require_compatible(const BlockValuesT<Backend> &left,
                     const BlockValuesT<Backend> &right,
                     const std::string &what)
  {
    require(left.layout()->compatible_with(*right.layout()), what);
  }

  template <typename Backend>
  inline double
  pair(const CovectorBlockT<Backend> &covector,
       const PrimalBlockT<Backend>   &primal)
  {
    require_compatible(covector, primal,
                       "Cannot pair incompatible primal and covector blocks");

    double value = 0.0;
    for (std::size_t block = 0; block < covector.n_blocks(); ++block)
      value += Backend::dot(covector.block(block), primal.block(block));
    return value;
  }

  template <typename Backend>
  inline CovectorBlockT<Backend>
  subtract(const CovectorBlockT<Backend> &left,
           const CovectorBlockT<Backend> &right)
  {
    require_compatible(left, right, "Cannot subtract incompatible covectors");
    CovectorBlockT<Backend> result = left;
    for (std::size_t block = 0; block < result.n_blocks(); ++block)
      result.add_scaled_block(block, -1.0, right.block(block));
    return result;
  }

  template <typename Backend>
  inline PrimalBlockT<Backend>
  extract_primal_block(const PrimalBlockT<Backend> &source,
                       const std::size_t            block,
                       const std::string           &label)
  {
    return PrimalBlockT<Backend>(
      source.layout()->single_block(block, label),
      std::vector<typename Backend::Vector>{source.block(block)});
  }

  template <typename Backend>
  inline CovectorBlockT<Backend>
  extract_covector_block(const CovectorBlockT<Backend> &source,
                         const std::size_t              block,
                         const std::string             &label)
  {
    return CovectorBlockT<Backend>(
      source.layout()->single_block(block, label),
      std::vector<typename Backend::Vector>{source.block(block)});
  }
} // namespace nmopt::contract
