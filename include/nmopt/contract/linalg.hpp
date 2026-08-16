#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::contract
{
  class ContractError : public std::runtime_error
  {
  public:
    explicit ContractError(const std::string &message)
      : std::runtime_error(message)
    {}
  };

  inline void
  require(const bool condition, const std::string &message)
  {
    if (!condition)
      throw ContractError(message);
  }

  class DenseVector
  {
  public:
    DenseVector() = default;

    explicit DenseVector(const std::size_t size, const double value = 0.0)
      : values_(size, value)
    {}

    DenseVector(std::initializer_list<double> values)
      : values_(values)
    {}

    explicit DenseVector(std::vector<double> values)
      : values_(std::move(values))
    {}

    std::size_t
    size() const
    {
      return values_.size();
    }

    double &
    operator[](const std::size_t index)
    {
      return values_.at(index);
    }

    const double &
    operator[](const std::size_t index) const
    {
      return values_.at(index);
    }

    void
    add_scaled(const double factor, const DenseVector &other)
    {
      require(size() == other.size(), "DenseVector sizes do not match");
      for (std::size_t i = 0; i < size(); ++i)
        values_[i] += factor * other[i];
    }

    void
    scale(const double factor)
    {
      for (auto &value : values_)
        value *= factor;
    }

    const std::vector<double> &
    values() const
    {
      return values_;
    }

  private:
    std::vector<double> values_;
  };

  inline double
  dot(const DenseVector &left, const DenseVector &right)
  {
    require(left.size() == right.size(), "DenseVector sizes do not match");

    double value = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i)
      value += left[i] * right[i];
    return value;
  }

  class DenseMatrix
  {
  public:
    DenseMatrix() = default;

    DenseMatrix(const std::size_t          rows,
                const std::size_t          columns,
                std::vector<double> entries)
      : rows_(rows)
      , columns_(columns)
      , entries_(std::move(entries))
    {
      require(rows_ > 0 && columns_ > 0, "DenseMatrix must be non-empty");
      require(entries_.size() == rows_ * columns_,
              "DenseMatrix entry count does not match its shape");
    }

    std::size_t
    rows() const
    {
      return rows_;
    }

    std::size_t
    columns() const
    {
      return columns_;
    }

    double &
    operator()(const std::size_t row, const std::size_t column)
    {
      return entries_.at(row * columns_ + column);
    }

    const double &
    operator()(const std::size_t row, const std::size_t column) const
    {
      return entries_.at(row * columns_ + column);
    }

    DenseVector
    vmult(const DenseVector &input) const
    {
      require(input.size() == columns_,
              "DenseMatrix input does not match column count");

      DenseVector output(rows_);
      for (std::size_t row = 0; row < rows_; ++row)
        for (std::size_t column = 0; column < columns_; ++column)
          output[row] += (*this)(row, column) * input[column];
      return output;
    }

    DenseVector
    transpose_vmult(const DenseVector &input) const
    {
      require(input.size() == rows_,
              "DenseMatrix transpose input does not match row count");

      DenseVector output(columns_);
      for (std::size_t row = 0; row < rows_; ++row)
        for (std::size_t column = 0; column < columns_; ++column)
          output[column] += (*this)(row, column) * input[row];
      return output;
    }

    DenseMatrix
    transpose() const
    {
      std::vector<double> entries(columns_ * rows_);
      for (std::size_t row = 0; row < rows_; ++row)
        for (std::size_t column = 0; column < columns_; ++column)
          entries[column * rows_ + row] = (*this)(row, column);
      return DenseMatrix(columns_, rows_, std::move(entries));
    }

    DenseVector
    solve(DenseVector right_hand_side) const
    {
      require(rows_ == columns_, "Only square DenseMatrix systems can be solved");
      require(right_hand_side.size() == rows_,
              "DenseMatrix right-hand side has the wrong size");

      DenseMatrix work = *this;
      constexpr double pivot_tolerance = 1e-14;

      for (std::size_t pivot_column = 0; pivot_column < columns_;
           ++pivot_column)
        {
          std::size_t pivot_row = pivot_column;
          for (std::size_t candidate = pivot_column + 1; candidate < rows_;
               ++candidate)
            if (std::abs(work(candidate, pivot_column)) >
                std::abs(work(pivot_row, pivot_column)))
              pivot_row = candidate;

          require(std::abs(work(pivot_row, pivot_column)) > pivot_tolerance,
                  "DenseMatrix is singular to the configured pivot tolerance");

          if (pivot_row != pivot_column)
            {
              for (std::size_t column = pivot_column; column < columns_;
                   ++column)
                std::swap(work(pivot_row, column), work(pivot_column, column));
              std::swap(right_hand_side[pivot_row],
                        right_hand_side[pivot_column]);
            }

          for (std::size_t row = pivot_column + 1; row < rows_; ++row)
            {
              const double factor =
                work(row, pivot_column) / work(pivot_column, pivot_column);
              work(row, pivot_column) = 0.0;
              for (std::size_t column = pivot_column + 1; column < columns_;
                   ++column)
                work(row, column) -= factor * work(pivot_column, column);
              right_hand_side[row] -= factor * right_hand_side[pivot_column];
            }
        }

      DenseVector solution(columns_);
      for (std::size_t backward = columns_; backward > 0; --backward)
        {
          const std::size_t row = backward - 1;
          double            value = right_hand_side[row];
          for (std::size_t column = row + 1; column < columns_; ++column)
            value -= work(row, column) * solution[column];
          solution[row] = value / work(row, row);
        }
      return solution;
    }

  private:
    std::size_t         rows_ = 0;
    std::size_t         columns_ = 0;
    std::vector<double> entries_;
  };

  // The v0 reference backend. Production backends supply the same small
  // vector capability surface without exposing their implementation here.
  struct DenseBackend
  {
    using Vector = DenseVector;

    static Vector
    zeros(const std::size_t size)
    {
      return Vector(size);
    }

    static std::size_t
    size(const Vector &vector)
    {
      return vector.size();
    }

    static double
    value(const Vector &vector, const std::size_t index)
    {
      return vector[index];
    }

    static void
    set_value(Vector &vector, const std::size_t index, const double value)
    {
      vector[index] = value;
    }

    static double
    dot(const Vector &left, const Vector &right)
    {
      return nmopt::contract::dot(left, right);
    }

    static void
    add_scaled(Vector &target, const double factor, const Vector &source)
    {
      target.add_scaled(factor, source);
    }

    static void
    scale(Vector &target, const double factor)
    {
      target.scale(factor);
    }
  };
} // namespace nmopt::contract
