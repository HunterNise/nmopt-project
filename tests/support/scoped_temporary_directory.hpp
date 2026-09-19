#pragma once

#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>

namespace nmopt::test_support
{
  class ScopedTemporaryDirectory final
  {
  public:
    explicit ScopedTemporaryDirectory(std::string prefix)
      : path_(create_unique_directory(std::move(prefix)))
    {}

    ScopedTemporaryDirectory(const ScopedTemporaryDirectory &) = delete;
    ScopedTemporaryDirectory &
    operator=(const ScopedTemporaryDirectory &) = delete;
    ScopedTemporaryDirectory(ScopedTemporaryDirectory &&) = delete;
    ScopedTemporaryDirectory &
    operator=(ScopedTemporaryDirectory &&) = delete;

    ~ScopedTemporaryDirectory() noexcept
    {
      std::error_code error;
      std::filesystem::remove_all(path_, error);
    }

    const std::filesystem::path &
    path() const noexcept
    {
      return path_;
    }

  private:
    static std::filesystem::path
    create_unique_directory(std::string prefix)
    {
      const auto parent = std::filesystem::temp_directory_path();
      const auto timestamp = std::chrono::duration_cast<std::chrono::nanoseconds>(
                               std::chrono::steady_clock::now().time_since_epoch())
                               .count();
      static std::atomic<std::uint64_t> sequence{0};
      const auto instance = sequence.fetch_add(1, std::memory_order_relaxed);

      for (std::uint64_t attempt = 0;; ++attempt)
        {
          const auto candidate =
            parent / (prefix + "-" + std::to_string(timestamp) + "-" +
                      std::to_string(instance) + "-" +
                      std::to_string(attempt));
          std::error_code error;
          if (std::filesystem::create_directory(candidate, error))
            return candidate;
          if (error && error != std::errc::file_exists)
            throw std::system_error(
              error, "could not create unique temporary test directory");
        }
    }

    std::filesystem::path path_;
  };
} // namespace nmopt::test_support
