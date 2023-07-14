/**
 * Random utility functions.
 */
#ifndef BTLLIB_UTIL_HPP
#define BTLLIB_UTIL_HPP

#include "btllib/cstring.hpp"
#include "btllib/status.hpp"

#include <cmath>
#include <condition_variable>
#include <mutex>
#include <string>
#include <vector>

namespace btllib
{

  static constexpr double PHRED_OFFSET = 33.0;

  /**
   * Split a string into component substrings with a delimiter.
   *
   * @param s String to split.
   * @param delim Delimiter to split with.
   *
   * @return Vector of substrings delimited by `delim`, excluding delimiters
   * themselves.
   */
  std::vector<std::string>
  split(const std::string &s, const std::string &delim);

  /**
   * Join a vector of strings into a single string with a delimiter.
   *
   * @param s Vector of strings to join.
   * @param delim Delimiter to join the strings with.
   *
   * @return String with all the components joined.
   */
  std::string
  join(const std::vector<std::string> &s, const std::string &delim);

  /**
   * Trim whitespace on the left side of the given string.
   *
   * @param s String to trim, edited in-place.
   *
   */
  void
  ltrim(std::string &s);
  void
  ltrim(btllib::CString &s);

  /**
   * Trim whitespace on the right side of the given string.
   *
   * @param s String to trim, edited in-place.
   *
   */
  void
  rtrim(std::string &s);
  void
  rtrim(btllib::CString &s);

  /**
   * Trim whitespace on the left and right side of the given string.
   *
   * @param s String to trim, edited in-place.
   *
   */
  void
  trim(std::string &s);
  void
  trim(btllib::CString &s);

  /**
   * Check whether the given string starts with a prefix.
   *
   * @param s String to check.
   * @param prefix Prefix to check for.
   *
   */
  bool
  startswith(std::string s, std::string prefix);

  /**
   * Check whether the given string ends with a suffix.
   *
   * @param s String to check.
   * @param suffix Suffix to check for.
   *
   */
  bool
  endswith(std::string s, std::string suffix);

  /**
   * Equivalent to the GNU implementation of basename,
   * but returns a string copy of the result.
   *
   * @param path The path to get basename from.
   *
   * @return The basename of the path.
   */
  std::string
  get_basename(const std::string &path);

  /**
   * Equivalent to the GNU implementation of dirname,
   * but returns a string copy of the result.
   *
   * @param path The path to get dirname from.
   *
   * @return The dirname of the path.
   */
  std::string
  get_dirname(const std::string &path);

  /**
   * Calculate the average phred score of a string,
   * depending on the start position and length.
   *
   * @param qual The quality string to calculate the average from.
   * @param start_pos The start position of the substring. Defaults to 0.
   * @param len The length of the substring. Defaults to 0. If 0, the whole string
   * is used.
   *
   * @return The average phred score of the substring.
   */
  double
  calc_phred_avg(const std::string &qual, size_t start_pos = 0, size_t len = 0);

  /**
   * Range minimum query class.
   * Finds the index of the minimum value in a range of any iterable that can be indexed via [x]
   * where x is the index corresponding to the element in the iterable.
   * @tparam T The type of the iterable.
   */
  template <class T>
  class RangeMinimumQuery
  {
  public:
    /**
     * Constructor for RangeMinimumQuery.
     * @param array The iterable to be queried.
     * @param size_of_array The size of the iterable.
     */
    RangeMinimumQuery(const T &array, size_t size_of_array);
    /**
     * Query the index of the minimum value in the range [start, end].
     * @param start The start index of the range.
     * @param end The end index of the range.
     * @return The index of the minimum value in the range.
     * @pre start <= end
     * @pre end < size_of_array
     * @pre start >= 0
     *
     * @post The return value is in the range [start, end]
     * @post The return value is the index of the minimum value in the range.
     */
    size_t query(size_t start, size_t end);

  private:
    //* The lookup table for the query.
    std::vector<std::vector<size_t>> lookup_table;
  };

  // adapted from https://www.geeksforgeeks.org/range-minimum-query-for-static-array/
  template <class T>
  RangeMinimumQuery<T>::RangeMinimumQuery(const T &array, size_t size_of_array)
  {
    size_t log_size = (size_t)(log2(size_of_array)) + 1; // NOLINT(bugprone-narrowing-conversions,cppcoreguidelines-narrowing-conversions)
    lookup_table.resize(log_size, std::vector<size_t>(size_of_array));
    for (size_t i = 0; i < size_of_array; i++)
    {
      lookup_table[0][i] = i;
    }
    for (size_t j = 1; j < log_size; j++)
    {
      for (size_t i = 0; i + (1 << j) - 1 < size_of_array; i++)
      {
        if (array[lookup_table[j - 1][i]] < array[lookup_table[j - 1][i + (1 << (j - 1))]])
        {
          lookup_table[j][i] = lookup_table[j - 1][i];
        }
        else
        {
          lookup_table[j][i] = lookup_table[j - 1][i + (1 << (j - 1))];
        }
      }
    }
  }

  template <class T>
  size_t RangeMinimumQuery<T>::query(size_t start, size_t end)
  {
    if (start > end)
    {
      log_error("RangeMinimumQuery::query: start > end");
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
    if (end >= lookup_table[0].size())
    {
      log_error("RangeMinimumQuery::query: end >= lookup_table[0].size()");
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
    if (start >= lookup_table[0].size())
    {
      log_error("RangeMinimumQuery::query: start >= lookup_table[0].size()");
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }

    auto j = (size_t)(log2(end - start + 1)); // NOLINT(bugprone-narrowing-conversions,cppcoreguidelines-narrowing-conversions)
    if (lookup_table[j][start] <= lookup_table[j][end - (1 << j) + 1])
    {
      return lookup_table[j][start];
    }
    return lookup_table[j][end - (1 << j) + 1];
  }

  // This exists in C++20, but we don't support that yet
  /// @cond HIDDEN_SYMBOLS
  class Barrier
  {

  public:
    Barrier(const unsigned count)
        : counter_default(count)
    {
    }

    void wait();

  private:
    std::mutex m;
    std::condition_variable cv;
    unsigned counter{0};
    unsigned counter_default;
    unsigned waiting{0};
  };
  /// @endcond

} // namespace btllib

#endif