#ifndef BTL_CHECK_HPP
#define BTL_CHECK_HPP

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

namespace btl {

template<typename Arg, typename... Args>
inline void
check_error(bool condition, Arg&& arg, Args&&... args)
{
  if (condition) {
    std::cerr << "[ERROR] " << std::forward<Arg>(arg);
    ((std::cerr << std::forward<Args>(args)), ...);
    std::cerr << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

template<typename Arg, typename... Args>
inline void
check_warning(bool condition, Arg&& arg, Args&&... args)
{
  if (condition) {
    std::cerr << "[WARNING] " << std::forward<Arg>(arg);
    ((std::cerr << std::forward<Args>(args)), ...);
    std::cerr << std::endl;
  }
}

inline void
check_stream(const std::ios& stream, const char* path)
{
  check_error(
    !stream.good(), "'", path, "' stream error: ", std::strerror(errno));
}

} // namespace btl

#endif