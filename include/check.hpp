#ifndef BTL_CHECK_H
#define BTL_CHECK_H

#include <cstdlib>
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
    std::cerr << '\n';
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
    std::cerr << '\n';
  }
}

} // namespace btl

#endif