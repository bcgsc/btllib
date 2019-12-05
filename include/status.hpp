#ifndef BTL_STATUS_HPP
#define BTL_STATUS_HPP

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

namespace btl {
  
template<typename Arg, typename... Args>
inline void
raise_error(Arg&& arg, Args&&... args) {
    std::cerr << "[ERROR] " << std::forward<Arg>(arg);
    ((std::cerr << std::forward<Args>(args)), ...);
    std::cerr << std::endl;
    std::exit(EXIT_FAILURE);
}

template<typename Arg, typename... Args>
inline void
raise_warning() {
    std::cerr << "[WARNING] " << std::forward<Arg>(arg);
    ((std::cerr << std::forward<Args>(args)), ...);
    std::cerr << std::endl;
}

template<typename Arg, typename... Args>
inline void
check_error(bool condition, Arg&& arg, Args&&... args)
{
    if (condition) {
        raise_error(args, args);
    }
}

template<typename Arg, typename... Args>
inline void
check_warning(bool condition, Arg&& arg, Args&&... args)
{
    if (condition) {
        raise_warning(arg, args);
    }
}

inline void
check_stream(const std::ios& stream, const char* name)
{
  check_error(
    !stream.good(), "'", name, "' stream error: ", std::strerror(errno));
}

} // namespace btl

#endif