#ifndef BTL_STATUS_HPP
#define BTL_STATUS_HPP

#include <cstdlib>
#include <cstring>
#include <chrono>
#include <iostream>
#include <string>

namespace btl {

inline void
log_info(const std::string& msg) {
    std::cerr << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now())
        << "[INFO] " << msg << std::endl;
}

inline void
log_warning(const std::string& msg) {
    std::cerr << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now())
        << "[WARNING] " << msg << std::endl;
}

inline void
log_error(const std::string& msg) {
    std::cerr << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now())
        << "[ERROR] " << msg << std::endl;
}

inline void
check_error(bool condition, const std::string& msg)
{
    if (condition) {
        log_error(msg);
        std::exit(EXIT_FAILURE);
    }
}

inline void
check_warning(bool condition, const std::string& msg)
{
    if (condition) {
        log_warning(msg);
    }
}

inline void
check_stream(const std::ios& stream, const std::string& name)
{
  check_error(
    !stream.good(), "'" + name + "' stream error: " + std::strerror(errno));
}

} // namespace btl

#endif