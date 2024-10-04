#include "btllib/status.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

#include <sys/stat.h>

namespace btllib {

std::string
get_time()
{
  time_t now;
  const auto timeret = time(&now);
  if (timeret == (time_t)(-1)) {
    std::cerr << "btllib: time() failed." << std::endl;
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
  char buf[sizeof("2011-10-08T07:07:09Z")];
  std::tm tm_result = {};
  localtime_r(&now, &tm_result);
  const auto ret = std::strftime(buf, sizeof(buf), "%F %T", &tm_result);
  if (ret < sizeof(buf) - 2) {
    std::cerr << "btllib: strftime failed." << std::endl;
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
  return std::string(buf);
}

void
log_info(const std::string& msg)
{
  std::string info_msg = '[' + get_time() + "]" + PRINT_COLOR_INFO + "[INFO] " +
                         PRINT_COLOR_END + msg + '\n';
  std::cerr << info_msg << std::flush;
}

void
log_warning(const std::string& msg)
{
  std::string warning_msg = '[' + get_time() + "]" + PRINT_COLOR_WARNING +
                            "[WARNING] " + PRINT_COLOR_END + msg + '\n';
  std::cerr << warning_msg << std::flush;
}

void
log_error(const std::string& msg)
{
  std::string error_msg = '[' + get_time() + "]" + PRINT_COLOR_ERROR +
                          "[ERROR] " + PRINT_COLOR_END + msg + '\n';
  std::cerr << error_msg << std::flush;
}

void
check_info(bool condition, const std::string& msg)
{
  if (condition) {
    log_info(msg);
  }
}

void
check_warning(bool condition, const std::string& msg)
{
  if (condition) {
    log_warning(msg);
  }
}

void
check_error(bool condition, const std::string& msg)
{
  if (condition) {
    log_error(msg);
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
}

std::string
get_strerror()
{
  static const size_t buflen = 1024;
  char buf[buflen];
// POSIX and GNU implementation of strerror_r differ, even in function signature
// and so we need to check which one is used
#if __APPLE__ ||                                                               \
  ((_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && !_GNU_SOURCE)
  strerror_r(errno, buf, buflen);
  return buf;
#else
  return strerror_r(errno, buf, buflen);
#endif
}

void
check_stream(const std::ios& stream, const std::string& name)
{
  if (!stream.good()) {
    log_error("'" + name + "' stream error: " + get_strerror());
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
}

void
check_file_accessibility(const std::string& filepath)
{
  struct stat buffer
  {};
  const auto ret = stat(filepath.c_str(), &buffer);
  btllib::check_error(ret != 0, get_strerror() + ": " + filepath);
}

} // namespace btllib