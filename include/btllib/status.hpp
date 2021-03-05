/**
 * Functions for logging and error checking.
 */
#ifndef BTLLIB_STATUS_HPP
#define BTLLIB_STATUS_HPP

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <string>

namespace btllib {

inline std::string
get_time();

/**
 * Log info level events.
 *
 * @param msg Message to print.
 */
inline void
log_info(const std::string& msg);

/**
 * Log warning level events.
 *
 * @param msg Message to print.
 */
inline void
log_warning(const std::string& msg);

/**
 * Log error level events.
 *
 * @param msg Message to print.
 */
inline void
log_error(const std::string& msg);

/**
 * Conditionally log info level events.
 *
 * @param condition If this is true, the message is printed.
 * @param msg Message to print.
 */
inline void
check_info(bool condition, const std::string& msg);

/**
 * Conditionally log warning level events.
 *
 * @param condition If this is true, the message is printed.
 * @param msg Message to print.
 */
inline void
check_warning(bool condition, const std::string& msg);

/**
 * Conditionally log error level events. The program exits if the condition is
 * true.
 *
 * @param condition If this is true, the message is printed and the program
 * exits.
 * @param msg Message to print.
 */
inline void
check_error(bool condition, const std::string& msg);

/**
 * Check whether the stream is good. Program prints an error message and exits
 * if not.
 *
 * @param stream Stream to check goodness of.
 * @param name Name of the stream, e.g. filepath or stdin
 */
inline void
check_stream(const std::ios& stream, const std::string& name);

inline std::string
get_time()
{
  time_t now;
  time(&now);
  char buf[sizeof("2011-10-08T07:07:09Z")];
  strftime(buf, sizeof buf, "%F %T", localtime(&now));
  return std::string(buf);
}

inline void
log_info(const std::string& msg)
{
  std::cerr << ('[' + get_time() + "] [INFO] " + msg + '\n') << std::flush;
}

inline void
log_warning(const std::string& msg)
{
  std::cerr << ('[' + get_time() + "] [WARNING] " + msg + '\n') << std::flush;
}

inline void
log_error(const std::string& msg)
{
  std::cerr << ('[' + get_time() + "] [ERROR] " + msg + '\n') << std::flush;
}

inline void
check_info(bool condition, const std::string& msg)
{
  if (condition) {
    log_info(msg);
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
check_error(bool condition, const std::string& msg)
{
  if (condition) {
    log_error(msg);
    std::exit(EXIT_FAILURE);
  }
}

inline void
check_stream(const std::ios& stream, const std::string& name)
{
  check_error(!stream.good(),
              "'" + name + "' stream error: " + std::strerror(errno));
}

} // namespace btllib

#endif
