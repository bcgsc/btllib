#ifndef BTLLIB_DATA_SAVELOAD_HPP
#define BTLLIB_DATA_SAVELOAD_HPP

#include "status.hpp"
#include "util.hpp"

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <csignal>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include <dlfcn.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

namespace btllib {

using open_t = int (*)(const char* __file, int __oflag, ...);
using creat_t = int (*)(const char* __file, mode_t __mode);
using fopen_t = FILE* (*)(const char* __restrict __filename,
                          const char* __restrict __modes);

static open_t orig_open = nullptr;
static open_t orig_open64 = nullptr;
static creat_t orig_creat = nullptr;
static fopen_t orig_fopen = nullptr;
static fopen_t orig_fopen64 = nullptr;

/** SIGCHLD handler. Reap child processes and report an error if any
 * fail. */
static void
sigchld_handler(const int sig)
{
  assert(sig == SIGCHLD);
  (void)sig;

  pid_t pid;
  int status;
  while ((pid = waitpid(-1, &status, WNOHANG)) > 0) {
    if (status != 0) {
      if (WIFEXITED(status)) { // NOLINT
        std::cerr << "PID " << pid << " exited with status "
                  << WEXITSTATUS(status) << std::endl; // NOLINT
      } else if (WIFSIGNALED(status)) {                // NOLINT
        std::cerr << "PID " << pid << " killed by signal "
                  << WTERMSIG(status) // NOLINT
                  << std::endl;
      } else {
        std::cerr << "PID " << pid << " exited with code " << status
                  << std::endl;
      }
      std::exit(EXIT_FAILURE);
    }
  }
  if (pid == -1 && errno != ECHILD) {
    std::perror("waitpid");
    std::exit(EXIT_FAILURE);
  }
}

static bool
data_saveload_init()
{
  struct sigaction action; // NOLINT
  action.sa_handler = sigchld_handler;
  sigemptyset(&action.sa_mask);
  action.sa_flags = SA_RESTART;
  sigaction(SIGCHLD, &action, nullptr);

  const char* err;

  dlerror();
  orig_open = (open_t)dlsym(RTLD_NEXT, "open");
  err = dlerror();
  check_error(orig_open == nullptr && err != nullptr,
              "dlsym open: " + (err != nullptr ? std::string(err) : ""));

  dlerror();
  orig_open64 = (open_t)dlsym(RTLD_NEXT, "open64");
  err = dlerror();
  check_error(orig_open64 == nullptr && err != nullptr,
              "dlsym open64: " + (err != nullptr ? std::string(err) : ""));

  dlerror();
  orig_creat = (creat_t)dlsym(RTLD_NEXT, "creat");
  err = dlerror();
  check_error(orig_creat == nullptr && err != nullptr,
              "dlsym creat: " + (err != nullptr ? std::string(err) : ""));

  dlerror();
  orig_fopen = (fopen_t)dlsym(RTLD_NEXT, "fopen");
  err = dlerror();
  check_error(orig_fopen == nullptr && err != nullptr,
              "dlsym fopen: " + (err != nullptr ? std::string(err) : ""));

  dlerror();
  orig_fopen64 = (fopen_t)dlsym(RTLD_NEXT, "fopen64");
  err = dlerror();
  check_error(orig_fopen64 == nullptr && err != nullptr,
              "dlsym fopen64: " + (err != nullptr ? std::string(err) : ""));

  return true;
}

static const bool data_saveload_initialized = data_saveload_init();

static inline std::string
get_saveload_cmd(const std::string& path, const bool save)
{
  struct Datatype
  {
    std::vector<std::string> prefixes;
    std::vector<std::string> suffixes;
    std::vector<std::string> cmds_check_existence;
    std::vector<std::string> load_cmds;
    std::vector<std::string> save_cmds;
  };

  static const Datatype DATATYPES[]{
    { { "http://", "https://", "ftp://" },
      {},
      { "wget --version" },
      { "wget -O-" },
      { "" } },
    { {}, { ".url" }, { "wget --version" }, { "wget -O- -i" }, { "" } },
    { {}, { ".ar" }, { "ar --version" }, { "ar -p" }, { "" } },
    { {}, { ".tar" }, { "tar --version" }, { "tar -xOf" }, { "" } },
    { {},
      { ".tar.z", ".tar.gz", ".tgz" },
      { "tar --version" },
      { "tar -xzOf" },
      { "" } },
    { {}, { ".tar.bz2" }, { "tar --version" }, { "tar -xjOf" }, { "" } },
    { {},
      { ".tar.xz" },
      { "tar --version && unxz --version" },
      { "tar --use-compress-program=unxz -xOf" },
      { "" } },
    { {},
      { ".gz", ".z" },
      { "pigz --version", "gzip --version" },
      { "pigz -dc", "gzip -dc" },
      { "pigz >", "gzip >" } },
    { {}, { ".bz2" }, { "bzip2 --version" }, { "bunzip2 -dc" }, { "bzip2 >" } },
    { {}, { ".xz" }, { "xz --version" }, { "unxz -dc" }, { "xz -T0 >" } },
    { {}, { ".7z" }, { "7z" }, { "7z -so e" }, { "7z -si a" } },
    { {}, { ".zip" }, { "zip --version" }, { "unzip -p" }, { "" } },
    { {},
      { ".bam", ".cram" },
      { "samtools --version" },
      { "samtools view -h" },
      { "samtool -Sb - >" } },
  };

  for (const auto& datatype : DATATYPES) {
    bool found_datatype = false;
    for (const auto& prefix : datatype.prefixes) {
      if (starts_with(path, prefix)) {
        found_datatype = true;
        break;
      }
    }
    for (const auto& suffix : datatype.suffixes) {
      if (ends_with(path, suffix)) {
        found_datatype = true;
        break;
      }
    }

    if (found_datatype) {
      bool found_cmd = false;
      int cmd_idx = 0;
      for (const auto& existence_cmd : datatype.cmds_check_existence) {
        bool good = true;
        auto sub_cmds = split(existence_cmd, "&&");
        std::for_each(sub_cmds.begin(), sub_cmds.end(), trim);
        for (const auto& sub_cmd : sub_cmds) {
          auto args = split(sub_cmd, " ");
          std::for_each(args.begin(), args.end(), trim);

          pid_t pid = fork();
          if (pid == 0) {
            int null_fd = orig_open("/dev/null", O_WRONLY, 0);
            dup2(null_fd, STDOUT_FILENO);
            dup2(null_fd, STDERR_FILENO);
            close(null_fd);

            switch (args.size()) {
              case 1:
                execlp(args[0].c_str(), args[0].c_str(), NULL);
              case 2:
                execlp(args[0].c_str(), args[0].c_str(), args[1].c_str(), NULL);
              case 3:
                execlp(args[0].c_str(),
                       args[0].c_str(),
                       args[1].c_str(),
                       args[2].c_str(),
                       NULL);
              case 4:
                execlp(args[0].c_str(),
                       args[0].c_str(),
                       args[1].c_str(),
                       args[2].c_str(),
                       args[3].c_str(),
                       NULL);
              default:
                log_error("Invalid number of arguments supplied to execlp (" +
                          std::to_string(args.size()) + ").");
                std::exit(EXIT_FAILURE);
            }
            log_error("execlp failed.");
            std::exit(EXIT_FAILURE);
          } else {
            int status;
            waitpid(pid, &status, 0);
            if (WIFEXITED(status) && WEXITSTATUS(status) != 0) { // NOLINT
              good = false;
              break;
            }
          }
        }
        if (good) {
          found_cmd = true;
          break;
        }
        cmd_idx++;
      }

      if (found_cmd) {
        std::string cmd;
        if (save) {
          cmd = datatype.save_cmds[cmd_idx];
        } else {
          cmd = datatype.load_cmds[cmd_idx];
        }
        if (cmd.empty()) {
          log_warning("Filetype recognized for '" + path +
                      "', but no tool available to work with it.");
          return "";
        }
        if (cmd.back() == '>') {
          cmd += path;
        } else {
          cmd += " ";
          cmd += path;
        }
        return cmd;
      }
      break;
    }
  }

  return "";
}

/** Open a pipe to save/load the specified file.
 * Not thread safe.
 * @return a file descriptor
 */
static int
data_saveload(const std::string& cmd, const bool save)
{
  static const int READ_END = 0;
  static const int WRITE_END = 1;

  int fd[2];
  check_error(pipe(fd) == -1, "Error opening a pipe.");

  auto args = split(cmd, " ");
  std::for_each(args.begin(), args.end(), trim);

  std::string stdout_to_file;
  decltype(args)::iterator it;
  for (it = args.begin(); it != args.end(); ++it) {
    if (it->front() == '>') {
      stdout_to_file = it->substr(1);
      break;
    }
  }
  if (it != args.end()) {
    args.erase(it);
  }

  pid_t pid = fork();
  check_error(pid == -1, "Error on fork.");

  if (pid == 0) {
    if (save) {
      dup2(fd[READ_END], STDIN_FILENO);
      close(fd[READ_END]);
      close(fd[WRITE_END]);

      if (!stdout_to_file.empty()) {
        int outfd =
          orig_open(stdout_to_file.c_str(),
                    O_WRONLY | O_CREAT,
                    S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
        dup2(outfd, STDOUT_FILENO);
        close(outfd);
      }

      switch (args.size()) {
        case 1:
          execlp(args[0].c_str(), args[0].c_str(), NULL);
        case 2:
          execlp(args[0].c_str(), args[0].c_str(), args[1].c_str(), NULL);
        case 3:
          execlp(args[0].c_str(),
                 args[0].c_str(),
                 args[1].c_str(),
                 args[2].c_str(),
                 NULL);
        case 4:
          execlp(args[0].c_str(),
                 args[0].c_str(),
                 args[1].c_str(),
                 args[2].c_str(),
                 args[3].c_str(),
                 NULL);
        default:
          log_error("Invalid number of arguments supplied to execlp (" +
                    std::to_string(args.size()) + ").");
          std::exit(EXIT_FAILURE);
      }
      log_error("execlp failed.");
      exit(EXIT_FAILURE);
    } else {
      dup2(fd[WRITE_END], STDOUT_FILENO);
      close(fd[READ_END]);
      close(fd[WRITE_END]);

      switch (args.size()) {
        case 1:
          execlp(args[0].c_str(), args[0].c_str(), NULL);
        case 2:
          execlp(args[0].c_str(), args[0].c_str(), args[1].c_str(), NULL);
        case 3:
          execlp(args[0].c_str(),
                 args[0].c_str(),
                 args[1].c_str(),
                 args[2].c_str(),
                 NULL);
        case 4:
          execlp(args[0].c_str(),
                 args[0].c_str(),
                 args[1].c_str(),
                 args[2].c_str(),
                 args[3].c_str(),
                 NULL);
        default:
          log_error("Invalid number of arguments supplied to execlp (" +
                    std::to_string(args.size()) + ").");
          std::exit(EXIT_FAILURE);
      }
      log_error("execlp failed.");
      exit(EXIT_FAILURE);
    }
  } else {
    if (save) {
      close(fd[READ_END]);
      return fd[WRITE_END];
    }
    close(fd[WRITE_END]);
    return fd[READ_END];
  }
  return -1;
}

/** Open a pipe to uncompress the specified file.
 * @return a FILE pointer
 */
static FILE*
fdata_saveload(const std::string& cmd, const bool save)
{
  int fd = data_saveload(cmd, save);
  if (fd == -1) {
    perror(cmd.c_str());
    std::exit(EXIT_FAILURE);
  }
  return fdopen(fd, save ? "w" : "r");
}

extern "C"
{

  int open(const char* __file, int __oflag, ...)
  {
    if ((bool(__oflag | O_RDONLY) || bool(__oflag | O_WRONLY)) &&
        !bool(__oflag | O_RDWR)) {
      bool save = bool(__oflag | O_WRONLY);
      auto cmd = get_saveload_cmd(__file, save);
      if (!cmd.empty()) {
        return data_saveload(cmd, save);
      }
    }
    va_list args;
    va_start(args, __oflag);
    int ret = orig_open(__file, __oflag, va_arg(args, mode_t));
    va_end(args);
    return ret;
  }

  int open64(const char* __file, int __oflag, ...)
  {
    if ((bool(__oflag | O_RDONLY) || bool(__oflag | O_WRONLY)) &&
        !bool(__oflag | O_RDWR)) {
      bool save = bool(__oflag | O_WRONLY);
      auto cmd = get_saveload_cmd(__file, save);
      if (!cmd.empty()) {
        return data_saveload(cmd, save);
      }
    }
    va_list args;
    va_start(args, __oflag);
    int ret = orig_open64(__file, __oflag, va_arg(args, mode_t));
    va_end(args);
    return ret;
  }

  int creat(const char* __file, mode_t __mode)
  {
    auto cmd = get_saveload_cmd(__file, true);
    if (!cmd.empty()) {
      return data_saveload(cmd, true);
    }
    return orig_creat(__file, __mode);
  }

  FILE* fopen(const char* __restrict __filename, const char* __restrict __modes)
  {
    if ((strcmp(__modes, "r") == 0) || (strcmp(__modes, "w") == 0)) {
      bool save = (strcmp(__modes, "w") == 0);
      auto cmd = get_saveload_cmd(__filename, save);
      if (!cmd.empty()) {
        return fdata_saveload(cmd, save);
      }
    }
    return orig_fopen(__filename, __modes);
  }

  FILE* fopen64(const char* __restrict __filename,
                const char* __restrict __modes)
  {
    if ((strcmp(__modes, "r") == 0) || (strcmp(__modes, "w") == 0)) {
      bool save = (strcmp(__modes, "w") == 0);
      auto cmd = get_saveload_cmd(__filename, save);
      if (!cmd.empty()) {
        return fdata_saveload(cmd, save);
      }
    }
    return orig_fopen64(__filename, __modes);
  }

} // extern "C"

} // namespace btllib

#endif
