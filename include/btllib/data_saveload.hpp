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
#include <string>

#include <dlfcn.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

namespace btllib {

inline FILE*
data_load(const std::string& source);
inline FILE*
data_save(const std::string& sink, bool append);

/** SIGCHLD handler. Reap child processes and report an error if any
 * fail. */
inline void
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

inline bool
data_saveload_init()
{
  struct sigaction action; // NOLINT
  action.sa_handler = sigchld_handler;
  sigemptyset(&action.sa_mask);
  action.sa_flags = SA_RESTART;
  sigaction(SIGCHLD, &action, nullptr);
  return true;
}

static const bool data_saveload_initialized = data_saveload_init();

enum SaveloadOp
{
  READ,
  WRITE,
  APPEND
};

inline std::string
get_saveload_cmd(const std::string& path, const SaveloadOp op)
{
  struct Datatype
  {
    std::vector<std::string> prefixes;
    std::vector<std::string> suffixes;
    std::vector<std::string> cmds_check_existence;
    std::vector<std::string> read_cmds;
    std::vector<std::string> write_cmds;
    std::vector<std::string> append_cmds;
  };

  // clang-format off
  static const Datatype DATATYPES[]{
    { { "http://", "https://", "ftp://" }, {}, { "which wget" }, { "wget -O-" }, { "" }, { "" } },
    { {}, { ".url" }, { "which wget" }, { "wget -O- -i" }, { "" }, { "" } },
    { {}, { ".ar" }, { "which ar" }, { "ar -p" }, { "" }, { "" } },
    { {}, { ".tar" }, { "which tar" }, { "tar -xOf" }, { "" }, { "" } },
    { {}, { ".tar.z", ".tar.gz", ".tgz" }, { "which tar" }, { "tar -xzOf" }, { "" }, { "" } },
    { {}, { ".tar.bz2" }, { "which tar" }, { "tar -xjOf" }, { "" }, { "" } },
    { {}, { ".tar.xz" }, { "which tar && which unxz" }, { "tar --use-compress-program=unxz -xOf" }, { "" }, { "" } },
    { {}, { ".gz", ".z" }, { "which pigz", "which gzip" }, { "pigz -dc", "gzip -dc" }, { "pigz >", "gzip >" }, { "pigz >>", "gzip >>" } },
    { {}, { ".bz2" }, { "which bzip2" }, { "bunzip2 -dc" }, { "bzip2 >" }, { "bzip2 >>" } },
    { {}, { ".xz" }, { "which xz" }, { "unxz -dc" }, { "xz -T0 >" }, { "xz -T0 >>" } },
    { {}, { ".7z" }, { "which 7z" }, { "7z -so e" }, { "7z -si a" }, { "7z -si a" } },
    { {}, { ".zip" }, { "which zip" }, { "unzip -p" }, { "" }, { "" } },
    { {}, { ".bam", ".cram" }, { "which samtools" }, { "samtools view -h" }, { "samtools -Sb - >" }, { "samtools -Sb - >>" } },
  };
  // clang-format on
  std::string default_cmd = "cat ";
  if (op == WRITE) {
    default_cmd += ">" + path;
  } else if (op == APPEND) {
    default_cmd += ">>" + path;
  } else {
    default_cmd += path;
  }

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
            int null_fd = open("/dev/null", O_WRONLY, 0);
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
        switch (op) {
          case READ:
            cmd = datatype.read_cmds[cmd_idx];
            break;
          case WRITE:
            cmd = datatype.write_cmds[cmd_idx];
            break;
          case APPEND:
            cmd = datatype.append_cmds[cmd_idx];
            break;
        }
        if (cmd.empty()) {
          log_warning("Filetype recognized for '" + path +
                      "', but no tool available to work with it.");
          return default_cmd;
        }
        if (cmd.back() == '>') {
          cmd += path;
        } else {
          cmd += " ";
          cmd += path;
        }
        return cmd;
      }
      log_warning("Filetype recognized for '" + path +
                  "', but no tool available to work with it.");
      return default_cmd;
    }
  }

  return default_cmd;
}

inline FILE*
run_saveload_cmd(const std::string& cmd, SaveloadOp op)
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
    if (op == WRITE || op == APPEND) {
      dup2(fd[READ_END], STDIN_FILENO);
      close(fd[READ_END]);
      close(fd[WRITE_END]);

      if (!stdout_to_file.empty()) {
        int outfd =
          open(stdout_to_file.c_str(),
               O_WRONLY | O_CREAT | (op == APPEND ? O_APPEND : 0),
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
    if (op == WRITE || op == APPEND) {
      close(fd[READ_END]);
      return fdopen(fd[WRITE_END], "w");
    }
    close(fd[WRITE_END]);
    return fdopen(fd[READ_END], "r");
  }
  return nullptr;
}

inline FILE*
data_load(const std::string& source)
{
  if (source == "-") {
    return stdin;
  }
  const auto cmd = get_saveload_cmd(source, READ);
  check_error(cmd.empty(), "Error loading from " + source);
  return run_saveload_cmd(cmd, READ);
}

inline FILE*
data_save(const std::string& sink, bool append)
{
  if (sink == "-") {
    return stdout;
  }
  const auto cmd = get_saveload_cmd(sink, append ? APPEND : WRITE);
  check_error(cmd.empty(), "Error saving to " + sink);
  return run_saveload_cmd(cmd, append ? APPEND : WRITE);
}

} // namespace btllib

#endif
