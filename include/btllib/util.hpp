#ifndef BTLLIB_UTIL_HPP
#define BTLLIB_UTIL_HPP

#include <algorithm>
#include <vector>
#include <string>

static inline std::vector<std::string>
split(const std::string& s, const std::string& delim) {
    std::vector<std::string> tokens;
    size_t pos1 = 0, pos2 = 0;
    while ((pos2 = s.find(delim, pos2)) != std::string::npos) {
        tokens.push_back(s.substr(pos1, pos2 - pos1));
        pos2 += delim.size();
        pos1 = pos2;
    }
    tokens.push_back(s.substr(pos1));
    return tokens;
}

static inline void
ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [] (int ch) {
        return !std::isspace(ch);
    }));
}

static inline
void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [] (int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

static inline
void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

static inline bool
starts_with(std::string s, std::string prefix) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    std::transform(prefix.begin(), prefix.end(), prefix.begin(), ::tolower);
    return s.find(prefix) == 0;
};

static inline bool
ends_with(std::string s, std::string suffix) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    std::transform(suffix.begin(), suffix.end(), suffix.begin(), ::tolower);
    auto pos = s.rfind(suffix);
    return (pos != std::string::npos) && (pos == s.size() - suffix.size());
};

#endif