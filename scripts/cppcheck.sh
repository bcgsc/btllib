#!/bin/bash

files=$(find include -type f | grep "\(.*\.h\)\|\(.*\.cpp\)\|\(.*\.cxx\)")
files+=" "
files+=$(find extras/python -type f | grep "\(.*\.h\)\|\(.*\.cpp\)\|\(.*\.cxx\)")
files+=" "
files+=$(find extras/java -type f | grep "\(.*\.h\)\|\(.*\.cpp\)\|\(.*\.cxx\)")

output=$(cppcheck $files --language=c++ --force)
echo "$output"

errors=$(echo $output | grep -i "error")
if [[ $errors ]]; then
    exit 1
fi
