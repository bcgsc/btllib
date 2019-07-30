#!/bin/bash

cd "${MESON_SOURCE_ROOT}"

files=$(scripts/get_include_files.sh)
files+=" "
files+=$(scripts/get_wrapper_files.sh)
output=$(cppcheck $files --language=c++ --force)
echo "$output"

errors=$(echo $output | grep -i "error")
if [[ $errors ]]; then
    exit 1
fi
