#!/bin/bash

if [ -z "${MESON_SOURCE_ROOT}" ]; then
    echo "[ERROR] This script can only be ran with meson!"
    exit 1
fi

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
