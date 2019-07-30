#!/bin/bash

cd "${MESON_SOURCE_ROOT}"

files=$(scripts/get_include_files.sh)
for file in $files; do
    diff <(cat $file) <(clang-format -style=mozilla $file)
done