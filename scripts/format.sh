#!/bin/bash

cd "${MESON_BUILD_ROOT}"

for file in "$@"; do
    diff <(cat $file) <(clang-format -style=mozilla $file)
done