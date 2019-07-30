#!/bin/bash

cd "${MESON_SOURCE_ROOT}"

files=$(find include -type f | grep "\(.*\.h$\)\|\(.*\.hpp$\)\|\(.*\.cpp$\)\|\(.*\.cxx$\)")
for file in $files; do
    echo -n "$(realpath $file) "
done