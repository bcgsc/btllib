#!/bin/bash

cd "${MESON_SOURCE_ROOT}"

files=$(find extras/python -type f | grep "\(.*\.h\)\|\(.*\.cpp\)\|\(.*\.cxx\)")
files+=" "
files+=$(find extras/java -type f | grep "\(.*\.h\)\|\(.*\.cpp\)\|\(.*\.cxx\)")
for file in $files; do
    echo -n "$(realpath $file) "
done