#!/bin/bash

cd "${MESON_SOURCE_ROOT}"
files=$(find include -type f | grep "\(.*\.h$\)\|\(.*\.cpp$\)\|\(.*\.cxx$\)")
echo $files