#!/bin/bash

set -e

# Remove old wrapper files
rm -f ${MESON_SOURCE_ROOT}/extras/java/*.java
rm -f ${MESON_SOURCE_ROOT}/extras/python/*.py

# Generate java swig files
cd ${MESON_SOURCE_ROOT}
ln -sf $PWD/include extras/java/
cd extras/java
swig -java -c++ -Iinclude btl.i
rm include

# Add java files to meson build
java_files=$(find . -maxdepth 1 -iname "*.java" | sed "s~\(.*\)~'\1'~" | tr '\n' , | sed 's~.$~~')
meson_build=$(cat meson.build)
echo "$meson_build" | sed "s~java_files.*=.*\[.*\]~java_files = [$java_files]~" >meson.build

# Generate python swig files
cd ${MESON_SOURCE_ROOT}
ln -sf $PWD/include extras/python/
cd extras/python
swig -python -py3 -c++ -Iinclude btl.i
rm include