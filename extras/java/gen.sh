#!/bin/bash

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
cd $SCRIPTPATH/../.. #"${MESON_SOURCE_ROOT}"
ln -sf $PWD/include extras/java/
cd extras/java
swig -java -c++ -Iinclude btl.i
rm include

java_files=$(find . -maxdepth 1 -iname "*.java" | sed "s~\(.*\)~'\1'~" | tr '\n' , | sed 's~.$~~')
meson_build=$(cat meson.build)
echo "$meson_build" | sed "s~java_files.*=.*\[.*\]~java_files = [$java_files]~" >meson.build