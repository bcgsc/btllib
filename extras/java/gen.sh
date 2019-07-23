#!/bin/bash

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
cd $SCRIPTPATH/../.. #"${MESON_SOURCE_ROOT}"
ln -sf $PWD/include extras/java/
cd extras/java
swig -java -c++ -Iinclude btl.i
rm include