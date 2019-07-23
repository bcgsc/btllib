#!/bin/bash

cd "${MESON_SOURCE_ROOT}"
ln -sf $PWD/include extras/java/
cd extras/java
swig -java -c++ -Iinclude btl.i
rm include