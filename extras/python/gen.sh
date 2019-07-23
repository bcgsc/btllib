#!/bin/bash

cd "${MESON_SOURCE_ROOT}"
ln -sf $PWD/include extras/python/
cd extras/python
swig -python -py3 -c++ -Iinclude btl.i
rm include