#!/bin/bash

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
cd $SCRIPTPATH/../.. #"${MESON_SOURCE_ROOT}"
ln -sf $PWD/include extras/python/
cd extras/python
swig -python -py3 -c++ -Iinclude btl.i
rm include