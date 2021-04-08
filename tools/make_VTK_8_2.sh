#!/usr/bin/env bash

mkdir build_vtk
cd build_vtk
wget --no-verbose https://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz
tar xf VTK-8.2.0.tar.gz
cd VTK-8.2.0
mkdir build
cd build
cmake .. -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release -DVTK_WRAP_PYTHON=OFF
make -j6
make install
cd ../../../
rm -rf build_vtk
ldconfig -v
