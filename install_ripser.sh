#!/bin/bash

set -e

git clone https://github.com/simonzhang00/ripser-plusplus.git

cd ripser-plusplus/ripserplusplus

mkdir build
cd build

cmake .. \
    -DCMAKE_C_COMPILER="${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-gcc" \
    -DCMAKE_CXX_COMPILER="${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-g++" \
    -DCMAKE_CUDA_HOST_COMPILER="${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-g++"

make -j$(nproc)
