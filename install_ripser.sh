#!/bin/bash
# Exit immediately if a command exits with a non-zero status
set -e

echo "=================================================="
echo "   Starting Automated Compilation for Ripser++"
echo "   (CUDA 10.1 + GCC 8 Compatibility Mode)"
echo "=================================================="

# Define installation paths
INSTALL_DIR="./"
REPO_DIR="${INSTALL_DIR}/ripser-plusplus"

# 1. Check and create the dedicated Conda build environment
if ! conda env list | grep -q "build_env"; then
    echo "Creating isolated build environment (build_env)..."
    conda create -n build_env -c conda-forge gcc_linux-64=8 gxx_linux-64=8 cmake make git -y
else
    echo "Detected existing build_env environment."
fi

# 2. Clone the repository
echo "Preparing source code directories..."
mkdir -p "${INSTALL_DIR}"
if [ ! -d "${REPO_DIR}" ]; then
    git clone --recursive https://github.com/simonzhang00/ripser-plusplus.git "${REPO_DIR}"
fi

cd "${REPO_DIR}/ripserplusplus"
mkdir -p build
cd build
rm -rf *

# 3. Automatically detect system CUDA runtime library path
echo "Searching for system libcudart.so..."
CUDA_LIB=$(find /usr /usr/local -name libcudart.so 2>/dev/null | head -n 1)

if [ -z "$CUDA_LIB" ]; then
    # Default back to your known system path if search fails
    CUDA_LIB="/usr/lib/x86_64-linux-gnu/libcudart.so"
fi
echo "Using CUDA Library Path: $CUDA_LIB"

# 4. Compile safely in the background using conda run to bypass C++17 flags
echo "Compiling project (Clearing Conda environment conflicts)..."
conda run --no-capture-output -n build_env bash -c "
    unset CXXFLAGS
    unset CFLAGS
    cmake -DCMAKE_C_COMPILER=\$CC \
          -DCMAKE_CXX_COMPILER=\$CXX \
          -DCMAKE_CUDA_HOST_COMPILER=\$CXX \
          -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
          -DCUDA_CUDART_LIBRARY=$CUDA_LIB \
          -DCMAKE_CXX_STANDARD=14 \
          -DCMAKE_CUDA_STANDARD=14 \
          ..
    make -j 4
"

echo "=================================================="
echo "Ripser++ automated compilation successful!"
echo "Executable: ${REPO_DIR}/ripserplusplus/build/ripser++"
echo "=================================================="

# 5. Verify the compilation
${REPO_DIR}/ripserplusplus/build/ripser++ --help
