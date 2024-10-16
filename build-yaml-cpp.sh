#!/usr/bin/env bash
 
set -e
 
# get location of script
ECOEVOLITY_BASE_DIR=""
this_dir=`dirname "$0"`
if [ "$this_dir" = "." ]
then
    ECOEVOLITY_BASE_DIR="$(pwd)"
else
    cd "$this_dir"
    ECOEVOLITY_BASE_DIR="$(pwd)"
fi

YAML_CPP_BUILD_DIR="${ECOEVOLITY_BASE_DIR}/yaml-cpp-build"
INSTALL_DIR="${YAML_CPP_BUILD_DIR}/installed"

if [ -d "$YAML_CPP_BUILD_DIR" ]
then 
    echo "ERROR: build directory '$YAML_CPP_BUILD_DIR' already exists."
    echo "To recompile, please remove this directory and re-run this script."
    exit 1
else
    mkdir -p "${YAML_CPP_BUILD_DIR}/build"
    mkdir -p "$INSTALL_DIR"
fi

YAML_CPP_DIR="${ECOEVOLITY_BASE_DIR}/src/external/yaml-cpp-master-ce056ac"

# number of cpus to use during compile
COMPILETHREADS=4
 
cd "${YAML_CPP_BUILD_DIR}/build"
cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" "$YAML_CPP_DIR"
make clean
make -j $COMPILETHREADS
make install

cd "$ECOEVOLITY_BASE_DIR"

echo
echo
echo yaml-cpp headers and binaries are in:
echo "    $INSTALL_DIR"

env_path="${YAML_CPP_BUILD_DIR}/yaml-cpp-env.sh"
echo "#!/bin/sh" > "$env_path"
echo export PATH="${INSTALL_DIR}/bin:\${PATH}" >> "$env_path"
echo export LD_LIBRARY_PATH="${INSTALL_DIR}/lib:\${LD_LIBRARY_PATH}" >> "$env_path"
echo export PKG_CONFIG_PATH="${INSTALL_DIR}/lib/pkgconfig:\${PKG_CONFIG_PATH}" >> "$env_path"
echo export YAML_CPP_PREFIX="${INSTALL_DIR}" >> "$env_path"

rm -r "${YAML_CPP_BUILD_DIR}/build"
