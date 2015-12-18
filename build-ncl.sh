#!/usr/bin/env bash
 
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

NCL_BUILD_DIR="${ECOEVOLITY_BASE_DIR}/ncl-build"
INSTALL_DIR="${NCL_BUILD_DIR}/installed"

if [ -d "$NCL_BUILD_DIR" ]
then 
    echo "ERROR: build directory '$NCL_BUILD_DIR' already exists."
    echo "To recompile, please remove this directory and re-run this script."
    exit 1
else
    mkdir -p "${NCL_BUILD_DIR}/build"
    mkdir -p "$INSTALL_DIR"
fi

NCL_TAR_BALL="ncl-2.1.18.tar.gz"
NCL_TAR_BALL_PATH="${NCL_BUILD_DIR}/${NCL_TAR_BALL}"

curl -o "$NCL_TAR_BALL_PATH" "http://tcpdiag.dl.sourceforge.net/project/ncl/NCL/ncl-2.1.18/${NCL_TAR_BALL}"
NCL_DIR="${NCL_TAR_BALL_PATH%.tar.gz}"
tar -xzf "$NCL_TAR_BALL_PATH" -C "$NCL_BUILD_DIR"

# number of cpus to use during compile
COMPILETHREADS=2
 
cd "${NCL_BUILD_DIR}/build"
$NCL_DIR/configure --prefix="$INSTALL_DIR"
make clean
make -j $COMPILETHREADS
make install

cd "$ECOEVOLITY_BASE_DIR"

echo
echo
echo NCL headers and binaries are in:
echo "    $INSTALL_DIR"

env_path="${NCL_BUILD_DIR}/ncl-env.sh"
echo "#!/bin/sh" > "$env_path"
echo export PATH="${INSTALL_DIR}/bin:${PATH}" >> "$env_path"
echo export LD_LIBRARY_PATH="${INSTALL_DIR}/lib:${LD_LIBRARY_PATH}" >> "$env_path"
echo export PKG_CONFIG_PATH="${INSTALL_DIR}/lib/pkgconfig:${PKG_CONFIG_PATH}" >> "$env_path"
echo export NCL_PREFIX="${INSTALL_DIR}" >> "$env_path"

rm "$NCL_TAR_BALL_PATH"
rm -r "${NCL_BUILD_DIR}/build"
rm -r "$NCL_DIR"
