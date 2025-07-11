#!/usr/bin/env bash
 
set -e

# make sure we get back to directory of caller
current_dir="$(pwd)"
function return_on_exit () {
    cd "$current_dir"
}
trap return_on_exit EXIT

# get location of script
dep_dir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
 
yaml_cpp_build_dir="${dep_dir}/yaml-cpp-build"
install_dir="${yaml_cpp_build_dir}/installed"

if [ -d "$yaml_cpp_build_dir" ]
then 
    echo "ERROR: build directory '$yaml_cpp_build_dir' already exists."
    echo "To recompile, please remove this directory and re-run this script."
    exit 1
else
    mkdir -p "${yaml_cpp_build_dir}/build"
    mkdir -p "$install_dir"
fi

ecoevolity_base_dir="$(dirname "$dep_dir")"
yaml_cpp_dir="${ecoevolity_base_dir}/src/external/yaml-cpp-master-ce056ac"

# number of cpus to use during compile
num_threads=4
 
cd "${yaml_cpp_build_dir}/build"
cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="$install_dir" "$yaml_cpp_dir"
make clean
make -j $num_threads
make install

echo
echo
echo yaml-cpp headers and binaries are in:
echo "    $install_dir"

env_path="${dep_dir}/env-yaml-cpp.sh"
echo export LD_LIBRARY_PATH="${install_dir}/lib:\${LD_LIBRARY_PATH}" > "$env_path"
echo export PKG_CONFIG_PATH="${install_dir}/share/pkgconfig:\${PKG_CONFIG_PATH}" >> "$env_path"
echo export YAML_CPP_PREFIX="${install_dir}" >> "$env_path"

rm -r "${yaml_cpp_build_dir}/build"
