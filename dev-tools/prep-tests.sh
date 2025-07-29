#!/usr/bin/env bash

set -e

# make sure we get back to directory of caller
current_dir="$(pwd)"
function return_on_exit () {
    cd "$current_dir"
}
trap return_on_exit EXIT

# get location of script
# NOTE: this does not work with our Slurm submission system
this_dir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo "Loading modules specified in 'modules-to-load.sh'..."
source "${this_dir}/modules-to-load.sh" >/dev/null 2>&1 || echo "  No modules loaded"

base_dir="$(dirname "$this_dir")"

test_build_dir="${base_dir}/test-suite-build"
if [ -n "$1" ]
then
    test_build_dir="$1"
fi

if [ -d "$test_build_dir" ]
then 
    echo "ERROR: build directory '$test_build_dir' already exists."
    echo "Please specify a different directory or remove the existing directory "
    echo "and re-run this script."
    exit 1
else
    if [ ! -d "$(dirname "$test_build_dir")" ]
    then
        echo "ERROR: Cannot make dir for building tests, because the parent directory does not exist: $(dirname "$test_build_dir")"
        exit 1
    fi
    mkdir "${test_build_dir}"
fi

ncl_env_path="${base_dir}/dependencies/env-ncl.sh"
if [ -e "$ncl_env_path" ]
then
    echo "Sourcing '$ncl_env_path'..."
    source "$ncl_env_path"
fi

yaml_cpp_env_path="${base_dir}/dependencies/env-yaml-cpp.sh"
# if [ -e "$yaml_cpp_env_path" ]
# then
#     echo "Sourcing '$yaml_cpp_env_path'..."
#     source "$yaml_cpp_env_path"
# fi

pll_env_path="${base_dir}/dependencies/env-pll.sh"
if [ -e "$pll_env_path" ]
then
    echo "Sourcing '$pll_env_path'..."
    source "$pll_env_path"
fi

# number of cpus to use during compile
nthreads=4
 
echo ""
echo "Building tests..."
echo ""
cd "${test_build_dir}"
cmake -DCMAKE_BUILD_TYPE=Debug -DFORCE_BUNDLED_NCL=OFF -DFORCE_BUNDLED_YAML_CPP=ON -DBUILD_WITH_THREADS=ON "$base_dir"
make clean
make -j $nthreads check

echo ""
echo "----------------------------------------------------------------------"
echo "Test executable is ready to run:"
echo "  ${test_build_dir}/test/test_ecoevolity"
echo "----------------------------------------------------------------------"
echo ""

env_path="${this_dir}/env-test-location.sh"
echo "export ecoevolity_test_exe_path=\"${test_build_dir}/test/test_ecoevolity\"" >> "$env_path"
