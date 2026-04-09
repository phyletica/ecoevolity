#!/usr/bin/env bash

set -e

# make sure we get back to directory of caller
current_dir="$(pwd)"
function return_on_exit () {
    cd "$current_dir"
}
trap return_on_exit EXIT

usage () {
    echo ""
    echo "usage: $0 [-h|--help] [-d|--test-dir-name <DIR-NAME>] [-e|--exhaustive]"
    echo "  -h|--help        Show help message and exit."
    echo "  -e|--exhaustive  Run FULL test suite. Default: Run short tests."
    echo "  -d|--test-dir-name  Name of build directory for tests. Default: 'test-suite-build'."
    echo ""
}
 
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

test_build_dir="${ECOEVOLITY_BASE_DIR}/test-suite-build"
exhaustive_arg="OFF"

if [ "$(echo "$@" | grep -c "=")" -gt 0 ]
then
    echo "ERROR: Do not use '=' for arguments. For example, use"
    echo "'--test-dir-name testfun' instead of '--test-dir-name=testfun'."
    exit 1
fi

while [ "$1" != "" ]
do
    case $1 in
        -h| --help)
            usage
            exit
            ;;
        -d| --test-dir-name)
            shift
            test_build_dir="$1"
            ;;
        -e| --exhaustive)
            exhaustive_arg="ON"
            ;;
        * )
            extra_args="$extra_args $1"
    esac
    shift
done

args="-DCMAKE_BUILD_TYPE=Debug -DFORCE_BUNDLED_NCL=OFF -DFORCE_BUNDLED_YAML_CPP=OFF -DBUILD_WITH_THREADS=ON -DEXHAUSTIVE_TESTING=${exhaustive_arg}"

if [ -n "$extra_args" ]
then
    args="${args} ${extra_args}"
fi

if [ -d "$test_build_dir" ]
then 
    echo "ERROR: build directory '$test_build_dir' already exists."
    echo "To recompile, please remove this directory and re-run this script."
    exit 1
else
    mkdir -p "${test_build_dir}"
fi

if [ -e "${ECOEVOLITY_BASE_DIR}/dependencies/env-ncl.sh" ]
then
    source "${ECOEVOLITY_BASE_DIR}/dependencies/env-ncl.sh"
fi

if [ -e "${ECOEVOLITY_BASE_DIR}/dependencies/env-yaml-cpp.sh" ]
then
    source "${ECOEVOLITY_BASE_DIR}/dependencies/env-yaml-cpp.sh"
fi

# number of cpus to use during compile
compile_threads=4
 
cd "${test_build_dir}"
cmake  $args ../
make clean
make -j $compile_threads check

cd test
./test_ecoevolity $@
