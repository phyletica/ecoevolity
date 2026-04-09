#!/usr/bin/env bash
 
set -e

# make sure we get back to directory of caller
current_dir="$(pwd)"
function return_on_exit () {
    cd "$current_dir"
}
trap return_on_exit EXIT

# number of cpus to use during compile
num_threads=4

# get location of script
dep_dir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

base_dir="$(dirname "$dep_dir")"

echo "Loading modules specified in '../modules-to-load.sh'..."
source "${base_dir}/modules-to-load.sh" >/dev/null 2>&1 || echo "  No modules loaded"

ncl_dir="${dep_dir}/ncl-build"
ncl_repo_dir="${ncl_dir}/ncl"
ncl_build_dir="${ncl_dir}/build"
ncl_install_dir="${ncl_dir}/installed"

if [ -d "$ncl_dir" ]
then 
    echo "ERROR: build directory '$ncl_dir' already exists."
    echo "To recompile, please remove this directory and re-run this script."
    exit 1
else
    mkdir "$ncl_dir"
    mkdir "$ncl_build_dir"
    mkdir "$ncl_install_dir"
fi

ncl_commit="64db8d97"

(
    cd "$ncl_dir"
    git clone git@github.com:mtholder/ncl.git "$ncl_repo_dir"
    cd "$ncl_repo_dir"
    git checkout -b project-env "$ncl_commit"

    sh bootstrap.sh
    cd "$ncl_build_dir"
    "${ncl_repo_dir}/configure" --prefix="$ncl_install_dir"
    make clean
    make -j $num_threads
    make install
)

# separate the static libraries
# mkdir "${ncl_install_dir}/lib/static"
# cp "${ncl_install_dir}/lib/ncl/libncl.a" "${ncl_install_dir}/lib/static"

echo
echo
echo NCL headers and binaries are in:
echo "    $ncl_install_dir"
# echo
# echo Static NCL library is located at:
# echo "    ${install_dir}/lib/static/libncl.a"

env_path="${dep_dir}/env-ncl.sh"
echo export PATH="${ncl_install_dir}/bin:\${PATH}" > "$env_path"
echo export LD_LIBRARY_PATH="${ncl_install_dir}/lib/ncl:\${LD_LIBRARY_PATH}" >> "$env_path"
echo export PKG_CONFIG_PATH="${ncl_install_dir}/lib/pkgconfig:\${PKG_CONFIG_PATH}" >> "$env_path"
echo export NCL_PREFIX="${ncl_install_dir}" >> "$env_path"

if [ -n "$ncl_build_dir" ] && [ -d "$ncl_build_dir" ]
then
    echo "Cleaning up by removing build directory '$ncl_build_dir'"
    rm -r "$ncl_build_dir"
fi
if [ -n "$ncl_repo_dir" ] && [ -d "$ncl_repo_dir" ]
then
    echo "Cleaning up by removing ncl repo '$ncl_repo_dir'"
    rm -rf "$ncl_repo_dir"
fi
