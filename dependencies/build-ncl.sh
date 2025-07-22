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

ncl_build_dir="${dep_dir}/ncl-build"
install_dir="${ncl_build_dir}/installed"

if [ -d "$ncl_build_dir" ]
then 
    echo "ERROR: build directory '$ncl_build_dir' already exists."
    echo "To recompile, please remove this directory and re-run this script."
    exit 1
else
    mkdir -p "${ncl_build_dir}/build"
    mkdir -p "$install_dir"
fi

ncl_tar_ball="ncl-2.1.18.tar.gz"
ncl_tar_ball_path="${ncl_build_dir}/${ncl_tar_ball}"

curl -L -o "$ncl_tar_ball_path" "https://sourceforge.net/projects/ncl/files/NCL/ncl-2.1.18/${ncl_tar_ball}"
ncl_dir="${ncl_tar_ball_path%.tar.gz}"
tar -xzf "$ncl_tar_ball_path" -C "$ncl_build_dir"

# number of cpus to use during compile
num_threads=4
 
cd "${ncl_build_dir}/build"
$ncl_dir/configure --prefix="$install_dir"
make clean
make -j $num_threads
make install

# separate the static libraries
# mkdir "${install_dir}/lib/static"
# cp "${install_dir}/lib/ncl/libncl.a" "${install_dir}/lib/static"

echo
echo
echo NCL headers and binaries are in:
echo "    $install_dir"
# echo
# echo Static NCL library is located at:
# echo "    ${install_dir}/lib/static/libncl.a"

env_path="${dep_dir}/env-ncl.sh"
echo export PATH="${install_dir}/bin:\${PATH}" > "$env_path"
echo export LD_LIBRARY_PATH="${install_dir}/lib/ncl:\${LD_LIBRARY_PATH}" >> "$env_path"
echo export PKG_CONFIG_PATH="${install_dir}/lib/pkgconfig:\${PKG_CONFIG_PATH}" >> "$env_path"
echo export NCL_PREFIX="${install_dir}" >> "$env_path"

rm "$ncl_tar_ball_path"
rm -r "${ncl_build_dir}/build"
rm -r "$ncl_dir"
