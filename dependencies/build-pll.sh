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

pll_stage_dir="${dep_dir}/pll-build"
pll_build_dir="${pll_stage_dir}/build"
install_dir="${pll_stage_dir}/installed"

if [ -d "$pll_stage_dir" ]
then 
    echo "ERROR: build directory '$pll_stage_dir' already exists."
    echo "To recompile, please remove this directory and re-run this script."
    exit 1
else
    mkdir "$pll_stage_dir"
    mkdir "$pll_build_dir"
    mkdir "$install_dir"
fi

pll_name="libpll-2"
pll_dir="${dep_dir}/${pll_name}"
git clone "https://github.com/xflouris/${pll_name}.git"

# number of cpus to use during compile
num_threads=4

cd "$pll_dir"
./autogen.sh
cd "${pll_build_dir}"
$pll_dir/configure --prefix="$install_dir"
make -j $num_threads
make install

echo
echo
echo PLL headers and binaries are in:
echo "    $install_dir"

env_path="${dep_dir}/env-pll.sh"
pkg_config_dir="${install_dir}/pkgconfig"
mkdir "$pkg_config_dir"
pkg_config_path="${pkg_config_dir}/pll.pc"
echo export LD_LIBRARY_PATH="${install_dir}/lib:\${LD_LIBRARY_PATH}" > "$env_path"
echo export PKG_CONFIG_PATH="${pkg_config_dir}:\${PKG_CONFIG_PATH}" >> "$env_path"
echo export PLL_PREFIX="${install_dir}" >> "$env_path"

echo prefix="$install_dir"           > "$pkg_config_path"
echo exec_prefix=\${prefix}         >> "$pkg_config_path"
echo includedir=\${prefix}/include  >> "$pkg_config_path"
echo libdir=\${exec_prefix}/lib     >> "$pkg_config_path"

echo ""                                     >> "$pkg_config_path"
echo Name: PLL                              >> "$pkg_config_path"
echo Description: Phylo likelihood library  >> "$pkg_config_path"
echo Version: 2                             >> "$pkg_config_path"
echo Libs: -L\${libdir} -lpll               >> "$pkg_config_path"
echo Cflags: -I\${includedir}               >> "$pkg_config_path"

rm -r "${pll_build_dir}"
