#!/bin/sh

usage () {
    echo ""
    echo "usage: $0 [-h|--help] [-p|--prefix <INSTALL-PREFIX>] [-b|--build-type <BUILD-TYPE>] [-s|--static] [-t|--threads]"
    echo "  -h|--help        Show help message and exit."
    echo "  -s|--static      Build statically linked binaries. Default: dynamic."
    echo "  -t|--threads     Compile with multithreading. Default: No multithreading."
    echo "  -p|--prefix      Install path. Default: './build'."
    echo "  -b|--build-type  Build type. Options: debug, release, relwithdebinfo."
    echo "                   Default: relwithdebinfo."
    echo ""
}

build_ecoevolity () {
    cmake_args=$@
    if [ -n "$ECOEVOLITY_BASE_DIR" ]
    then
        base_dir="$ECOEVOLITY_BASE_DIR"
    else
        echo "ERROR: ECOEVOLITY_BASE_DIR was not defined prior to calling "
        echo "build_ecoevolity."
        exit 1
    fi
    if [ -n "$ECOEVOLITY_BUILD_DIR" ]
    then
        build_dir="$ECOEVOLITY_BUILD_DIR"
    else
        build_dir="${base_dir}/build"
    fi
    if [ -d "$build_dir" ]
    then
        echo "ERROR: build directory '$build_dir' already exists."
        echo "To reconfigure, please remove this directory and re-run this script."
        exit 1
    else
        mkdir -p "$build_dir"
    fi
    
    # configure make files and build
    cd "$build_dir"
    echo "Configuring make files..."
    echo "${base_dir}/ $cmake_args" | xargs cmake || exit 1
    echo "Building..."
    make clean || exit 1
    make -j $COMPILETHREADS || exit 1
    echo "Installing..."
    make install || exit 1
    cd "$base_dir"
    echo "Done!"
}


# get location of script
ECOEVOLITY_BASE_DIR=""
this_dir=`dirname "$0"`
if [ "$this_dir" = "." ]
then
    export ECOEVOLITY_BASE_DIR="$(pwd)"
else
    cd "$this_dir"
    export ECOEVOLITY_BASE_DIR="$(pwd)"
fi
export ECOEVOLITY_BUILD_DIR="${ECOEVOLITY_BASE_DIR}/build"

# make sure submodules are here and up to date
# git submodule init
# git submodule update

# process args
extra_args=""
static=""
threads=""
absolute_root_size=""
install_prefix="$ECOEVOLITY_BUILD_DIR"
universal_mac_build=""
linux_dist_build=""
mbit=64
build_type="RELWITHDEBINFO"
COMPILETHREADS="4"

if [ "$(echo "$@" | grep -c "=")" -gt 0 ]
then
    echo "ERROR: Do not use '=' for arguments. For example, use"
    echo "'--prefix /usr/local' instead of '--prefix=/usr/local'."
    exit 1
fi

while [ "$1" != "" ]
do
    case $1 in
        -h| --help)
            usage
            exit
            ;;
        -p| --prefix)
            shift
            install_prefix="$1"
            ;;
        -b| --build-type)
            shift
            build_type="$1"
            ;;
        -s| --static)
            static=1
            ;;
        -t| --threads)
            threads=1
            ;;
        -j)
            shift
            COMPILETHREADS="$1"
            ;;
        -m)
            shift
            mbit="$1"
            ;;
        --universal-mac-build)
            universal_mac_build=1
            static=1
            ;;
        --linux-dist-build)
            linux_dist_build=1
            static=1
            ;;
        --absolute-root-size)
            absolute_root_size=1
            ;;
        * )
            extra_args="$extra_args $1"
    esac
    shift
done
export COMPILETHREADS
args="-DCMAKE_BUILD_TYPE=${build_type} -DFORCE_BUNDLED_NCL=yes -DFORCE_BUNDLED_YAML_CPP=yes -DCMAKE_INSTALL_PREFIX=${install_prefix}"
if [ -n "$static" ]
then
    args="${args} -DSTATIC_LINKING=yes"
fi
if [ -n "$threads" ]
then
    args="${args} -DBUILD_WITH_THREADS=yes"
fi
if [ -n "$absolute_root_size" ]
then
    args="${args} -DBUILD_WITH_ABSOLUTE_ROOT_SIZE=yes"
fi
if [ -n "$extra_args" ]
then
    args="${args} ${extra_args}"
fi
if [ -n "$linux_dist_build" ]
then
    args="${args} -DCMAKE_CXX_FLAGS=-m${mbit} -DCMAKE_LD_FLAGS=-m${mbit}"
fi
if [ -n "$universal_mac_build" ]
then
    args="${args} -DCMAKE_CXX_FLAGS=\"-arch ppc -arch i386 -arch x86_64\""
fi

# check for build directory
build_ecoevolity $args
