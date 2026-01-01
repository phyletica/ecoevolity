BootStrap: docker
From: ubuntu:24.04

%post

    export DEBIAN_FRONTEND=noninteractive
    
    echo "Installing required packages..."
    
    # install build dependencies
    apt-get update -y
    apt-get install -y build-essential bash-completion git cmake

%environment
    export PATH=$PATH:/scif/apps/ecoevolity/bin

%appinstall ecoevolity
    git clone https://github.com/phyletica/ecoevolity.git /ecoevolity
    cd /ecoevolity
    git checkout {{ version }}
    ./build.sh --threads --prefix ${SCIF_APPROOT} --build-type relwithdebinfo
