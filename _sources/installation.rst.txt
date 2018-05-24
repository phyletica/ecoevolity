.. _installation:

############
Installation
############

..  .. contents::
        :local:
        :depth: 2


.. _prerequisites:

*************
Prerequisites
*************

Compiling |eco|_ requires |cmake|_ and a new-ish C++ compiler (one that
supports the C++11 standard).
If you are on a Linux machine, use your package manager to update/install g++
and cmake.
On a fresh install of Ubuntu 16.04 you can simply::

    $ sudo apt install cmake g++

If you are on a Mac, you can install Xcode command line tools and download and
install |cmake|_ from https://cmake.org/.

We also strongly recommend using |git|_ to acquire the source code.
If you use |git|_, the |eco|_ tools will report the version of
|eco|_ you are using much more precisely, which will make your
work more reproducible.

While not required, the |pyco|_ Python package can be useful for summarizing
the output |eco|_.
This will be used in tutorials.
:ref:`See below <pycoevolity_install>`
for how to install |pyco|_


***********************
Getting the source code
***********************

If you have |git|_ installed, clone the repository::

    $ git clone https://github.com/phyletica/ecoevolity.git

If you prefer not to use |git|_, you can download an archive of the 
`source code here <https://github.com/phyletica/ecoevolity/archive/master.zip>`_.
Once downloaded, you'll need to unzip the archive.


***********************
Basic build and install 
***********************

Next, move into the downloaded directory and run the build script::

    $ cd ecoevolity
    $ ./build.sh

If the build was successful, the |eco|_ executables should now be in the
``./build/bin`` directory, and you should be able to run::

    $ ./build/bin/ecoevolity -h

and see the |eco| help menu, the beginning of which should looks something
like::

    ======================================================================
                                  Ecoevolity
                      Estimating evolutionary coevality
            Version 0.2.0 (dev d4a4d48: 2018-05-23T15:22:40-05:00)
    ======================================================================
    
    Usage: ecoevolity [OPTIONS] YAML-CONFIG-FILE
    
    Ecoevolity: Estimating evolutionary coevality
    
    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit

The executables in ``./build/bin`` are ready to use, but you'll probably want
to put them in your PATH (a list of directories that your shell looks in to
find the commands you type on the command line). You can do this via::

    $ sudo cp ./build/bin/* /usr/local/bin

If this worked, you're good to go; you can try ``ecoeovlity -h`` to be sure.

If it didn't work, you probably don't have admin privileges.
If so, you can create a bin folder in your home folder and put the tools
there::

    $ mkdir -p "${HOME}/bin"
    $ cp ./build/bin/* "${HOME}/bin"

Then, you can add this directory to your PATH (if it's not already there; you
can check with ``echo $PATH``)::

    $ export PATH="${PATH}:${HOME}/bin"

Note, this update to PATH is only for your current terminal window.  If you
want this to be permanent (work for all future terminal windows), add ``export
PATH="${PATH}:${HOME}/bin"`` to your ``.bashrc`` or ``.bash_profile`` file in
your home directory.


********************
Install during build
********************

If you want to build and install in one go, you just need to specify where you
want the installation to go, for example::

    $ sudo ./build.sh --prefix /usr/local


*****************************
Building the threaded version
*****************************

If you want to install a version of |eco|_ that performs the likelihood
calculations across multiple threads, you just need to add the ``--threads``
flag::

    $ ./build.sh --threads

In my opinion, you're usually better off running multiple independent chains
rather than multithreading, but the option is there.


.. _pycoevolity_install:

**********************
Installing pycoevolity
**********************

|Pyco|_ is a Python package for summarizing the output of |eco|_.
It should work with Python 2 or 3, and can be installed via::

    $ pip install git+git://github.com/phyletica/pycoevolity.git
