.. _installation:

############
Installation
############


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


************
Installation
************

If you have |git|_ installed, clone the repository::

    $ git clone https://github.com/phyletica/ecoevolity.git

Next, move into the downloaded directory::

    $ cd ecoevolity

To install globally::

    $ sudo ./build.sh --install

To install the threaded version globally::

    $ sudo ./build.sh --install --threads

If you do not have admin privileges, you can install the serial version to your
home directory by::

    $ ./build.sh --install --prefix "$HOME"

or the threaded version to your home directory::

    $ ./build.sh --install --threads --prefix "$HOME"

If the install was successful, and the install directory is in your PATH, you
should be able to view the help menu of |eco|_::

    $ ecoevolity -h

You should see output that looks something like::

    ======================================================================
                                  Ecoevolity
                      Estimating evolutionary coevality
          Version 0.1.0 (master 07eab4d: 2017-10-30T16:01:34-05:00)
    ======================================================================
    
    Usage: ecoevolity [OPTIONS] YAML-CONFIG-FILE
    
    Ecoevolity: Estimating evolutionary coevality
    
    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      --seed=SEED           Seed for random number generator. Default: Set from clock.
      --ignore-data         Ignore data to sample from the prior distribution.
                            Default: Use data to sample from the posterior distribution
      --nthreads=NTHREADS   Number of threads to use for likelihood calculations.
                            Default: 1 (no multithreading). If you are using the
                            '--ignore-data' option, no likelihood calculations
                            will be performed, and so no multithreading is used.
      --prefix=PREFIX       Optional string to prefix all output files.
      --relax-constant-sites
                            By default, if you specify 'constant_sites_removed =
                            true' and constant sites are found, Ecoevolity throws
                            an error. With this option, Ecoevolity will
                            automatically ignore the constant sites and only issue
                            a warning (and correct for constant sites in the
                            likelihood calculation). Please make sure you
                            understand what you are doing when you use this option.
      --relax-missing-sites
                            By default, if a column is found for which there is no
                            data for at least one population, Ecoevolity throws an
                            error. With this option, Ecoevolity will automatically
                            ignore such sites and only issue a warning.
      --relax-triallelic-sites
                            By default, if a DNA site is found for which there is
                            more than two nucleotide states, Ecoevolity throws an
                            error. With this option, Ecoevolity will automatically
                            recode such sites as biallelic and only issue a
                            warning. These sites are recoded by assigning state 0
                            to the first nucleotide found and state 1 to all
                            others. If you do not wish to recode such sites and
                            prefer to ignore them, please remove all sites with
                            more than two nucleotide states from your DNA
                            alignments. NOTE: only alignments of nucleotides are
                            affected by this option, not alignments of standard
                            characters (i.e., 0, 1, 2).
      --dry-run             Do not run analysis; only process and report settings.
