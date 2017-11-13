.. _installation:

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
            Version 0.1.0 (dev a5f964c: 2017-11-09T14:15:52-06:00)
    ======================================================================
    
    Usage: ecoevolity [OPTIONS] YAML-CONFIG-FILE
    
    Ecoevolity: Estimating evolutionary coevality
    
    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      --seed=SEED           Seed for random number generator. Default: Set from clock.
      --ignore-data         Ignore data to sample from the prior distribution.
                            Default: Use data to sample from the posterior distribution
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
    
