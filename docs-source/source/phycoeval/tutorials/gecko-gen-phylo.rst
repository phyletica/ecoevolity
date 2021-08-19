|phyco_logo_long|

.. _gecko_phy_tutorial:

********************************
Gecko generalized phylogenetics 
********************************

Objective
=========

Our goal is to learn how to set up, run, and interpret an |phyco| analysis in
order to infer phylogenies with shared and multifurcating divergences.

Background Information
======================

Biogeographers are often interested in understanding how environmental changes
affect diversification.
When such processes simultaneously affected multiple ancestral species,
this can produce patterns of shared divergences.
|Phyco| provides a tool for testing such predictions using a fully phylogenetic
Bayesian model choice approach
:cite:`Oaks2021phycoeval`.

For more detailed background information about the method implemented in |phyco|,
please :ref:`see here <phycobackground>`.

Software Used in this Activity
==============================

For this activity, we will be using |phyco| (part of the |eco| package) and
|pyco|.
If you haven't already done so, please follow the
:ref:`instructions for installing ecoevolity <installation>`.
Either installing |eco| directly or using the Docker image will work
for this tutorial.

While not required, you may want to install the programs
|Tracer|_ (https://github.com/beast-dev/tracer/releases)
and
|Figtree|_ (https://github.com/rambaut/figtree/releases).
|Tracer|_
is a nice tool for visualizing the mixing and convergence behavior of Markov
chain Monte Carlo (MCMC) analyses,
and
|Figtree|_ is a nice tool for visualizing phylogenetic trees.

The Data
========

We will be analyzing 50 RADseq loci from 7 insular populations of bent-toed
geckos (*Cyrtodactylus*) from the Philippines.
This is a small subset of the data analyzed in :cite:`Oaks2021phycoeval`.

We will use ``git`` to download these data, along with a |phyco|
configuration file.
To do this open your preferred shell (command-line) environment (e.g., Bash),
and type the following two commands to download a directory ("folder")
with the data
and move into it::

    git clone https://github.com/joaks1/phycoeval-example-data.git
    cd phycoeval-example-data

.. admonition:: Extra steps for Docker users

    If you are using the Docker container to run |phyco|, follow these two extra
    steps to be able to analyze the example files you just downloaded.
    First, make sure you are **NOT** inside the Docker container before
    proceeding; you want to be at your computer's command-line console.
    After you've used ``cd`` to navigate into the ``phycoeval-example-data``
    directory we created above, execute the following command to enter the
    Docker container::
    
        docker run -v "$(pwd)":/portal -it phyletica/ecoevolity-docker bash
    
    .. note::
    
        Depending on your system and how Docker is configured, you may need to use
        ``sudo`` to run Docker commands. If you received a "permission denied"
        message when you ran the command above, try::
        
            sudo docker run -v "$(pwd)":/portal -it phyletica/ecoevolity-docker bash
    
    Then, once inside the Docker container, type::
    
        cd portal

    For the remainder of the tutorial, enter (or copy and paste) any
    command-line instructions into the command line of the Docker container.
    However, you can view and edit the input and output files on your computer
    (outside the container) with your plain text editor of choice.

When you list the contents of the ``phycoeval-example-data`` using ``ls``,
you should see the following files:

#.  ``Cyrtodactylus-tutorial-data.nex``

    This file contains nexus-formatted RADseq data from 7 *Cyrtodactylus*
    lizards from 7 populations.
    For detailed information about how to format nexus files for |phyco|,
    :ref:`click here <phycodata>`.

#.  ``phycoeval-config.yml``

    This is a |yaml|_-formatted configuration file that will tell
    |phyco| how to analyze the data.

#.  ``README.md``

    You can ignore this file.


Biallelic characters
--------------------

|Phyco| assumes all of your characters have at most two states (biallelic).
If we provide |phyco| with nucleotide data, it will automatically recode the data
as biallelic by considering the first character it finds in each column as
state ``0``, and if it finds a second state it considers it state ``1``.
If it finds a third state in a column, it will report an error message and
quit.
For characters (columns) with more than two states, you have two options:

#.  remove these columns from your matrix, or
#.  recode them as biallelic.

|Phyco| will do the latter for you (more on this in a bit).

Linked characters
-----------------

|Phyco| also assumes all your characters are unlinked.
Based on analyses of data simulated with 
linked characters :cite:`Oaks2021phycoeval`, we recommend that you
analyze all of your characters (including the constant ones) and violate the
assumption of unlinked characters.
In short, |phyco| performs better when you use all of the sites (including the
constant ones) compared to reducing the data to only one variable character per
locus.
The example data sets we'll be analyzing consist of 50 loci each comprising
about 90 linked sites.


Running |phyco|
===============

Setting up the configuration file
---------------------------------

Once we have our nexus-formatted data ready, the next step in an |phyco|
analysis is setting up the configuration file.
|Phyco| requires a |yaml|_-formatted configuration file.
|yaml|_ is a human-friendly data standard that allows you to provide |phyco|
the information it needs in a format that is easy for you to read and edit.

.. note::

    The website |yamllint| is a nice tool for debugging |yaml|_ syntax.
    You can copy and paste your config file there to check if you're using
    valid |yaml| syntax.


For detailed information about all the settings that can be included in an
|phyco| config file,
:ref:`click here <phycoconfigfile>`.

In our directory, there is a |yaml|_-formatted config file named
``phycoeval-config.yml``, which contains the following information::

    ---
    data:
        ploidy: 2
        constant_sites_removed: false
        alignment:
            genotypes_are_diploid: true
            markers_are_dominant: false
            population_name_is_prefix: false
            population_name_delimiter: ' '
            path: Cyrtodactylus-tutorial-data.nex
    tree_model:
        tree_space: generalized
        starting_tree: comb
        tree_prior:
            uniform_root_and_betas:
                parameters:
                    root_height:
                        estimate: true
                        prior:
                            exponential_distribution:
                                mean: 0.01
    branch_parameters:
        population_size:
            equal_population_sizes: true
            value: 0.0005
            estimate: true
            prior:
                gamma_distribution:
                    shape: 4.0
                    mean: 0.0005
    mutation_parameters:
        freq_1:
            value: 0.5
            estimate: false
        mutation_rate:
            value: 1.0
            estimate: false
    mcmc_settings:
        chain_length: 7500
        sample_frequency: 5


Changing the prior on the root age
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
