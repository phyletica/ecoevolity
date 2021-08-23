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


Setting up the |phyco| configuration file
=========================================

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

In the ``phycoeval-example-data`` directory, there is a |yaml|_-formatted
config file named ``phycoeval-config.yml``, which contains the following
information::

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
    mutation_parameters:
        mutation_rate:
            value: 1.0
            estimate: false
        freq_1:
            value: 0.5
            estimate: false
    branch_parameters:
        population_size:
            equal_population_sizes: true
            value: 0.0005
            estimate: true
            prior:
                gamma_distribution:
                    shape: 4.0
                    mean: 0.0005
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
    mcmc_settings:
        chain_length: 7500
        sample_frequency: 5

Let's break this down to go over all the settings specified in this config.
First, the ``data`` section.

data
----

::

    data:
        ploidy: 2
        constant_sites_removed: false
        alignment:
            genotypes_are_diploid: true
            markers_are_dominant: false
            population_name_is_prefix: false
            population_name_delimiter: ' '
            path: Cyrtodactylus-tutorial-data.nex

This section tells |phyco| all about our data, including:

#.  the geckos are diploid
    (``ploidy: 2``),
#.  we have not removed constant sites from the alignment,
#.  each cell in the alignment represents the genotype of a diploid individual
    (``genotypes_are_diploid: true``),
#.  our nucleotide data are not dominant (i.e., we can distinguish heterozygous
    genotypes for a site from either of the two homozygous genotypes),
#.  the last part of every sequence label specifies the population/species it
    belongs to (``population_name_is_prefix: false``),
    and
#.  the data are found in a file named ``Cyrtodactylus-tutorial-data.nex`` in
    the current directory

mutation_parameters
-------------------

The ``mutation_parameters`` section specifies that we wish to fix the rate of
mutation to 1::

        mutation_rate:
            value: 1.0
            estimate: false

This means that time will be measured in expected substitutions per site,
and effective population sizes will be scaled by the mutation rate
(:math:`\epopsize\murate`).

We also constrain the equilibrium frequencies of the two possible states for
each character to be equal::

        freq_1:
            value: 0.5
            estimate: false

If you are analyzing nucleotide data, this is almost certainly what you want to
do.
|Phyco| assumes biallelic data, and there are many ways to convert 4-state
nucleotide data into 2-states.
If we don't constrain the frequencies to be equal, our results might vary
depending on how choose to do this conversion.

branch_parameters
-----------------

The ``branch_parameters`` section specifies how we wish to model
the effective population sizes of branches across the tree::

        population_size:
            equal_population_sizes: true
            value: 0.0005
            estimate: true
            prior:
                gamma_distribution:
                    shape: 4.0
                    mean: 0.0005

The ``equal_population_sizes: true`` setting specifies that
all the branches share the same effective population size.
We specify a starting value for this parameter (``value: 0.0005``),
and that we wish to estimate it from the data.
Lastly, we specify a gamma-distributed prior distribution on the population
size with a shape and mean of 4 and 0.0005, respectively.
Because the mutation rate is fixed to 1 (see ``mutation_parameters`` above),
the effective population size is scaled by the mutation rate,
so we are putting a prior on 
:math:`\epopsize\murate`.


tree_model
----------

Next, the ``tree_model`` section specifies how we want |phyco| to model trees::

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

This tells |phyco| to allow "generalized" trees with shared and multifurcating
divergences (i.e., all possible non-reticulating tree models),
and to begin the MCMC chain with the ``comb`` tree (the tree where
all tips diverge from a single internal node).
The only ``tree_prior`` implemented in |phyco| is the
``uniform_root_and_betas`` tree prior,
which is
:ref:`described in more detal here <uniform_root_and_betas>`.

Lastly, this section specified how to model the age of the root node
(``root_height``).
This config specifies that we want to allow the age of the root to vary
(``estimate: true``),
and we are placing an exponentially distributed prior distribution on
it with a mean of 0.01.
Given that the mutation rate is fixed to 1 (see ``mutation_parameters`` section
above), this is in units of expected substitutions per site.

How to calibrate time?
======================

.. collapse:: Scaling to absolute time...
    
    Currently in |phyco| there are only two ways to calibrate time to be in units
    other than substitutions per site (e.g., years or millions of years).
    You can place an informative prior distribution on either the
    root age or the mutation rate (or both).
    If you do this, it is important to remember that the ``root_height``,
    ``mutation_rate``, and ``population_size`` are all interrelated.
    So, changing the setting for one of them, probably requires adjusting
    all three of these settings. Let's walk through an example to
    make this clearer.
    
    Let's say we have prior data about the mutation rate of our Philippine
    bent-toed geckos.
    Based on these prior data, we are very confident the rate of mutation of our
    RADseq loci is 0.001 substitutions per site per million years.
    So, we decide to put an informative prior on the rate of mutation per million
    years::
    
            mutation_rate:
                value: 0.001
                estimate: true
                prior:
                    gamma_distribution:
                        shape: 100.0
                        mean: 0.001
    
    That's great, but now we need to change our settings for the root age and
    population size accordingly.
    Time is now in millions of years, rather than expected substitutions
    per site like it was above when :math:`\murate = 1`.
    Above, when :math:`\murate = 1`, our prior expectation for the age of the root
    was 0.01 substitutions per site.
    Now, we expect the rate of mutation is 0.001 substitutions per million years,
    so our new expectation for the age of the root is
    :math:`0.01/0.001 = 10` million years.
    So we should change our prior on the ``root_height`` to something like::
    
                        root_height:
                            estimate: true
                            prior:
                                exponential_distribution:
                                    mean: 10.0
    
    Also, when :math:`\murate = 1`,
    our prior expectation for the mutation-scaled effective population size
    (:math:`\epopsize\murate`) was 0.0005.
    Now that we expect a rate of mutation of 0.001 substitutions per million
    years, we have to adjust our prior on the population size:
    
    .. math::
    
        \epopsize \times \murate &= 0.0005 \\
        \epopsize \times 0.001   &= 0.0005 \\
        \epopsize &= \frac{0.0005}{0.001} \\
        \epopsize &= 0.5
    
    So, we should also change our prior on ``population_size`` to something like::
    
            population_size:
                equal_population_sizes: true
                value: 0.5
                estimate: true
                prior:
                    gamma_distribution:
                        shape: 4.0
                        mean: 0.5
    
    For more about specifying settings for population size parameter(s),
    :ref:`see here <phycopopsize>`.

Running |phyco|
===============
