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


Getting the example files for this tutorial
===========================================

We will use ``git`` to download the data and |phyco| configuration file we need
for this tutorial.
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

    This file contains nexus-formatted RADseq data comprising 50 loci from 7
    bent-toed geckos (*Cyrtodactylus*) from insular populations in the
    Philippines.
    For detailed information about how to format nexus files for |phyco|,
    :ref:`click here <phycodata>`.

#.  ``phycoeval-config.yml``

    This is a |yaml|_-formatted configuration file that will tell
    |phyco| how to analyze the data.

#.  ``README.md``

    You can ignore this file.


The Data
========

We will be analyzing 50 RADseq loci from 7 insular populations of bent-toed
geckos (*Cyrtodactylus*) from the Philippines.
This is a small subset of the data analyzed in :cite:`Oaks2021phycoeval`.
|Phyco|_ assumes our characters have two-states (bialleleic) and are unlinked.

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
----------------------

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

Now that we have our nexus-formatted data file and
|yaml|_-formatted config file ready,
now we will use |phyco| to infer the phylogeny of our 7 geckos, including
potentially shared or multifurcating divergences.
First let's take a look at the help menu of |phyco|::

    phycoeval -h

If |phyco| is correctly installed, you should get output that looks something
like::

    ======================================================================
                                  Phycoeval
                      Estimating phylogenetic coevality
    
                                   Part of:
                                  Ecoevolity
           Version 0.3.2 (docs 7be2cdd: 2021-08-17T15:46:18-05:00)
    ======================================================================
    
    Usage: phycoeval [OPTIONS] YAML-CONFIG-FILE
    
    Phycoeval: Estimating phylogenetic coevality
    
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
                            true' and constant sites are found, phycoeval throws
                            an error. With this option, phycoeval will
                            automatically ignore the constant sites and only issue
                            a warning (and correct for constant sites in the
                            likelihood calculation). Please make sure you
                            understand what you are doing when you use this option.
      --relax-missing-sites
                            By default, if a column is found for which there is no
                            data for at least one population, phycoeval throws an
                            error. With this option, phycoeval will automatically
                            ignore such sites and only issue a warning.
      --relax-triallelic-sites
                            By default, if a DNA site is found for which there is
                            more than two nucleotide states, phycoeval throws an
                            error. With this option, phycoeval will automatically
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

In addition to the settings in our YAML config file, this help menu shows us
what options for |phyco| we can specify on the command line.
**NOTE**: if you did not compile |eco| to allow multi-threading, you will not
see the ``--nthreads`` example.


let's try running an analysis with |phyco|::

    phycoeval phycoeval-config.yml

Oops, we got an error::

    #######################################################################
    ###############################  ERROR  ###############################
    37 sites from the alignment in:
        'Cyrtodactylus-tutorial-data.nex'
    have no data for at least one population.
    #######################################################################

This error message is telling us that we have 37 sites (columns) in our
character matrix for which we have no data for at least one population.
Such sites can be common in RADseq loci, because most assemblers enforce
thresholds on missing data at the locus level (not for each site).
Rather than removing these sites ourselves, we can tell |phyco| to ignore
them.

Before we do that, let's talk about another common error you might see with
your own data that will look something like (you won't get this error for the
example data)::

    #######################################################################
    ###############################  ERROR  ###############################
    7 sites from the alignment in:
        'Cyrtodactylus-tutorial-data.nex'
    have more than two character states.
    #######################################################################

This error message is telling us that we have 7 sites (columns) in our
character matrix where there are more than two states represented across our
samples.
As discussed above, at this point, we have 2 options:

#.  Remove any characters with more than two states.
#.  Recode these sites as biallelic.

The latter |phyco| will do for us if we specify the
``--relax-triallelic-sites`` option in our command.
When we use this option, |phyco| will consider the first state in a column as
``0``, and any other state found in the column as ``1``.

OK, let's try running |phyco| again, but tell it to ignore the
37 sites for which at least one population has no data::

    phycoeval --relax-missing-sites phycoeval-config.yml

Now, you should be running!
How long the analysis takes to run will depend on your computer, but
it took about 2 minutes on my laptop.

Checking the pre-MCMC screen output
-----------------------------------

While the analysis is running, let's scroll up in our console and look at the
information |phyco| reported before the MCMC chain started sampling.
This output goes by very quickly, but it's **very** important to read it
over to make sure |phyco| is doing what you think it's doing.
At the very top, you will see something like::

    ======================================================================
                                  Phycoeval
                      Estimating phylogenetic coevality
    
                                   Part of:
                                  Ecoevolity
           Version 0.3.2 (docs 7be2cdd: 2021-08-17T15:46:18-05:00)
    ======================================================================
    
    Seed: 1384420509
    Using data in order to sample from the posterior distribution...
    Config path: phycoeval-config.yml

This gives us detailed information about the version of |eco| we are using
(right down to the specific commit SHA).
It also tells us what number was used to seed the random number generator; we
would need to provide |phyco| this seed number using the ``--seed``
command-line option in order to reproduce our results.
The output is also telling us that |phyco| is using the data to sample from the
posterior, as opposed to ignoring the data to sample from the prior.
We can do the latter by specifying the ``--ignore-data`` command-line option.

Next, |phyco| reports a fully specified YAML-formatted configuration file to
the screen.
It's always good to look this over to make sure |phyco| is configured as you
expect.

Next, you will see a summary of the data::

    ----------------------------------------------------------------------
    Summary of data from 'Cyrtodactylus-tutorial-data.nex':
        Genotypes: diploid
        Markers are dominant? false
        Number of populations: 7
        Number of sites: 4538
        Number of variable sites: 102
        Number of patterns: 29
        Patterns folded? true
        Population label (max # of alleles sampled):
            philippinicusSibuyan (2)
            philippinicusPanay (2)
            philippinicusNegros (2)
            philippinicusMindoroELR (2)
            philippinicusMindoroRMB (2)
            philippinicusTablas (2)
            philippinicusLuzonCamarinesNorte (2)
    ----------------------------------------------------------------------


It is **very important** to look over this summary to make sure |phyco|
"sees" your data the way you expect it should.
Make sure the number of sites, the number of populations, and the number of
alleles per population match your expectations!
It is very easy to have typos in your population labels that result in extra
tip populations/species you didn't intend.

Next, |phyco| reports how many threads it is using and where it is
logging information about sampled trees, parameter values, and MCMC
operator statistics::

    Number of threads: 1
    Tree log path: phycoeval-config-trees-run-1.nex
    State log path: phycoeval-config-state-run-1.log
    Operator log path: phycoeval-config-operator-run-1.log


Next, |phyco| reports a summary of the state of the MCMC chain to the screen for
every 10th sample it's logging to the output files.
After the chain finishes, |phyco| also reports statistics associated with all
of the MCMC operators; this information is also logged to the file specified
above as the ``Operator log path``.

The Output Files
----------------

After the analysis finishes, type ``ls`` at the command line::

    ls

You should see three log files that were created by |phyco|:

#.  ``phycoeval-config-trees-run-1.nex``
#.  ``phycoeval-config-state-run-1.log``
#.  ``phycoeval-config-operator-run-1.log``

If you open the ``operator`` log with a plain text editor, you'll see this file
contains information about the MCMC operators.
If you have convergence or mixing issues, this information can be useful for
troubleshooting.

If you open the ``state`` file, you'll see it contains the log likelihood,
prior, and parameter values associated with each MCMC sample.


Running additional chains
=========================

Before we start summarizing the results, let's run a second chain so we can
assess convergence and boost our posterior sample size.
If you have multiple cores on your computer and you've compiled |eco|_ to
allow multi-threading, let's tell |phyco| to use two threads for the
second MCMC analysis::

    phycoeval --nthreads 2 --relax-missing-sites phycoeval-config.yml

Now's a good time to go grab a coffee, while you wait for the second chain to
finish running.
With 2 threads, it takes about 70 seconds on my laptop.

After the chain finishes, if you use ``ls`` you should see the
following ``run-2`` log files:

#.  ``phycoeval-config-trees-run-2.nex``
#.  ``phycoeval-config-state-run-2.log``
#.  ``phycoeval-config-operator-run-2.log``

Two MCMC chains is enough to proceed with this tutorial, but for
"real" analyses, I recommend doing more.
It increases your sample size from the posterior distribution and allows you to
better assess convergence.

I ran 2 additional chains (4 total) on my computer, however, you can proceed
with the steps below with only 2.


Summarizing the Results
=======================

Assessing convergence and mixing
--------------------------------

Before summarizing our posterior sample of trees and other parameter values, we
will use the |sumphyco| command-line tool to assess whether our MCMC chains
converged and mixed well, and if so, how many samples should we ignore from the
beginning of each chain before they converged (i.e., "burn-in" samples).
If you enter the following command::

    sumphycoeval -h

you should see the help menu of |sumphyco| printed to your terminal screen,
showing all the command-line options available for this tool.
The first option we will use to help assess convergence and mixing is the
``-c`` option.
This tells |sumphyco| that we want it to report statistics for assessing
convergence and mixing for different levels of burn-in (i.e., ignoring
different numbers of samples from the beginning of each MCMC chain).
The number we provide after ``-c`` tells |sumphyco| the interval (in
numbers of samples) we want between the different burn-in values.
Enter the following command to have |sumphyco| create a tab-delimited table
of convergence statistics where each row is ignoring a larger number of samples
from the beginning of each chain, incrementing by 100 samples::

    sumphycoeval -c 100 phycoeval-config-trees-run-?.nex > convergence-stats.tsv

After this command finishes, you can open the ``convergence-stats.tsv`` with
a text editor or any spreadsheet program (e.g., Excel) to view the table of
convergence statistics.
It should look something like:

.. csv-table::
    :file: ../../_static/convergence-stats.tsv
    :widths: auto
    :width: 99%
    :header-rows: 1
    :delim: tab
    :align: left
    :class: scrollwide

Your numbers will be different, because MCMC is stochastic and you may have run
more or fewer chains (I ran 4).

|Sumphyco| reports the effective
sample size (ESS) and potential scale reduction factor (PSRF)
for the tree length, root height (age), and the effective
popualtion size of the root branch.
The ESS estimates how many effectively independent samples the
MCMC chains collected for a parameter; the larger the number
the better.
The PSRF measures whether independent MCMC chains are sampling overlapping
values for a parameter (i.e., are they sampling from the same distribution);
values close to one (e.g., < 1.2) indicate the chains converged and are
sampling from the same distribution.
|Sumphyco| also reports the average standard deviation of split frequencies
(ASDSF), which should be close to zero (e.g., < 0.01) if the chains converged
and are sampling the same region of tree space.

Based on the table created by |sumphyco| above, it looks like my MCMC chains
converged quickly and mixed reasonably well.
Ignoring the first 201 samples seems like a good choice as it corresponds with
high ESS values, PSRF values close to 1, and PSRF < 0.01.
If you have |Tracer|_ installed, you can use it to open the
``phycoeval-config-state-run-?.log``
files and look at the convergence and mixing behavior of the chains.

Summarizing the posterior samples of trees and other parameters
---------------------------------------------------------------

We have assessed the convergence and mixing of our
MCMC chains and selected the number of samples we want to ignore
from the beginning of each chain as "burn-in";
I selected 201 samples, but you might choose a different number based on the
metrics from your chains.
Now, we are ready to summarize our post burn-in MCMC samples.
We will use |sumphyco| again to do this.
Enter the following command, but feel free to adjust your burn-in ``-b`` value::

    sumphycoeval -b 201 --map-tree-out cyrt-map-tree.nex phycoeval-config-trees-run-?.nex > posterior-summary.yml

.. note::

    The units for the burn-in option (``-b``/``--burnin``) is the number of
    samples collected by each MCMC chain, **NOT** the number of MCMC
    generations.

This tells |sumphyco| to summarize the sampled trees (ignoring the first 201
from each chain), write the maximum *a posteriori* (MAP) tree to a file named
``cyrt-map-tree.nex``, and write a |yaml|_-formatted summary of the posterior
to a file named ``posterior-summary.yml``.

YAML posterior summary
^^^^^^^^^^^^^^^^^^^^^^

By default, |sumphyco| writes a |yaml|_-formatted summary of the posterior
sample of trees to standard output (i.e., to the terminal's screen).
In the command above, we used ``> posterior-summary.yml`` to redirect
this summary output to a file name ``posterior-summary.yml``.
I chose to use the |yaml|_ format for the posterior summary, because
it is human-friendly to read and is a standard format, so most
modern programming languages have a package for parsing it.

The |yaml|_-formatted posterior summary is very rich with information.  We will
walk through some of the ``posterior-summary.yml`` file to see what kinds of
information it contains.

**leaf_label_map**
""""""""""""""""""

It begins with a ``leaf_label_map`` that tells you the numbers that correspond
to the leaf (tip) labels in the trees::

    leaf_label_map:
        0: philippinicusLuzonCamarinesNorte
        1: philippinicusMindoroELR
        2: philippinicusMindoroRMB
        3: philippinicusNegros
        4: philippinicusPanay
        5: philippinicusSibuyan
        6: philippinicusTablas

**summary_of_tree_sources**
"""""""""""""""""""""""""""

These numbers will be used to refer to the leaf labels throughout the rest of
the file.
Next, it provides a summary about where all the tree samples came from::

    summary_of_tree_sources:
        total_number_of_trees_sampled: 5200
        sources:
            -
              path: phycoeval-config-trees-run-1.nex
              number_of_trees_skipped: 201
              number_of_trees_sampled: 1300
            -
              path: phycoeval-config-trees-run-2.nex
              number_of_trees_skipped: 201
              number_of_trees_sampled: 1300
            -
              path: phycoeval-config-trees-run-3.nex
              number_of_trees_skipped: 201
              number_of_trees_sampled: 1300
            -
              path: phycoeval-config-trees-run-4.nex
              number_of_trees_skipped: 201
              number_of_trees_sampled: 1300

**summary_of_split_freq_std_deviations**
""""""""""""""""""""""""""""""""""""""""

Next is a summary of the standard deviations of the frequencies of all tree
splits (clades) with a frequency greater than 10% (``min_frequency: 0.1``; this
can be changed with the ``--min-split-freq`` option in |sumphyco|)::

    summary_of_split_freq_std_deviations:
        min_frequency: 0.100000000000000006
        n: 8
        mean: 0.00934072992015095964
        max: 0.0230469891327866362

More specifically, for each split (clade) that is present in more
than 10% of the posterior samples, the standard deviation of the
frequencies across all the MCMC chain is calculated.
Above, we are shown the mean and max of all these standard deviations of split
frequencies (SDSF).
If all the chains are sampling the same distribution of trees, their split
frequencies should be similar (barring sampling error), so we want the SDSFs to
be small (a mean SDSF < 0.1 is a rule of thumb that is often used; but, that's
not a magic number!).

**tree_length**
"""""""""""""""

Next, the ``tree_length`` section provides a summary of the total length (i.e.,
the sum of all branch lengths) of sampled trees::

    tree_length:
        n: 5200
        ess: 1667.86323370249988
        mean: 0.0110116006001833976
        median: 0.0109497159206294985
        std_dev: 0.00124402282794587181
        range: [0.00728641572310199936, 0.0161335463633999997]
        eti_95: [0.00877078878371835048, 0.0136830748419434732]
        hpdi_95: [0.00857532062927799862, 0.0134531947701339999]
        psrf: 0.999966895173833081

Some of these summaries are straightfoward, but let's flesh out some of them,
because you will see them many times in the summary file:

.. _sum_stat_defs:

``n``
    The number of posterior samples.

``ess``
    The effective sample size. Due to correlations across neighboring MCMC
    samples, the effective number of samples is often less than ``n``.

``eti_95``
    The equal-tailed 95% credible interval.
    I.e., the values associated with the 2.5 and 97.5 percentiles.

``hpdi_95``
    The highest posterior density 95% credible interval.
    The narrowest interval that includes 95% of the distribution (95% of the
    samples when trying to approximate it from a finite sample).

``psrf``
    The potential scale reduction factor.
    This measures whether independent MCMC chains are sampling overlapping
    values for a parameter (i.e., are they sampling from the same
    distribution); values close to one (e.g., < 1.2) suggest the chains
    converged and are sampling from the same distribution.

**summary_of_map_topologies**
"""""""""""""""""""""""""""""

Up next, there's a summary of the maximum *a posteriori* tree
topologies (tree models):

summary_of_map_topologies:
    - count: 280
      frequency: 0.0538461538461538491
      newick: ( ... )[& ... ]:0

In my summary there is only 1, but there might be more than one if multiple
topologies tie for the most frequently sampled.
The ``count`` is the number of times the tree was sampled, and the
``frequency`` is simply the ``count`` divided by the total number of samples
(i.e., an approximation of the posterior probability of the topology).
Lastly, the annotated MAP topology is provided in ``newick`` format;
this is the same tree we told |sumphyco| to write to the
``cyrt-map-tree.nex``, which we discuss more
:ref:`in the section below <annotated_map_tree>`.

**topologies**
""""""""""""""

Up next, the ``topologies`` section provides a list of *all* sampled topologies
sorted from most to least frequent::

    topologies:
        -
          count: 280
          frequency: 0.0538461538461538491
          cumulative_frequency: 0.0538461538461538491
          number_of_heights: 3
          heights:
              - number_of_nodes: 1
                splits:
                    - leaf_indices: [3, 4]
                n: 280
                ess: 133.894317591770431
                mean: 0.00063164050550695394
                median: 0.000625316888421000071
                std_dev: 0.000188093550720526111
                range: [0.000186275333359999992, 0.00132545817561999994]
                eti_95: [0.000295922009107825012, 0.0010075712220599998]
                hpdi_95: [0.000294101649274000024, 0.00100045882127]
              - number_of_nodes: 2
                splits:
                    - leaf_indices: [2, 3, 4]
                      node:
                          descendant_splits:
                              - [2]
                              - [3, 4]
                    - leaf_indices: [5, 6]
                n: 280
                ess: 250.986736852775152
                mean: 0.000940334872247221092
                median: 0.000939304089077999931
                std_dev: 0.000173866213213687921
                range: [0.000495016237200999977, 0.00143800553865000005]
                eti_95: [0.00063093969598985002, 0.00128013828409874984]
                hpdi_95: [0.000631647972762000012, 0.00127991331285999999]
              .
              .
              .

For each topology, the number of times it was sampled (``count``), frequency
(i.e., approximate posterior probability), and cumulative frequency is given.
Each topology is described as a list of its divergence times (``heights``).
For each divergence, the number of nodes that map to it
(``number_of_nodes``) is given, and the ``splits`` (clades) defined
by those nodes are listed.
If a split (clade) contains more than 2 leaves, then the splits that descend
from it are also listed.
Then the sampled values of the divergence time (or height)
are summarized;
:ref:`see above descriptions <sum_stat_defs>`
of some of the more cryptic summary statistics.

**heights**
"""""""""""

Next, the ``heights`` section provides a list of all sampled divergence events
(or ``heights``) in order of decreasing frequency (decreasing approximate
posterior probability)::

    heights:
        -
          number_of_nodes: 2
          splits:
              - leaf_indices: [2, 3, 4]
              - leaf_indices: [5, 6]
          count: 2586
          frequency: 0.497307692307692284
          n: 2586
          ess: 1609.86703630561283
          mean: 0.000936390624723727988
          median: 0.000927218230700499989
          std_dev: 0.000171689597604143731
          range: [0.000479753482237999993, 0.00162345554608000001]
          eti_95: [0.000633384379520000027, 0.00129245641293750015]
          hpdi_95: [0.000631576245212999983, 0.00128864682125999996]
        -
          number_of_nodes: 1
          splits:
              - leaf_indices: [2, 3, 4]
          count: 2559
          frequency: 0.492115384615384621
          n: 2559
          ess: 1165.99838502598004
          mean: 0.00101518354159117708
          median: 0.00100011380661999992
          std_dev: 0.000218746084096721499
          range: [0.000475600605589000022, 0.00193458727650999999]
          eti_95: [0.000629584362187399953, 0.00148788643064199939]
          hpdi_95: [0.000601960582020000029, 0.00144496080627999992]

For each divergence event (height), the number of nodes mapped to it
(``number_of_nodes``) is given, along with the list of splits (clades)
associated with those nodes.
The ``count`` and ``frequency`` of the divergence event is given
followed by 
the summary of the sampled values of the height (divergence time);
:ref:`see above descriptions <sum_stat_defs>`
of some of the more cryptic summary statistics.

**splits**
""""""""""

The next section, ``splits``, summarizes all of the splits that were
sampled.
Because |phyco| only estimates rooted trees, we can think of a split (branch)
and clade interchangeably, and define it by all the leaves that descend from it
(i.e., it "splits" those leaves off from the rest of the tree).
Each split (branch) on a tree will have a length and effective population size,
and the height or divergence time of a split will always refer to the height of
the branch's tipward node (i.e., the node that represents the most recent
common ancestor of the clade created by the split).

The ``splits`` section is organized into three main subsections: ``root``,
``leaves``, and ``nontrivial_splits``.
The root and all leaf branches are always present in every tree, and so
these will always have a frequency of 1.
All other splits in a tree are considered ``nontrivial_splits``,
because they are not always present in every tree.
For every split, there is a summary of the sampled values of the splits
divergence time (``height``), effective population size (``pop_size``), and
length (``length``).
For the root, the length is always zero, and for the leaves the height is
always zero.

For every sampled split that has more than two descendants, there is also a
summary of the ``nodes`` associated with the split.
We define a node by the splits that descend from it. 
For example, in the tree below, we can represent the "split" (branch) marked by
the "X" by listing the leaves that descend from it: [A, B, C].

::

                /---- E
                |
      /---------|
      |         |
      |         \---- D
    --|
      |
      |     /-------- C
      |     |
      \--X--|
            |   /---- B
            |   |
            \---|
                |
                \---- A

However, the "node" associated with "X" is more specific, because
it a list of the splits that descend from it: [[A, B], [C]].
So, for Split X, there are four possible nodes:

-   [A, B], [C]
-   [A, C], [B]
-   [B, C], [A]
-   [A], [B], [C]

For each node listed for a split in the |yaml|_
file, you will also find a summaries for its ``height``,
``length``, and ``pop_size``.
These are different from the summaries for the split, because it only
summarizes the samples that had that particular node configuration.

**number_of_heights_summary**
"""""""""""""""""""""""""""""

::

    number_of_heights_summary:
        n: 5200
        ess: 107.736499736150066
        mean: 3.83980769230769337
        median: 4
        std_dev: 0.812666578741826351
        range: [2, 6]
        eti_95: [2, 5]
        hpdi_95: [2, 5]

**number_of_heights**
"""""""""""""""""""""

::

    numbers_of_heights:
        -
          number_of_heights: 4
          count: 2527
          frequency: 0.485961538461538445
          cumulative_frequency: 0.485961538461538445
        -
          number_of_heights: 3
          count: 1494
          frequency: 0.287307692307692319
          cumulative_frequency: 0.77326923076923082
        -
          number_of_heights: 5
          count: 881
          frequency: 0.169423076923076926
          cumulative_frequency: 0.942692307692307718
        -
          number_of_heights: 2
          count: 204
          frequency: 0.0392307692307692288
          cumulative_frequency: 0.981923076923076898
        -
          number_of_heights: 6
          count: 94
          frequency: 0.0180769230769230772
          cumulative_frequency: 1


.. _annotated_map_tree:

Annotated MAP tree
^^^^^^^^^^^^^^^^^^

The ``--map-tree-out cyrt-map-tree.nex`` option above told |sumphyco| to
write the maximum *a posteriori* (MAP) tree to a file named
``cyrt-map-tree.nex``.
The MAP tree is simply the tree topology most frequently sampled by the MCMC
chains.

If you open the trees with something like FigTree, and display ``height_index``
for the node labels, you can see the shared divs (nodes with the same number).
``index_freq`` shows the posterior probability of the shared divs.
``node_freq`` shows the posterior probability of multifurcations, and
``split_freq`` shows the normal posterior probabilities you get from a standard
Bayesian phylo analysis.

The divergence times are summarized by divergence events (height indices;
index_height_mean or index_height_median)
