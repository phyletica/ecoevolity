
.. _gecko_div_tutorial:

*****************
Gecko Divergences
*****************

Objective
=========

Our goal is to learn how to set up, run, and interpret an |eco| analysis in
order to test models of shared divergence times among taxa.

Background Information
======================

.. include:: /snippets/brief-background.rst

For more detailed background information about the method implemented in |eco|,
please :ref:`see here <background>`.

Software Used in this Activity
==============================

For this activity, we will be using |eco| and |pyco|.
If you haven't already done so, please follow the
:ref:`instructions for installing these packages <installation>`.
Either installing |eco| and |pyco| directly or using the Docker image will work
for this tutorial.

While not required, you may want to install the program
|Tracer|_ (http://tree.bio.ed.ac.uk/software/tracer/);
it's a nice tool for visualizing the mixing and convergence behavior of Markov
chain Monte Carlo (MCMC) analyses.

The Data
========

We will be analyzing "RADseq" data from three pairs of gecko populations from
the Philippine Islands.
Each pair of populations inhabit two oceanic islands in the Philippines:

#.  Babuyan Claro and Calayan
#.  Maestre De Campo and Masbate
#.  Dalupiri and Camiguin Norte

Based on previous data :cite:`Siler2012,Siler2014kikuchii`,
each of these pairs diverged very recently.
However, because these islands were never connected, each pair's distribution
across both islands must have been due to over-water dispersal.
Such rare dispersal events were likely idiosyncratic to each pair of island
populations, and thus we do *not* expect that any of these pairs co-diverged.
For each pair of populations, the data comprise 200 RADseq loci randomly
sub-sampled from the loci analyzed by
Oaks :cite:`Oaks2018ecoevolity`.

If you have not already done so, please download the data.
If you have |git|_ installed, and you are in your command-line console,
you can get the data via::

    git clone https://github.com/phyletica/ecoevolity-example-data.git

If you prefer not to use |git|_, you can download an archive of the example
data by
`clicking here <https://github.com/phyletica/ecoevolity-example-data/archive/master.zip>`_.
If you do the latter, you'll have to unzip the archive.

Once downloaded, navigate into the example data directory using the
command line::

    cd ecoevolity-example-data

.. admonition:: Extra steps for Docker users

    If you are using the Docker container to run |eco|, follow these two extra
    steps to be able to analyze the example files you just downloaded.
    After you've used ``cd`` to navigate into the example data directory,
    execute the following command to enter the Docker container::
    
        docker run -v "$(pwd)":/portal -it phyletica/ecoevolity-docker bash
    
    .. note::
    
        Depending on your system and how Docker is configured, you may need to use
        ``sudo`` to run Docker commands. If you received a "permission denied"
        message when you ran the command above, try::
        
            sudo docker run -v "$(pwd)":/portal -it phyletica/ecoevolity-docker bash
    
    Then, once inside the Docker container, type::
    
        cd portal

    For the remainder of the tutorial, enter (or copy and paste) any command
    line instructions at the command line of the Docker container.
    However, you can view and edit the input and output files on your computer
    (outside the container) with your plain text editor of choice.

Inside the example data directory, you should find a configuration file
``ecoevolity-config.yml`` and a ``nexus-files`` directory containing the
nexus-formatted RADseq data for each of our three pairs of gecko populations.
For detailed information about how to format nexus files for |eco|,
:ref:`click here <data>`.

Biallelic characters
--------------------

|Eco| assumes all of your characters have at most two states (biallelic).
If we provide |eco| with nucleotide data, it will automatically recode the data
as biallelic by considering the first character it finds in each column as
state ``0``, and if it finds a second state it considers it state ``1``.
If it finds a third state in a column, it will report an error message and
quit.
For characters (columns) with more than two states, you have two options:

#.  remove these columns from your matrix, or
#.  recode them as biallelic.

|Eco| will do the latter for you (more on this in a bit).

Linked characters
-----------------

|Eco| also assumes all your characters are unlinked.
Based on simulations of loci of 100, 500, and 1000 linked characters
:cite:`Oaks2018ecoevolity`, we strongly recommend that you analyze all of your
characters (including the constant ones) and violate the assumption of unlinked
characters.
In short, |Eco| performs better when you use all of the sites (especially the
constant ones) compared to reducing the data to only one variable character per
locus.
The example gecko data sets consist of 200 loci each comprising about 90 linked
sites.

Orthology
---------

Within each pair, |eco| assumes each column represents an orthologous character
across your samples from both populations.
However, the characters do not need to be orthologous among your different
pairs.
In other words, you can sample different loci from different pairs of
populations.
Thus, if you are assembling your loci from raw sequence reads, it might be
worth assembling each population pair separately to maximize the number of loci
you get for each pair.

Running |eco|
=============

Setting up the configuration file
---------------------------------

Once we have our nexus files ready, the next step in an |eco| analysis is
setting up the configuration file.
|Eco| requires a |yaml|_-formatted configuration file.
|yaml|_ is a human-friendly data standard that allows you to provide |eco| the
information it needs in a format that is easy for you to read and edit.

.. note::

    The website |yamllint| is a nice tool for debugging |yaml|_ syntax.
    You can copy and paste your config file there to check if you're using
    valid |yaml| syntax.


For detailed information about all the settings that can be included in an
|eco| config file,
:ref:<`click here <configfile>`.
For the purposes of this tutorial, we'll edit a few elements of the example
configuration file we downloaded.

Choosing a prior on the divergence models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Open the ``ecoevolity-config.yml`` file in the example data directory
using a plain text editor.
Near the top, you will see a section that looks like::

    event_model_prior:
        dirichlet_process:
            parameters:
                concentration:
                    value:      1.414216
                    estimate:   false

This tells |eco| that for the prior on :term:`event models`, we want to use a
Dirichlet process with its concentration parameter fixed to 1.414216.
The concentration parameter allows us to control how much sharing of divergence
times we expect *a priori*.
For an introduction to the Dirichlet process, check out
:ref:`this section of the background page <prior_on_divergence_models>`
or
`this blog post <http://phyletica.org/dirichlet-process/>`_.

One of the |eco| tools, called |dpprobs|, was designed to help you get a feel
for the Dirichlet process and choose a prior on the concentration parameter.

When you type::

    dpprobs -h

on the command line, you should see the help menu for |dpprobs| written to your
screen.
Let's explore the current setting in the configuration, which specifies a
concentration parameter of 1.414216::

    dpprobs -p concentration 1.414216 3

This command tells the program that (1) the parameter (``-p``) we are providing
is the concentration parameter, (2) the value of that parameter is 1.414216,
and (3) there are three comparisons.
You should get output that looks like::

    ======================================================================
                                   DPprobs
                  Simulating Dirichlet process probabilities
    
                                   Part of:
                                  Ecoevolity
          Version 0.2.0 (master 4f8811a: 2018-04-30T13:55:12-05:00)
    ======================================================================
    
    Seed = 117483248
    Number of samples = 100000
    Number of elements = 3
    Concentration = 1.41422
    Mean number of categories = 2
    
    Estimated probabilities of the number of categories:
    ----------------------------------------------------------------------
    p(ncats = 1) = 0.24394      (n = 1)
    p(ncats = 2) = 0.51361      (n = 3)
    p(ncats = 3) = 0.24245      (n = 1)
    ----------------------------------------------------------------------

    Runtime: 0 seconds.

This tells us that the prior mean number of divergence events under a Dirichlet
process with the specified concentration parameter and 3 elements is 2.
It also shows us Monte Carlo esimates for the prior probability for each
possible number of divergence events.

Now, let's find out what value of the concentration parameter corresponds
with a prior mean number of divergences events of 2.5.
We'll also request 1 million simulations under the Dirichlet process to
get more precise Monte Carlo estimates of the prior probabilities of
all possible numbers of divergence events::

    dpprobs -p mean -n 1000000 2.5 3

Our output now looks like::

    ======================================================================
                                   DPprobs
                  Simulating Dirichlet process probabilities
    
                                   Part of:
                                  Ecoevolity
          Version 0.2.0 (master 4f8811a: 2018-04-30T13:55:12-05:00)
    ======================================================================
    
    Seed = 268086367
    Number of samples = 1000000
    Number of elements = 3
    Concentration = 4.37229
    Mean number of categories = 2.5
    
    Estimated probabilities of the number of categories:
    ----------------------------------------------------------------------
    p(ncats = 1) = 0.05873      (n = 1)
    p(ncats = 2) = 0.383424     (n = 3)
    p(ncats = 3) = 0.557846     (n = 1)
    ----------------------------------------------------------------------
    
    Runtime: 0 seconds.


This tells us a concentration parameter of about 4.37 corresponds with a prior
expectation of 2.5 divergence time events (when there are 3 comparisons).
The output also shows us that the prior probability of the model with three
divergences is almost 56\%.

If, instead of fixing the concentration parameter, we prefer to place a
gamma-distributed prior on it, we can specify the shape parameter we would
like to use for that gamma prior::

    dpprobs -p mean -n 1000000 --shape 2.0 2.5 3

Output::

    ======================================================================
                                   DPprobs
                  Simulating Dirichlet process probabilities
    
                                   Part of:
                                  Ecoevolity
          Version 0.2.0 (master 4f8811a: 2018-04-30T13:55:12-05:00)
    ======================================================================
    
    Seed = 1423802726
    Number of samples = 1000000
    Number of elements = 3
    Concentration ~ gamma(shape = 2, scale = 2.18615)
    Mean number of categories given mean concentration = 2.5
    
    Estimated probabilities of the number of categories:
    ----------------------------------------------------------------------
    p(ncats = 1) = 0.123055     (n = 1)
    p(ncats = 2) = 0.392924     (n = 3)
    p(ncats = 3) = 0.484021     (n = 1)
    ----------------------------------------------------------------------
    
    Runtime: 0 seconds.

There is no "correct" or "best" value or prior to use for the concentration
parameter.
It depends on your system and prior expectations.
Often, when testing for shared diveregences, I like to choose a value or prior
for the concentration parameter that puts 50% of the prior probability on the
maximum number of divergence events (i.e., the divergence model with no shared
divergences).
That way, if the results support shared divergences, I can be more confident
that the data are driving the result (however, there's no guarantee; lack of
information and a prior on divergence times that puts a lot of density in
regions with low likelihood can also drive such a result).
But, that is just a personal preference (i.e., there is no fundamental
mathematical justification for it).

Okay, let's increase teh concentration parameter to favor more divergence
events *a priori*.
Change::

    event_model_prior:
        dirichlet_process:
            parameters:
                concentration:
                    value:      1.414216
                    estimate:   false

To::

    event_model_prior:
        dirichlet_process:
            parameters:
                concentration:
                    value:      3.564
                    estimate:   false

Under this new setting, what is the prior expectation for the number
of divergence event, and what is the prior probability for the model
with 3 divergence events?
Use |dpprobs| to find out!

.. note::

    We recommend that you analyze your data under multiple settings for the
    concentration parameter, so that you can assess how sensitive the results
    are to your prior assumptions.


Choosing a prior on divergence times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, let's change the prior on the divergence times.
Currently, the config file specifies an exponential distribution with a rate of
10::

    event_time_prior:
        exponential_distribution:
            rate: 10.0

This corresponds to a prior mean divergence time of 1/rate = 1/10.0 = 0.1.
What are the units of time?
To determine this, we need to look at the settings for the mutation rate.
If we look in the ``global_comparison_settings`` section, we see::

            mutation_rate:
                value: 1.0
                estimate: false

Given that this global comparison setting is not overridden within any of the
``comparison`` settings, this tells us that we are fixing the mutation rate of
all three pairs of gecko populations to 1.0.
As a result, our time units are expected substitutions per site.

Let's assume that, based on prior information, we are quite confident that
these geckos have experience fewer than 1 substitution per
every 10 characters since the populations diverged.
Let's change the prior to a gamma distribution with a shape of 2.0 and a
scale of 0.005::

    event_time_prior:
        gamma_distribution:
            shape: 2.0
            scale: 0.005

This corresponds with a prior mean time of divergence of
:math:`\textrm{shape} \times \textrm{scale} = 2(0.005) = 0.01`
substitutions per site, and variance of 
:math:`\textrm{shape} \times \textrm{scale}^2 = 2(0.005^2) = 0.00005`.

.. note::

    We strongly recommend that you analyze your data under multiple different
    priors on the divergence times, so that you can assess how sensitive the
    results are to your prior assumptions.
    The prior on the divergence times can have a large affect on the marginal
    likelihoods of the divergence models, and thus the model posterior
    probabilities
    :cite:`Oaks2012,Oaks2014reply,Oaks2018marginal`.
    While |eco| is much less sensitive to this than the ABC methods, it is
    still good to examine the affect of the divergence time prior.


Running the analysis
--------------------

Now that we have our data files and config file ready, let's go ahead
and try running an analysis with |eco|::

    ecoevolity ecoevolity-config.yml

Oops, we got an error::

    #######################################################################
    ###############################  ERROR  ###############################
    1 sites from the alignment in:
        'nexus-files/BabuyanClaro-Calayan.nex'
    have more than two character states.
    #######################################################################

This error message is telling us that we have a site (column) in one of our
character matrices where there are more than two states represented across our
samples.
As discussed above, at this point, we have 2 options: (1) remove any
characters that have more than two states, or (2) recode these sites
as biallelic.
The latter, |eco| will do for us if we specify the ``--relax-triallelic-sites``
option in our command.
When we use this option, |eco| will consider the first state in a column as
``0``, and any other state found in that column as ``1``.
In all the data sets we've analyzed so far, there has been no discernible
difference in the results between removing or recoding the triallelic sites.
However, if you have such sites in your data, we recommend trying both options
to see if it affects the results.
Ok, let's go ahead and let |eco| recode our characters with more than two
states as binary::


    ecoevolity --relax-triallelic-sites ecoevolity-config.yml

Shoot, another error!::

    #######################################################################
    ###############################  ERROR  ###############################
    25 sites from the alignment in:
        'nexus-files/BabuyanClaro-Calayan.nex'
    have no data for at least one population.
    #######################################################################

This error message is informing us that for some of our characters, we have no
data for at least one of the two populations.
Such sites can be common in RADseq loci, because most assemblers seem to
enforce thresholds on missing data at the locus level.
Rather than removing these sites ourselves, we can tell |eco| to ignore them::

    ecoevolity --relax-triallelic-sites --relax-missing-sites ecoevolity-config.yml

Now, you should be running!

While the analysis is running, let's scroll up in our console and look at the
information |eco| reported before the MCMC chain started sampling.
This output goes by very quickly, but it's **very** important to read it
over to make sure that |eco| is doing what you think it's doing.

The output
==========

Assess stationarity with pyco-sumchains and |Tracer|

Usually run lots of chains to assess convergence

Run sumcoevolity (needs config)

Run pyco-sumtimes

Run pyco-sumevents (uses output of sumcoevolity)







Plotting the results
====================

