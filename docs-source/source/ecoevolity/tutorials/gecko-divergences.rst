|eco_logo_long|

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
|Tracer|_ (https://github.com/beast-dev/tracer/releases);
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
populations, and thus we do *not* expect any of these pairs co-diverged.
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
If you do the latter, you'll have to unzip the archive before proceeding.

Once downloaded, navigate into the example data directory using the
command line::

    cd ecoevolity-example-data

.. admonition:: Extra steps for Docker users

    If you are using the Docker container to run |eco|, follow these two extra
    steps to be able to analyze the example files you just downloaded.
    First, make sure you are **NOT** inside the Docker container before
    proceeding; you want to be at your computer's command-line console.
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

    For the remainder of the tutorial, enter (or copy and paste) any
    command-line instructions into the command line of the Docker container.
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
Based on simulations of data sets comprising loci with 100, 500, and 1000
linked characters :cite:`Oaks2018ecoevolity`, we strongly recommend that you
analyze all of your characters (including the constant ones) and violate the
assumption of unlinked characters.
In short, |Eco| performs better when you use all of the sites (especially the
constant ones) compared to reducing the data to only one variable character per
locus.
The example data sets we'll be analyzing consist of 200 loci each comprising
about 90 linked sites.

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
:ref:`click here <configfile>`.
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

This tells |eco| that, for the prior on :term:`event models<Event model>`, we want to use a
Dirichlet process with its concentration parameter fixed to 1.414216.
The concentration parameter allows us to control how much sharing of divergence
times we expect *a priori*.
For a gentle introduction to the Dirichlet process, check out
:ref:`this section of the background page <prior_on_divergence_models>`
or
`this blog post <http://phyletica.org/dirichlet-process/>`_.

One of the |eco| tools, called |dpprobs|, was designed to help you get a feel
for the Dirichlet process and choose a prior on the concentration parameter.

When you type::

    dpprobs -h

on the command line, you should see the help menu for |dpprobs| written to your
screen.
Let's explore the current setting in the configuration file::

    dpprobs -p concentration 1.414216 3

This command tells the program (1) the parameter (``-p``) we are providing
is the concentration parameter, (2) the value of that parameter is 1.414216,
and (3) there are three comparisons.
You should get output that looks similar to::

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
process with the specified concentration parameter and 3 comparisons is 2.
It also shows us Monte Carlo approximations for the prior probability for each
possible number of divergence events.

Now, let's find out what value of the concentration parameter corresponds
with a prior mean number of divergences events of 2.5.
We'll also request 1 million simulations under the Dirichlet process to get
more precise Monte Carlo estimates of the prior probabilities of the possible
numbers of divergence events::

    dpprobs -p mean -n 1000000 2.5 3

Our output now looks something like::

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
expectation of 2.5 divergence events (when there are 3 comparisons).
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

The output now shows us which gamma distribution has a mean concentration that
corresponds with a Dirichlet process with an expected number of divergence
events of 2.5.
There is no "correct" or "best" value or prior to use for the concentration
parameter.
It depends on your system and prior expectations.
Often, when testing for shared divergences, I like to choose a value or prior
for the concentration parameter that puts 50% of the prior probability on the
maximum number of divergence events (i.e., the divergence model with no shared
divergences).
That way, if the results support shared divergences, I can be more confident
the data are driving the result.
But, this is just my personal preference.

Okay, let's edit our config to increase the concentration parameter to favor
more divergence events *a priori*.
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

Under this new setting, what is the prior expectation for the number of
divergence events, and what is the prior probability for the model with 3
divergence events?
Use |cdpprobs| to find out!

.. note::

    We recommend that you analyze your data under multiple settings for the
    concentration parameter, so you can assess how sensitive the results
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
``comparison`` settings, this tells us we are fixing the mutation rate to
1.0 for all three pairs of gecko populations.
As a result, our time units are expected substitutions per site.

Let's assume that, based on prior information, we are quite confident these
geckos have experienced fewer than 1 substitution per every 10 characters since
the populations diverged.
Let's change the prior to a gamma distribution with a shape of 2.0 and a
scale of 0.005.
Change::


    event_time_prior:
        exponential_distribution:
            rate: 10.0

To::

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
    priors on the divergence times, so you can assess how sensitive the
    results are to your prior assumptions.
    The prior on the divergence times can have a large affect on the marginal
    likelihoods of the divergence models, and thus the model posterior
    probabilities
    :cite:`Oaks2012,Oaks2014reply,Oaks2018marginal`.
    While |eco| tends to be much less sensitive to this than the ABC methods,
    it is still good to examine the affect of the divergence time prior.


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
characters with more than two states, or (2) recode these sites
as biallelic.
The latter |eco| will do for us if we specify the ``--relax-triallelic-sites``
option in our command.
When we use this option, |eco| will consider the first state in a column as
``0``, and any other state found in the column as ``1``.
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
data for at least one of the two populations of the Babuyan Claro-Calayan pair.
Such sites can be common in RADseq loci, because most assemblers seem to
enforce thresholds on missing data at the locus level.
Rather than removing these sites ourselves, we can tell |eco| to ignore them::

    ecoevolity --relax-triallelic-sites --relax-missing-sites ecoevolity-config.yml

Now, you should be running!

Checking the pre-MCMC screen output
-----------------------------------

While the analysis is running, let's scroll up in our console and look at the
information |eco| reported before the MCMC chain started sampling.
This output goes by very quickly, but it's **very** important to read it
over to make sure |eco| is doing what you think it's doing.
At the very top, you will see something like::

    ======================================================================
                                  Ecoevolity
                      Estimating evolutionary coevality
          Version 0.2.0 (master c59888c: 2018-05-27T13:26:11-05:00)
    ======================================================================
    
    Seed: 713687808
    Using data in order to sample from the posterior distribution...
    Config path: ecoevolity-config.yml

This gives us detailed information about the version of |eco| we are using
(right down to the specific commit SHA).
It also tells us what number was used to seed the random number generator; we
would need to provide |eco| this seed number in order to reproduce our results.
The output is also telling us that |eco| is using the data to sample from the
posterior, as opposed to ignoring the data to sample from the prior.

Next, |eco| reports a fully specified configuration file to the screen.
It's always good to look this over to make sure |eco| is configured as you
expect.
After the config output, you should see a series of warning messages alerting
you to the fact that |eco| is recoding some polyallelic sites to biallelic, and
also ignoring some sites that lack data from at least one population.
Without using the ``--relax-triallelic-sites`` and  ``--relax-missing-sites``
options these warnings would have been errors that stopped |eco| from running.
Always look over these warnings to make sure |eco| is not removing or
recoding anything you didn't expect.

Next, you will see a summary of the data from your comparisons::

    ----------------------------------------------------------------------
    Summary of data from 3 comparisons:
        Summary of data from 'nexus-files/BabuyanClaro-Calayan.nex':
            Genotypes: diploid
            Markers are dominant? false
            Number of populations: 2
            Number of sites: 27340
            Number of variable sites: 102
            Number of patterns: 71
            Patterns folded? true
            Population label (max # of alleles sampled):
                BabuyanClaro0 (10)
                Calayan0 (10)
        Summary of data from 'nexus-files/MaestreDeCampo-Masbate.nex':
            Genotypes: diploid
            Markers are dominant? false
            Number of populations: 2
            Number of sites: 27311
            Number of variable sites: 257
            Number of patterns: 51
            Patterns folded? true
            Population label (max # of alleles sampled):
                MaestreDeCampo1 (6)
                Masbate1 (6)
        Summary of data from 'nexus-files/Dalupiri-CamiguinNorte.nex':
            Genotypes: diploid
            Markers are dominant? false
            Number of populations: 2
            Number of sites: 27285
            Number of variable sites: 112
            Number of patterns: 76
            Patterns folded? true
            Population label (max # of alleles sampled):
                CamiguinNorte2 (10)
                Dalupiri2 (10)
    ----------------------------------------------------------------------

It is **very important** to look over this summary to make sure |eco|
"sees" your data the way you expect it should.
Make sure the number of sites, the number of populations, and the number of
alleles per population match your expectations!

Next, |eco| reports a summary of the state of the MCMC chain to the screen for
every 10th sample it's logging to an output file.
After the chain finishes, |eco| also reports stats associated with all of the
MCMC operators; this information is also logged to a file.


The Output Files
================

After the analysis finishes, type ``ls`` at the command line::

    ls

You should see two log files that were created by |eco|:

#. ``ecoevolity-config-operator-run-1.log``
#. ``ecoevolity-config-state-run-1.log``

If you open the ``operator`` log with a plain text editor, you'll see this file
contains information about the MCMC operators.
If you have convergence or mixing issues, this information can be useful for
troubleshooting.

If you open the ``state`` file, you'll see it contains the MCMC samples.

Summarizing and Plotting the Results
====================================

Before we start summarizing the results, let's run a second chain so we can
assess convergence and boost our posterior sample size::

    ecoevolity --relax-triallelic-sites --relax-missing-sites ecoevolity-config.yml

Now's a good time to go grab a coffee, while you wait for the second chain to
finish running.
After the chain finishes, if you use ``ls`` you should see another ``operator``
and ``state`` log file:

#. ``ecoevolity-config-operator-run-2.log``
#. ``ecoevolity-config-state-run-2.log``

Assessing mixing and convergence
--------------------------------

Now, let's use the |pyco-sumchains| tool of the |pyco| packages to help assess
the convergence of our chains and choose what number of samples we want to remove
as "burn in"::

    pyco-sumchains -s 100 ecoevolity-config-state-run-1.log ecoevolity-config-state-run-2.log

|pyco-sumchains| removes an increasing number of samples
from the beginning of the chains (i.e., burn in) and reports the effective
sample size (ESS) and potential scale reduction factor (PSRF) for every
continuous parameter and log likelihood score.
The ``-s 100`` option tells |pyco-sumchains| to remove 100 samples at a time.
After the table of ESS and PSRF values, it also reports a summary that can
be useful::

    Continuous parameter with lowest mean ESS: pop_size_Dalupiri2
    Burnin value that maximized ESS of pop_size_Dalupiri2: 100 samples (ESS = 544.00)

This tells us which parameter had the lowest mean ESS value
across all the burn-in values.
It also tells us what burn-in value was associated with the highest ESS for
this parameter.

If you want the ESS and PSRF table in a text file that you can open with Excel,
simply use ``>`` to redirect the standard output::

    pyco-sumchains -s 100 ecoevolity-config-state-run-1.log ecoevolity-config-state-run-2.log > pyco-sumchains-table.txt

Once |pyco-sumchains| finishes you can open the ``pyco-sumchains-table.txt``
file with Excel to see how the ESS and PSRF of each parameter changes as the
burn in increases.

If you have |Tracer|_ installed, you can use it to look at the convergence and
mixing behavior of the chains.

From my |pyco-sumchains| output and inspecting the chains in |Tracer|_,
ignoring the fist 101 samples as burn-in seems conservative.
That's the burn-in value I will use below when summarizing and plotting the
results.

Summarizing divergence-model posterior probabilities
----------------------------------------------------

Next, let's use the |sumco| tool to calculate the posterior probabilities and
Bayes factors for divergence models::

    sumcoevolity -b 101 -c ecoevolity-config.yml -n 1000000 ecoevolity-config-state-run-1.log ecoevolity-config-state-run-2.log

Let's break down what these arguments mean:

*   ``-b 101`` tells |sumco| to ignore the first 101 samples in our log files
*   ``-c ecoevolity-config.yml`` tells |sumco| where our config file is.

    |sumco| will use the config file to run simulations under the Dirichlet
    process in order to get prior probabilities of divergence models for
    calculating Bayes factors.
*   ``-n 1000000`` tells |sumco| to use 1 million simulations under the
    Dirichlet process to approximate the prior probabilities

|sumco| will produce two files:

#.  ``sumcoevolity-results-model.txt``
#.  ``sumcoevolity-results-nevents.txt``

Let's open the ``sumcoevolity-results-model.txt`` first.
The content should looks something like::

    model	post_prob	cumulative_post_prob	prior_prob	bf
    0,1,2	>0.999545	1	0.50013	>2197.86

This is showing us a list of models that were sampled by the chains, ranked by their
posterior probability.
In my case, there was only a single model sampled: ``0,1,2``.
The notation is a list of divergence-event indices for each of our population pairs.
In my case, all three pairs have a unique index, and so all three diverged at
their own event (the 3-divergence model).
The output also includes the cumulative posterior probability so you can easily
evaluate credible sets.
It also shows the prior probability and Bayes factor (``bf``).
The Bayes factor is the factor by which the posterior odds of the model has
changed from our prior odds.
In this case, both the posterior and prior odds are comparing the probability
of the model in the row versus all other possible models.

Now, let's open the ``sumcoevolity-results-nevents.txt``::

    number_of_events	post_prob	cumulative_post_prob	prior_prob	bf
    3	>0.999545	1	0.50013	>2197.86
    1	<0.000454545	1	0.078258	<0.00535618
    2	<0.000454545	1	0.421612	<0.000623851

This content is similar to the previous file, but now we are looking at the
number of divergence events, rather than specific models.
Each row is one of the possible numbers of events, ranked by their posterior
probability.
Just like the last file, the cumulative posterior probability, prior
probability, and Bayes factor is given for each number of divergence events.


Plotting posterior probabilities of the number of events
--------------------------------------------------------

Now, we can use the |pyco-sumevents| tool of the |pyco| package to plot the
information in the ``sumcoevolity-results-nevents.txt`` file::

    pyco-sumevents sumcoevolity-results-nevents.txt

This should create a plot (in a variety of file formats) that looks something
like the following.

.. figure:: /_static/pycoevolity-nevents.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: Pycoevolity nevent plot

This is just a bar plot representation of the information in the
``sumcoevolity-results-nevents.txt``.
If you checkout the help menu using ``pyco-sumevents -h``,
you'll see there are several options to customize your plot.
Also, |pyco-sumevents| spits out the R script it used to create the
plot, if you want to customize it further.


Plotting marginal divergence times
----------------------------------

Now let's plot the marginal posterior distributions for the divergences
times of our three pairs of gecko populations::

    pyco-sumtimes -b 101 -z ecoevolity-config-state-run-1.log ecoevolity-config-state-run-2.log

This should create a plot (in a variety of file formats) that looks something
like the following.

.. figure:: /_static/pycoevolity-times.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: Pycoevolity time plot

Again, ``pyco-sumtimes -h`` will show you options for customizing your plot,
and |pyco-sumtimes| also outputs the R script that generated the plot.

Sampling the prior distribution
-------------------------------

With |eco|, it is very easy to run analyses that ignore the data, and
thus sample from the joint prior distribution.
All we need to do is add the ``--ignore-data`` option to our command from
above.
We'll also use the ``--prefix`` option, which will add a prefix to all the
output files, allowing us to keep track of the output from sampling the prior::

    ecoevolity --ignore-data --prefix "prior-" --relax-triallelic-sites --relax-missing-sites ecoevolity-config.yml

This should output two files:

*   prior-ecoevolity-config-operator-run-1.log
*   prior-ecoevolity-config-state-run-1.log

It should be *really* fast, making it very quick and easy to thoroughly
sample from the joint prior via MCMC.
This is a great way to "sanity check" your analyses;
quickly generate samples from the prior and make sure they match your
assumptions.
Go ahead and run the above command again to run another independent chain to
sample the prior::

    ecoevolity --ignore-data --prefix "prior-" --relax-triallelic-sites --relax-missing-sites ecoevolity-config.yml

Now, you can use |pyco-sumchains| and |Tracer|_ to make sure both MCMC chains
that sampled from the prior converged and mixed well, and to determine an
appropriate number of samples to remove as burn-in::

    pyco-sumchains -s 100 prior-ecoevolity-config-state-run-?.log

For my output, a burn-in of 101 seemed more than sufficient.
Now, to summarize the prior samples run::

    sumcoevolity -b 101 --prefix "prior-" -c ecoevolity-config.yml -f prior-ecoevolity-config-state-run-?.log

This should create 2 files:

*   ``prior-sumcoevolity-results-nevents.txt``
*   ``prior-sumcoevolity-results-model.txt``

Open these up and check to see if the "posterior" probabilities are similar to
the prior probabilities.
I put "posterior" in quotes, because the numbers in these columns are actually
the prior probabilities approximated by the MCMC chains that ignored the data.


Comparison with ABC
===================

Here are my results from |eco| after running 5 independent MCMC chains:

.. figure:: /_static/pycoevolity-nevents-pretty.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: Pycoevolity nevent plot from 5 chains

.. figure:: /_static/pycoevolity-times-pretty.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: Pycoevolity time plot from 5 chains

The average runtime for the chains was 7.6 minutes, so about :math:`7.6 \times
5 = 38` minutes of total computing time.
For comparison, below are results I obtained from analyzing the same data with
the ABC method ``dpp-msbayes`` using almost identical priors and simulating
500,000 samples from the prior.
This analysis took 6 days to run with simulations spread across 8 processors,
so about :math:`6 \times 8 = 48` days of computing time.

.. figure:: /_static/dppmsbayes-nevents.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: dpp-msbayes nevent plot

.. figure:: /_static/dppmsbayes-times.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: dpp-msbayes time plot
