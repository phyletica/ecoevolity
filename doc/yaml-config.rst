.. _configfile:

######################
The Configuration File
######################

Most of the settings for an |eco|_ analysis you specify in a configuration
file.
We use a |yaml|_ configuration file.
|yaml|_ is an "anti-markup" language and data standard that is much more human
friendly than alternatives like XML.
|yamllint| is a nice tool for debugging |yaml|_ syntax.

The simplest possible configuration file for |eco|_ would look something like::

    ---
    comparisons:
    - comparison:
        path: "alignments/G-crombota-rossi-BabuyanClaro-Calayan.nex"
    - comparison:
        path: "alignments/G-mindorensis-mindorensis-Lubang-Luzon.nex"
    - comparison:
        path: "alignments/G-mindorensis-mindorensis-MaestreDeCampo-Masbate.nex"
    - comparison:
        path: "alignments/G-sp_a-sp_b-Dalupiri-CamiguinNorte.nex"

This specifies the path to the nexus file for each of your population pairs.
The defaults for all other settings would be used, which I **strongly
discourage**.

Now, let's look at an example that is more thoroughly specified::

    ---
    # This is a comment
    event_model_prior:
        dirichlet_process:
            parameters:
                concentration:
                    value: 7.5
                    estimate: false 
    
    event_time_prior:
        exponential_distribution:
            rate: 100.0
    
    mcmc_settings:
        chain_length: 75000
        sample_frequency: 50
    
    operator_settings:
        auto_optimize: true
        auto_optimize_delay: 1000
        operators:
            ConcentrationScaler:
                weight: 3.0
                scale: 0.1
            ModelOperator:
                weight: 10.0
                number_of_auxiliary_categories: 4
            TimeSizeRateMixer:
                weight: 5.0
                scale: 0.02
            TimeSizeRateScaler:
                weight: 0.0
                scale: 0.02
            EventTimeScaler:
                weight: 1.0
                scale: 0.02
    
    global_comparison_settings:
        ploidy: 2
        genotypes_are_diploid: true
        markers_are_dominant: false
        population_name_delimiter: " "
        population_name_is_prefix: false
        constant_sites_removed: false
        equal_population_sizes: false
        parameters:
            population_size:
                estimate: true
                prior:
                    gamma_distribution:
                        shape: 4.0
                        scale: 0.001
            root_relative_population_size:
                estimate: true
                prior:
                    gamma_distribution:
                        shape: 100.0
                        scale: 0.01
                        offset: 0.0
            freq_1:
                value: 0.5
                estimate: false
    
            mutation_rate:
                value: 1.0
                estimate: false
        operators:
            RootPopulationSizeScaler:
                weight: 1.0
                scale: 0.05
            LeafPopulationSizeScaler:
                weight: 1.0
                scale: 0.05
            TimeRootSizeMixer:
                weight: 3.0
                scale: 0.05
    
    comparisons:
    - comparison:
        path: "../alignments/G-crombota-rossi-BabuyanClaro-Calayan.nex"
    - comparison:
        path: "../alignments/G-mindorensis-mindorensis-Lubang-Luzon.nex"
    - comparison:
        path: "../alignments/G-mindorensis-mindorensis-MaestreDeCampo-Masbate.nex"
    - comparison:
        path: "../alignments/G-sp_a-sp_b-Dalupiri-CamiguinNorte.nex"


Everything is hierarchically nested by the indent spacing.
For example,
``event_model_prior``, ``event_time_prior``, ``mcmc_settings``,
``operator_settings``, ``global_comparison_settings``, and ``comparison``
are all at the highest level of the hierarchy, and have various settings nested
within them.
Across a given level of the hierarchy, order does not matter. E.g., the
top-level groups of settings listed above can be arranged in any order.
Anything proceeded by a '#' is a comment that is ignored by |eco|_.

*****************
event_model_prior
*****************

The ``event_model_prior`` sets up the prior probabilities for all the different
ways we can cluster the comparisons together (or not).
Thus, the term "event model" is being used to refer to each of these
possibilities.
Currently, ``dirichlet_process`` (aka "DPP") is the only option, and it has a
single setting: the ``concentration`` parameter.

The settings::

    event_model_prior:
        dirichlet_process:
            parameters:
                concentration:
                    value: 7.5
                    estimate: false 

specify that the concentration parameter of the Dirichlet process should be
fixed to a value of 7.5.
Alternatively, you can put a gamma-distributed prior on the concentration
parameter, for example::

    event_model_prior:
        dirichlet_process:
            parameters:
                concentration:
                    value: 7.5
                    estimate: true 
                    prior:
                        gamma_distribution:
                            shape: 2.0
                            scale: 3.25

will allow the concentration parameter to be estimated.
Generally, if you have a large number of comparisons (say 6 or more), it can be
helpful to allow the concentration parameter to vary.
NOTE: the prior on the concentration parameter must be a gamma distribution.

If you've installed |eco|_, you should have a command line tool called
``dpprobs`` that can help you choose a value for the concentration parameter.
Typing::

  $ dpprobs -h

on the command line should provide the help menu, but the basic usage is::

  $ dpprobs -p concentration 7.5 4

which requests the prior probabilities for the number of possible divergence
events when the concentration parameter is 7.5 and there are 4 comparisons.

I often set the concentration so that 50% of the prior probability is
on the maximum number of events/categories.
The idea is that by placing most of the prior probability on the model with no
shared divergences (or expansions), if my results indicate a shared event, I
can be more confident that he data are driving that result.
But, that is just my preference (i.e., there is no fundamental mathematical
justification it).


****************
event_time_prior
****************

The ``event_time_prior`` specifies the prior on the divergence times.
For example::

    event_time_prior:
        exponential_distribution:
            rate: 100.0

specifies an exponential distribution with a rate of 100.0 (thus
the mean of the exponential prior is 1.0/100.0 = 0.01).
A gamma or uniform distribution can also be used.

If the mutation rate of one or more comparison is set to 1.0, then units of
time are expected number of substitutions per site.
If actual mutation rates are specified for your comparisons, then the
units of time will be in whatever unit the rates are in.
For example, if you give mutation rates in substitutions per site per million
years, time will be in units of millions of years.


*************
mcmc_settings
*************

The ``mcmc_settings`` simply specify how long to run the MCMC
chain, and how often to sample. For example::

    mcmc_settings:
        chain_length: 75000
        sample_frequency: 50

tells |eco|_ to run the chain for 75,000 generations, recording a sample every
50th generation (75000/50 = 1500 total samples).
This is a good starting point. For most datasets I have analyzed so far,
this has been sufficient. If the chain is having mixing problems, then you
can increase.
I recommend running a several independent chains (analyses) to confirm they are
converging.


*****************
operator_settings
*****************

The ``operator_settings`` control the behavior of the "global" MCMC operators
that update the values of the model's parameters.
As opposed to "global," there are also MCMC operators that only operate on the
parameters of a specific comparison; these are discussed further below.

::

    operator_settings:
        auto_optimize: true
        auto_optimize_delay: 1000

This specifies that the MCMC operators should automatically adjust their tuning
parameters to try an optimize mixing. The ``auto_optimize_delay: 1000`` tells
the operators to wait until they have been used 1000 times before they start
auto-tuning (this gives them data on their acceptance rate).
Generally, auto optimization should always be used, and a delay of 1000 seems
to work well.

::

    operator_settings:
        operators:
            ConcentrationScaler:
                weight: 3.0
                scale: 0.1
            ModelOperator:
                weight: 10.0
                number_of_auxiliary_categories: 4
            TimeSizeRateMixer:
                weight: 5.0
                scale: 0.02
            TimeSizeRateScaler:
                weight: 0.0
                scale: 0.02
            EventTimeScaler:
                weight: 1.0
                scale: 0.02

These settings control the "global" MCMC operator.  The weights are relative
and control how often each operator is used.
For example, an operator with ``weight: 2`` will be used twice as often (on
average) than an operator with ``weight: 1``.
``Mixer`` and ``Scaler`` operators have a ``scale`` parameter, which
controls how large of changes it will propose for the value of 
model parameters.
If auto optimization is turned on, these are starting values, and the values of
the ``scale`` parameters will be adjusted during the MCMC chain to try to
optimize mixing.

The ``ModelOperator`` as a setting for the ``number_of_auxiliary_categories``.
This controls how many "extra" event categories the Gibbs sampler has when
proposing changes to the number of events, and which comparisons are assigned
to each event.
More "extra" categories can improve mixing, but slows down the MCMC chain;
fewer categories will take less time, at the risk of poorer mixing.
I do not suspect you would ever need more than 4, but you may very well be able
to use 3 or 2, and still have good mixing.
The default is 4.

Generally, the default values for the ``operator_settings`` are sensible
and will work fine for many datasets.
I recommend trying the defaults first (i.e., simply do not specify the
``operator_settings`` section in the configuration file), and make adjustments
if your chains do not mix and/or converge well.

NOTE: you do not have to specify all settings within a group. For
example you can use::

    operator_settings:
        operators:
            ModelOperator:
                number_of_auxiliary_categories: 2

in your config file to change the number of auxiliary categories to 2, but
leave all other operator settings at their default values.


**************************
global_comparison_settings
**************************

The ``global_comparison_settings`` is an optional section that can
be useful for specifying settings to be applied to all of your
comparisons, unless otherwise specified.
All of the settings within ``global_comparison_settings`` can
also be specified for each comparison. For example::

    global_comparison_settings:
        ploidy: 2
        equal_population_sizes: false
    
    comparisons:
    - comparison:
        path: "species1.nex"
    - comparison:
        path: "species2.nex"
    - comparison:
        path: "species3.nex"
        ploidy: 1
        equal_population_sizes: true
    - comparison:
        path: "species4.nex"

specifies that Species 1, 2, and 4 are diploid organisms for which
you want to estimate the root (ancestor) and leaf (descendant) effective
population sizes separately, and Species 3 is haploid and you want to constrain
the root and leaf populations to be the same size.

In other words, ``global_comparison_settings`` allow you to specify default
settings for your comparisons that you can override.

::

        ploidy: 2

This is the ploidy of the organisms (i.e., 1 = haploid, 2 = diploid).

::

        genotypes_are_diploid: true

This tells |eco|_ how you have encoded your characters.
Does each cell of your character matrix represent the state of a diploid
individual?
If so, ``genotypes_are_diploid`` should be ``true``.
If each cell represents the state of a particular gene copy/genome, then
``genotypes_are_diploid`` should be ``false``.
If you have a code(s) to represent a heterozygote, then
``genotypes_are_diploid`` should definitely be ``true``.

::

        markers_are_dominant: false

This specifies whether your markers are dominant.
If the same code is used to designate the character state of a heterozygote and
one of the possible homozygotes, then ``markers_are_dominant`` should be
``true``.

::

        population_name_delimiter: " "
        population_name_is_prefix: false

For each row (individual or gene copy) of your character matrix, you need to
tell |eco|_ which population it was sampled from.
You do this by using prefixes or suffixes for each row label in your nexus
file, and the prefix or suffix needs to be delimited by a character.
So, if your matrix looks like::

    Begin data;
        Dimensions ntax=20 nchar=40000;
        Format datatype=standard symbols="01" missing=? gap=-;
        Matrix
    'population1-jro-001'  0010...
    'population1-jro-002'  0010...
    'population2-jro-003'  0000...
    'population2-jro-004'  0011...
    .
    .
    .

then you should specify::

        population_name_delimiter: "-"
        population_name_is_prefix: true

and |eco|_ will know that the data in the first two rows came from
"population1" and the data in the third and forth rows came from "population2".

.. note:: Underscore gotcha!

    The nexus format standard interprets underscores as spaces, unless the
    labels are quoted. So if you have::

        Begin data;
            Dimensions ntax=20 nchar=40000;
            Format datatype=standard symbols="01" missing=? gap=-;
            Matrix
        population1_jro-001  0010...
        population1_jro-002  0010...
        population2_jro-003  0000...
        .
        .
        .

    you need to specify::

        population_name_delimiter: " " # just a space!
        population_name_is_prefix: true


::

        constant_sites_removed: false

This tells |eco|_ whether or not you have removed all constant
characters/sites.

::

        equal_population_sizes: false

If ``true`` the effective sizes of the root and leaf populations are
constrained be equal (i.e., one population-size parameter for the comparison).
If ``false`` the effective sizes of the root and leaf populations are estimated
separately.


The parameters section
======================

The following section controls the settings for the parameters for each
comparison. Again, these can be specified in the ``global_comparison_settings``
section and/or for each comparison::

        parameters:
            population_size:
                estimate: true
                prior:
                    gamma_distribution:
                        shape: 4.0
                        scale: 0.001
            root_relative_population_size:
                estimate: true
                prior:
                    gamma_distribution:
                        shape: 100.0
                        scale: 0.01
                        offset: 0.0
            freq_1:
                value: 0.5
                estimate: false
            mutation_rate:
                value: 1.0
                estimate: false

This allows you to specify whether or not you want estimate each parameter, and
if so, what prior to use.

mutation_rate
-------------

The ``mutation_rate`` settings are for the mutation rate
(:math:`\mutationRate`) of the comparison.
How you scale this is up to you, but you need to make sure you are consistent
in how you scale time and effective population sizes.
For example, if you set the mutation rate to 1, then time and effective
population sizes will be scaled by the mutation rate.
Specifically, time will be in units of :math:`\eventTime\mutationRate` (i.e.,
expected substitutions per site), and effective population size will be measured
in units of :math:`\effectivePopSize\mutationRate`.
Alternatively, if you specify an actual rate of mutation per site per
generation, then time will be in units of generations,
and population size will in units of the effective number of diploid
individuals or gene copies (:math:`effectivePopSize`) if the ploidy is 2 or 1,
respectively.
Differences in generation times among pairs can also be accounted for
via the ``mutation_rate`` parameters, with the appropriate scaling
of the effective population sizes.
To help ensure the population sizes are scaled correctly, it can help to
remember that :math:`\ploidy 2 \effectivePopSize\mutationRate` should equal the
expected differences per base between two randomly selected genomes from a
population.

population_size
---------------

The ``population_size`` settings are for the effective sizes of the leaf
populations of the pair (i.e., :math:`\effectivePopSize`, but see discussion of
scaling the mutation rate above).

root_relative_population_size
-----------------------------

The population size of the root (ancestral) population is parameterized a bit
differently.
You specify a prior on the effective population size of the root *relative* to
the mean population size of the leaf (descendant) populations.
For example::

            root_relative_population_size:
                value: 1.0
                estimate: false

Constrains the root population to always have an effective population
size that is equal to the mean size of the leaf populations.
Thus, it is not an estimated (free) parameter; it is a deterministic function
of the leaf population sizes.
Similarly, ::

            root_relative_population_size:
                value: 2.0
                estimate: false

constrains the effective population size of the root to be twice
the mean effective population size of the leaves.
Alternatively, ::

            root_relative_population_size:
                estimate: true
                prior:
                    gamma_distribution:
                        shape: 100.0
                        scale: 0.01

allows the effective population size of the root to be estimated, and centers
the prior on its relative size on 1 (i.e., centers the prior expectation for
the actual root effective population size on the mean of the leaf sizes);
the mean of a gamma distribution is the product of the
shape and scale parameters: :math:`100 \times 0.01 = 1`
Similarly ::

            root_relative_population_size:
                estimate: true
                prior:
                    gamma_distribution:
                        shape: 100.0
                        scale: 0.02

allows the effective population size of the root to be estimated, and centers
the prior on its relative size on 2 (i.e., centers the prior expectation for
the actual root effective population size on twice the mean of the leaf sizes.

The hope of this parameterization is to allow you to specify a more informative
prior on the root effective population size.
There is usually a lot of prior uncertainty in the actual value of the root
population size, but we might have good reason to expect that it is similar to
the mean of the leaf sizes.

freq_1
------

The ``freq_1`` parameter is the equilibrium frequency of the "1" allele (or 1
minus the frequency of the "0" allele).
If you are using nucleotide data, I recommend that you fix the frequencies
of the "0" and "1" states to be equal::

            freq_1:
                value: 0.5
                estimate: false

This is because there is no natural way to recode the 4 nucleotide states to
binary.
Thus, if you try to estimate frequencies of the binary states, your results
will be sensitive to the vagaries of how you decide to code your nucleotides as
binary.

However, if the characters you are using are truly binary, then it might make
sense to estimate the frequencies of the two states.
Another option is::

            freq_1:
                value: empirical
                estimate: false

which fixes the frequencies of the two states to their empirical frequencies
(i.e., the frequencies at which they appear in your data).


The operators section
=====================

::

        operators:
            RootPopulationSizeScaler:
                weight: 1.0
                scale: 0.05
            LeafPopulationSizeScaler:
                weight: 1.0
                scale: 0.05
            TimeRootSizeMixer:
                weight: 3.0
                scale: 0.05
