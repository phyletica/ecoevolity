|phyco_logo_long|

.. _phycoconfigfile:

######################
The Configuration File
######################

.. contents:: Configuration File Components
    :local:
    :depth: 3

You specify most of the settings for a |phyco|_ analysis in a |yaml|_-formatted
configuration file.
|yaml|_ is an "anti-markup" language and data standard that is much more human
friendly than alternatives like XML.

.. note::

    The website |yamllint| is a nice tool for debugging |yaml|_ syntax.
    You can copy and paste your config file their to check of you're
    using valid |yaml| syntax.

The simplest possible configuration file for |phyco|_ would look something like::

    ---
    data:
        constant_sites_removed: false
        alignment:
            path: alignment/Gekko.nex


This specifies the path to the nexus file containing the aligned character data
for the species we want to infer a phylogeny for, and specifies whether we have
removed constant sites from the alignment.
The defaults for all other settings would be used in this case, which we
**strongly discourage**.

Now, let's look at an example that is more thoroughly specified::

    ---
    data:
        ploidy: 2
        constant_sites_removed: false
        alignment:
            genotypes_are_diploid: true
            markers_are_dominant: false
            population_name_is_prefix: false
            population_name_delimiter: ' '
            path: alignment/Gekko.nex
    # This is a comment!
    tree_model:
        tree_space: generalized
        starting_tree: comb
        tree_prior:
            uniform_root_and_betas:
                parameters:
                    root_height:
                        estimate: true
                        prior:
                            gamma_distribution:
                                shape: 4.0
                                mean: 0.01
                    alpha_of_node_height_beta_prior:
                        value: 1.0
                        estimate: false
    branch_parameters:
        population_size:
            equal_population_sizes: true
            value: 0.0002
            estimate: true
            prior:
                gamma_distribution:
                    shape: 4.0
                    mean: 0.0002
    mutation_parameters:
        freq_1:
            value: 0.5
            estimate: false
        mutation_rate:
            value: 1.0
            estimate: false
    mcmc_settings:
        chain_length: 15000
        sample_frequency: 10

All the settings are hierarchically nested by the indent spacing.
For example,
``data``, ``tree_model``, ``branch_parameters``,
``mutation_parameters``, and ``mcmc_settings``
are all at the highest level of the hierarchy, and have various settings nested
within them.
Across a given level of the hierarchy, order does not matter. E.g., the
top-level groups of settings listed above can be arranged in any order.
Anything proceeded by a '#' is a comment that is ignored by |phyco|_.

.. _param_syntax:

************************
General parameter syntax
************************

Before we get into the specific syntax for each setting in the config file,
let's first
look at the general syntax that applies
to all parameters in the ``tree_model``, ``branch_parameters``, and
``mutation_parameters`` sections of the config.

The general syntax for a parameter is::

    parameter_name:
        estimate: true # or false
        value: 1.0
        prior:
            a_valid_distribution:
                distribution_parameter: 1.0

So, for any parameter, you can specify:

#.  Whether or not the parameter should be fixed (``estimate: false``) or
    estimated (``estimate: true``).
#.  A value for the parameter. This is only the starting value if ``estimate``
    is ``true``, or is the fixed value if ``estimate`` is ``false``.
    **If a value is not specified, the starting value is drawn from the
    prior**.
#.  The prior probability distribution. This is ignored if ``estimate`` is
    ``false``.


**********
tree_model
**********

Next, let's walk through how to specify the tree model in the configration
file.

::

    tree_model:
        tree_space: generalized
        starting_tree: comb
        tree_prior:
            uniform_root_and_betas:
                parameters:
                    root_height:
                        estimate: true
                        prior:
                            gamma_distribution:
                                shape: 4.0
                                mean: 0.01
                    alpha_of_node_height_beta_prior:
                        value: 1.0
                        estimate: false

tree_space
==========

There are three options for the ``tree_space`` setting:

#.  ``tree_space: generalized``

    This tells |phyco| to use a generalized tree distribution that allows for
    shared and/or multifurcating divergences.

#.  ``tree_space: bifurcating``

    This specifies a "standard" tree distribution that assumes all
    divergences are independent and bifurcating.

#.  ``tree_space: fixed``

    This tells |phyco| that you want to fix the tree topology.  To do this, you
    must also provide a starting tree with the topology you wish to hold
    constant during the analyses.


starting_tree
=============

|phyco| will first check to see if ``starting_tree`` is set to one of two
options:

#.  ``starting_tree: comb``

    This specifies that the MCMC chain should start with the tree topology with
    only one multifurcating internal node.

#.  ``starting_tree: random``

    This specifies that the MCMC chain should start with the random topology
    with only independent, bifurcating divergences; i.e., a random topology
    with :math:`\nTips - 1` internal nodes, where :math:`\nTips` is the number
    of tips (populations/species).

If the setting for ``starting_tree`` is not ``comb`` or ``random``,
|phyco| will next check to see if the string provided is a valid path
to a file.
If it is, the starting tree will be read from the file, which can be in
phylip/newick format or nexus format.

If the setting for ``starting_tree`` is not ``comb``, ``random``, nor a path,
|phyco| will try to interpret the setting as a newick-formatted tree string.

So, in summary, you can specify a specific starting tree using a newick string
or a path to a newick or nexus formatted tree, or you can specify the MCMC
chain start with the comb tree or a random bifurcating tree.

tree_prior
==========

Currently, the only ``tree_prior`` implemented in |phyco| is
``uniform_root_and_betas``.
This prior assumes that all topologies are equally probable *a priori*,
the age of the root follows a parametric distribution (e.g., a gamma or
exponential distribution), and all other divergence times follow scaled beta
distributions.
When, ``tree_space`` is set to ``generalized``, all possible non-reticulating
topologies (including topologies with shared and/or multifurcating divergences)
are *a priori* equally probable.
When, ``tree_space`` is set to ``bifurcating``, all possible topologies with
independent, bifurcating divergences (i.e., :math:`\nTips - 1` internal nodes)
are equally probable.

The prior settings for the root age (height) and the beta distributions
for all other divergence times are specified under the ``parameters``
section fo the ``uniform_root_and_betas`` tree prior::

            uniform_root_and_betas:
                parameters:
                    root_height:
                        estimate: true
                        prior:
                            gamma_distribution:
                                shape: 4.0
                                mean: 0.01
                    alpha_of_node_height_beta_prior:
                        value: 1.0
                        estimate: false

root_height
-----------

The ``root_height`` parameter follows the 
:ref:`general parameter syntax discussed above <param_syntax>`.
Valid continuous probability distributions for the ``prior`` include:

#.  ``exponential_distribution``

    For the ``exponential_distribution``, you can specify the
    ``mean`` or ``rate`` (the rate = 1/mean), and an ``offset``.
    By default, the ``offset`` is zero, and the lower bound of the
    exponential distribution is zero. If you specify a positive
    offset, this shifts the entire distribution so that the lower
    bound is equal to ``offset``.
    **NOTE**: If you specify a ``mean`` and an ``offset``, the mean is
    interpreted as the mean **before** shifting the distribution.

#.  ``gamma_distribution``

    For the ``gamma_distribution``, you can specify the
    ``shape`` and ``mean`` **OR** the ``shape`` and ``scale``,
    where the
    :math:`\textrm{mean} = \textrm{shape} \times \textrm{scale}`.
    As with the exponential distribution, you can also specify
    an ``offset``, which will shift the gamma distribution
    to have a lower bound equal to ``offset``.
    **NOTE**: When you specify an ``offset``, the
    ``scale`` OR ``mean`` is interpreted as the
    scale or mean **before** shifting the distribution.

#.  ``uniform_distribution``

    The ``uniform_distribution`` has two parameters that can be specified: the
    ``min`` and ``max``.


alpha_of_node_height_beta_prior
-------------------------------

The ``uniform_root_and_betas`` tree prior assumes (*a priori*)
that all non-root divergence times follow a scaled beta distribution
between the present and the age of the youngest parent node of any
node mapped to the divergence time
(:ref:`see the figure below for an example <scaled_beta_divs_fig>`.

.. _scaled_beta_divs_fig:

.. figure:: /_static/rj-move-tree-cropped.svg
    :align: center
    :width: 99%
    :alt: Tree prior

    An example of the scaled beta distributions on non-root divergence times.

These are "scaled" betas, because, normally, a beta distribution has a lower
and upper bound of 0 and 1, respectively.
The beta distributions are scaled to be proper probability distributions (i.e.,
integrate to one) between zero (the "present") and the age of the youngest
parent node.
By default, |phyco| sets both the alpha and beta parameters of these beta
distributions to one, which makes them equivalent to uniform distributions
between zero and the age of the youngest parent node.
However, |phyco| allows you to control the settings for the alpha parameter of
these beta distributions via the ``alpha_of_node_height_beta_prior`` option::

                    alpha_of_node_height_beta_prior:
                        value: 1.0
                        estimate: false

Above, we are simply fixing the alpha parameter of the beta distributions on
the non-root divergence times to 1 (this is the default for |phyco|).
In general, the default is probably fine for most applications, and estimating
``alpha_of_node_height_beta_prior`` has not been well tested with simulated
data.


*****************
branch_parameters
*****************

::

    branch_parameters:
        population_size:
            equal_population_sizes: true
            value: 0.0002
            estimate: true
            prior:
                gamma_distribution:
                    shape: 4.0
                    mean: 0.0002

population_size
===============

The ``population_size`` parameter follows the 
:ref:`general parameter syntax discussed above <param_syntax>`,
with the exception of one additional setting: ``equal_population_sizes``.
If ``equal_population_sizes`` is set to ``true``, then all the branches in the
tree will share the same ``population_size`` parameter.
If ``equal_population_sizes`` is set to ``false``, then each branch in the tree
gets its own ``population_size`` parameter.

If you set the mutation rate to 1, then the effective population sizes
will be scaled by the mutation rate (:math:`\epopsize\murate`).
Alternatively, if you specify an actual rate of mutation per site per
generation, then the population size will be in units of the effective number
of diploid individuals or haploid genomes (:math:`\epopsize`), if the ploidy is 2
or 1, respectively.

.. _popsizenote:

.. note::

    **Important**: In |eco|, the ``population_size`` is related to, but **not**
    equal to :math:`\theta` (:math:`4\epopsize\murate`; the genetic
    diversity, or more precisely, the expected number of differences per base
    between two randomly selected haploid genomes).
    The relationship between ``population_size`` (represented by
    :math:`\epopsize`) and :math:`\theta` is:

    .. math::
        :label: phycothetarelationship
    
        \textrm{ploidy} \times 2\epopsize\murate = \theta \\
        \epopsize = \frac{\theta}{\textrm{ploidy} \times 2\murate}.

    Thus, if you have prior expectation that :math:`\theta = 0.002` and you've
    set the ``mutation_rate`` to 1.0 (i.e., :math:`\murate = 1`) and
    ``ploidy`` to 2, then your prior expectation for ``population_size`` is,

    .. math::
    
        \epopsize = \frac{0.002}{2 \times 2(1)} = \frac{0.002}{4} = 0.0005

The relationship above is also very important to keep in mind if you
specify a mutation rate in years.
For example, let's assume that for your taxa, you believe 0.004 is a
reasonable value for the average number of differences per base between two
randomly selected gene copies, and you've told |phyco| that the ``ploidy =
2``.
If you've set the mutation rate to ``1e-8``, then you can
use Equation :eq:`phycothetarelationship` to adjust your prior
expectation for ``population_size`` accordingly:

.. math::

    \epopsize &= \frac{\theta}{\textrm{ploidy} \times 2\murate} \\
              &= \frac{0.004}{2 \times 2(0.00000001)} = \frac{0.002}{2} = \textrm{100,000}


*******************
mutation_parameters
*******************

::

    mutation_parameters:
        freq_1:
            value: 0.5
            estimate: false
        mutation_rate:
            value: 1.0
            estimate: false

mutation_rate
=============

The ``mutation_rate`` settings are for rate of mutation across the tree.
.. How you scale this is up to you, but you need to make sure you are consistent
.. in how you scale time and effective population sizes.
If you set the mutation rate to 1, then time and effective
population sizes will be scaled by the mutation rate.
Specifically, time will be in units of :math:`\divtime\murate` (i.e.,
expected substitutions per site), and effective population size will be measured
in :math:`\epopsize\murate`.
Alternatively, if you specify an actual rate of mutation per site per
generation, then time will be in units of generations,
and population size will be in units of the effective number of diploid
individuals or gene copies (:math:`\epopsize`) if the ploidy is 2 or 1,
respectively.
To help ensure the population sizes are scaled correctly, it can help to
remember that :math:`\textrm{ploidy} \times 2\epopsize\murate` should
equal the expected differences per base between two randomly selected genomes
from a population.

.. _phycofreq_1:

freq_1
======

The ``freq_1`` parameter is the equilibrium frequency of the "1" allele (or 1
minus the frequency of the "0" allele).
**If you are using nucleotide data, we recommend that you fix the frequencies
of the "0" and "1" states to be equal**::

            freq_1:
                value: 0.5
                estimate: false

This is because there is no natural way to recode the 4 nucleotide states to
two states.
Thus, if you try to estimate frequencies of the two states, your results will
be sensitive to the vagaries related to how you decided to code your
nucleotides as binary.

However, if the characters you are using are truly biallelic, then it might
make sense to estimate the frequencies of the two states::

            freq_1:
                estimate: true
                prior:
                    beta_distribution:
                        alpha: 2.5 
                        beta: 1.2

Another option is::

            freq_1:
                value: empirical
                estimate: false

which fixes the frequencies of the two states to their empirical frequencies
(i.e., the frequencies at which they appear in your data).
Again, you shouldn't do this for nucleotide data.

.. note::
    
    The ``empirical`` option for the value only works for the ``freq_1``
    parameter.  You should get an error if you try to use it for any other
    parameters.


****
data
****

This section of the configuration file tells |phyco| about the data we are
asking it to analyze, and where it can find those data.

::

        ploidy: 2

This is the ploidy of the organisms (e.g., 1 = haploid, 2 = diploid), and
determines the meaning of the ``population_size`` parameters covered
under ``branch_parameters`` below.
If ``ploidy`` is set to one, then ``population_size`` will be the
(potentially mutation-scaled) haploid effective population size
(i.e., the effective number of haploid genomes).
If your taxa are diploid and you set ``ploidy`` to two, then
``population_size`` will be the (potentially mutation-scaled) effective number
of diploid individuals.

::

        constant_sites_removed: false

This tells |phyco| whether or not you have removed all constant
characters/sites.
If ``true``, |phyco| will correct the likelihood for having only
sampled variable characters.

The ``alignment`` section provides |phyco| with information about
the data you wish to analyse::

        alignment:
            genotypes_are_diploid: true
            markers_are_dominant: false
            population_name_is_prefix: false
            population_name_delimiter: ' '
            path: alignment/Gekko.nex

Below, we cover each of the possible settings nested under ``alignment``.


.. _phycogenotypesarediploid:

::

        genotypes_are_diploid: true

This tells |eco|_ how you have encoded your characters.
Does each cell of your character matrix represent the state of both alleles of
a diploid individual?
If so, ``genotypes_are_diploid`` should be ``true``.
If each cell represents the state of a particular gene copy, then
``genotypes_are_diploid`` should be ``false``.
If you have a code(s) to represent a heterozygote, then
``genotypes_are_diploid`` should definitely be ``true``.

::

        markers_are_dominant: false

This specifies whether your markers are dominant.
If the same code is used to designate the character state of a heterozygote and
one of the two possible homozygotes, then ``markers_are_dominant`` should be
``true``.
If you can tell the difference among a heterozygote and both homozygotes, then
``markers_are_dominant`` should be ``false``.

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
    'population1-lizard-001'  0010...
    'population1-lizard-002'  0010...
    'population2-lizard-003'  0000...
    'population2-lizard-004'  0011...
    .
    .
    .

then you should specify::

        population_name_delimiter: "-"
        population_name_is_prefix: true

and |eco|_ will know that the data in the first two rows came from
"population1" and the data in the third and forth rows came from "population2".

.. _phycounderscoregotcha:
.. note:: **Underscore gotcha!**

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



*************
mcmc_settings
*************

The ``mcmc_settings`` simply specify how long to run the MCMC
chain, and how often to record a sample from it. For example::

    mcmc_settings:
        chain_length: 15000
        sample_frequency: 10

tells |phyco| to run the chain for 15,000 generations, recording a sample every
10th generation (15,000/10 = 1500 total samples).
This is a good starting point. For most datasets we have analyzed so far,
this has been sufficient. If the chain is having mixing problems, then you
can try increasing these numbers.
We recommend running several independent chains (analyses) to:

#.  Confirm the chains are converging.
#.  Increase the number of samples from the posterior distribution (assuming
    the chains converged).
