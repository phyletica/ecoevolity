|phyco_logo_long|

.. _phycobackground:

##########
Background
##########

Many processes of diversification can simultaneously affect multiple
species.
For example, below is an animation of three species of lizards co-occurring on
an island that is fragmented twice by rising sea levels.

.. _shared_divs_bifurcating_gif:

.. figure:: /_static/slides-bifurcating.gif
    :align: center
    :width: 99%
    :alt: Shared divergences

    Island fragmentation causing shared divergences.

This creates two bouts of shared divergences across the tree, indicated
by the dashed lines above.
In addition to biogeography, there are many other examples of processes
of diversification that generate patterns of shared divergences.
Instead of lizards on islands, let's imaging three members of a gene family
residing along a region of a chromosome that gets duplicated.
This would create shared divergences across the phylogenetic history of the
gene family.
In epidemiology, when multiple infected individuals spread the pathogen to
others at a social gathering, this will create shared divergences in the
"transmission tree" of the pathogen.

In the animation below, the rising sea levels fragment the island
into three islands.

.. _shared_divs_multifurcating_gif:

.. figure:: /_static/slides-multifurcating.gif
    :align: center
    :width: 99%
    :alt: Shared multifurcating divergences

    Island fragmentation causing multifurcating, shared divergences.

In addition to causing a shared divergence across the tree, this also causes
multifurcations (a lineage diverging into three or more descendants).
Similarly, when an infected individual spreads a pathogen to two or
more others at a social gathering, this will create a multifurcating
divergence in the transmission tree.

Current phylogenetic methods for inferring rooted trees assume all divergences
are independent and bifurcating.
In other words, if we have :math:`\nTips` tips, current methods only consider
trees with :math:`\nTips - 1` independent, bifurcating divergences.
When shared and/or multifurcating divergences were common in the system we
want to study, such tree models are over-parameterized as
:ref:`illustrated in the figure below <true_v_current_tree_model>`.
**What is even worse, by assuming all divergences are independent and
bifurcating, current phylogenetic methods do not allow us to test for patterns
of shared or multifurcating divergences predicted by processes of
diversification that are of interest across the life sciences**.

.. _true_v_current_tree_model:

.. figure:: /_static/gecko-trees-flipped-cropped.svg
    :align: center
    :width: 99%
    :alt: True versus current tree model

    When shared or multifurcating divergences have occurred, current
    phylogenetic models are over-parameterized.

