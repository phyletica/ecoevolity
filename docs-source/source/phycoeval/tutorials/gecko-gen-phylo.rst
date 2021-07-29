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

.. include:: ../../snippets/brief-background.rst

For more detailed background information about the method implemented in |phyco|,
please :ref:`see here <phycobackground>`.

Software Used in this Activity
==============================

For this activity, we will be using |phyco| (part of the |eco| package) and
|pyco|.
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

We will be analyzing "RADseq" data from ...

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
