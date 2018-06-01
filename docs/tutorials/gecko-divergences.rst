
.. _gecko_div_tutorial:

*****************
Gecko Divergences
*****************

Objective
=========

Our goal is to learn how to set up, run, and interpret an |eco| analysis in
order to test models of shared divergence times among pairs of populations.

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

While not required, you may want to install the program
|Tracer|_ (http://tree.bio.ed.ac.uk/software/tracer/).

The Data
========

We will be analyzing "RADseq" data from three pairs of gecko populations from
the Philippine Islands.
Each pair of population inhabit two oceanic islands in the Philippines:

#.  Babuyan Claro and Calayan
#.  Maestre De Campo and Masbate
#.  Dalupiri and Camiguin Norte

Based on previous data :cite:`Siler2012,Siler2014kikuchii`,
each of these pairs diverged very recently.
However, because these islands were never connected, the distribution
across both islands must have been due to over-water dispersal.
Such rare dispersal events were likely idiosyncratic to each pair of island
populations, and thus we do *not* expect to see shared divergences.
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
    steps to be able to access the example files you just downloaded from within
    the container.
    After you've used ``cd`` to navigate into the example data directory,
    execute the following command to enter the Docker container::
    
        docker run -v "$(pwd)":/portal -it phyletica/ecoevolity-docker bash
    
    .. note::
    
        Depending on your system and how Docker is configured, you may need to use
        ``sudo`` to run Docker commands. If you received a "permission denied"
        message when you ran the command above, try::
        
            sudo docker run -v "$(pwd)":/portal -it phyletica/ecoevolity-docker bash
    
    Then, once inside the docker container, type::
    
        cd portal

Inside the example data directory, you should find a configuration file
``ecoevolity-config.yml`` and a ``nexus-files`` directory containing the
nexus-formatted RADseq data for each of our three pairs of gecko populations.
For detailed information about the nexus-formatted data files needed for
|eco|,
:ref:`click here <data>`.

|Eco| assumes all of your characters have at most two states (biallelic), so
you might be wondering why the nexus files contain nucleotide characters
(i.e., four states).
Very astute of you.
If we provide |eco| with nucleotide data, it will automatically recode the data
as biallelic by considering the first character it finds in each column as
state ``0``, and if it finds a second state it considers it state ``1``.
If it finds a third state in a column, it will report an error message and
quit.
For characters (columns) with more than two states, you have two options:

#.  remove these columns from your alignment, or
#.  recode them as biallelic (more on this in a bit).

|Eco| also assumes all your characters are unlinked, but the gecko data are 200
loci each comprising about 90 linked sites.
Based on simulations of loci of 100, 500, and 1000 linked characters
:cite:`Oaks2018ecoevolity`, we strongly recommend that you analyze all of your
characters (including the constant ones) and violate the assumption of unlinked
characters.
In short, |Eco| performs better when you use all of the sites (especially the
constant ones) compared to reducing the data to only one variable character per
locus.

Within each pair, we assume each column represents an orthologous character
across your samples from both populations.
However, your characters do not need to be orthologous among your different
pairs.

Running |eco|
=============

Setting up the configuration file
---------------------------------

Once we have our nexus files ready, the next step in an |eco| analysis
is setting up the configuration file.
|Eco| requires a |yaml|_-formatted configuration file.
|yaml|_ is a human-friendly data standard that allows you to provide |eco| the
information it needs in a format that is easy for you to read and edit.

.. note::

    The website |yamllint| is a nice tool for debugging |yaml|_ syntax.
    You can copy and paste your config file their to check of you're
    using valid |yaml| syntax.


For all the details about |eco| config files,
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

This tells |eco| that for the prior on event models, we want to 
use a Dirichlet process with its concentration parameter
fixed to 1.414216

Choosing a prior on divergence times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running the analysis
--------------------

The output
==========


Plotting the results
====================

