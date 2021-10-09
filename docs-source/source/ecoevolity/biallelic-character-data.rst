|eco_logo_long|

.. _data:

########
The Data
########
To use |eco|_, we need to have sampled genetic data from populations we are
interested in comparing.
|Eco| assumes your genetic characters are:

*   :ref:`orthologous <orthologous>`
*   :ref:`biallelic <biallelic>`
*   :ref:`effectively unlinked <unlinked>`

.. _orthologous:

*********
Orthology
*********

**Orthologous** means that you are looking at the "same" character across your
sampled genomes; the character was inherited from a common ancestor without any
duplications of the locus.
**Orthology is only assumed within each comparison**.
For example, if you are comparing the divergence times between two pairs of
species, the characters do not have to be orthologous across both pairs of
species.
You can take advantage of this if you are assembling your data from raw reads.
For example, if you have "RADseq" data, you can assemble your loci for each
comparison separately to increase the amount of data you end up with.


.. _biallelic:

*********
Biallelic
*********

**Biallelic** means that each character has two possible states.
But, most genomic data are composed of nucleotides, which have four states.
This can be accommodated by
:ref:`coding the nucleotides as biallelic <codingnucs>`.

.. _codingnucs:

Coding nucleotide characters
============================

|Eco| will read your nucleotide characters, and recode them as biallelic.
If |eco| finds a polyallelic character with more than two states across your
samples, it will report an error and stop running.
However, you can tell |eco| to recode these polyallelic characters as either
being the first state (0) or a different state (1).
Your other option is to exclude these sites.
For the data we have analyzed so far, there have been no discernible
differences in the results between recoding or excluding characters with more
than 2 states.

.. note::

    If you are analyzing nucleotide data as biallelic characters, we strongly
    recommend that you do not try and estimate the frequencies of the two
    states.
    Instead, :ref:`fix the frequencies to be equal <freq_1>`.
    This is because there are many ways to recode the 4 nucleotide states to two
    states.
    Thus, if you try to estimate frequencies of the two states, your results can be
    sensitive to how you decided to code your nucleotides as binary.

.. _unlinked:

********
Unlinked
********

**Unlinked** means that each character is assumed to have evolved along a gene
tree that is independent from the others (conditional on the population
history).
In other words, it assumes that the characters are far enough apart from one
another in the genome that they segregate independently during meiosis.

Many genomic data sets consist of many loci, each of which is a stretch of
linked nucleotides.
Examples include "RADseq" or sequence-capture data.
:ref:`What should we do with such data <linkage>`?


.. _linkage:

What to do with linked characters?
==================================

What if you have loci comprising sequences of linked nucleotides?
Based on simulations of loci of 100, 500, and 1000 linked characters
:cite:`Oaks2018ecoevolity`, we strongly recommend that you analyze all of your
characters (including the constant ones) and violate the assumption of unlinked
characters.
|Eco| performs much better when you use all of the sites (especially the
constant ones) compared to reducing the data to only one variable character per
locus.
So, your best bet is to put all of your loci for a comparison together into a
nexus-formatted character matrix.

.. _nexusfile:

********************
Formatting your data
********************

Currently, |eco| accepts two input formats for genetic characters:
Nexus and |yaml|_.
The character data for each of your comparisons will go into a separate
file in one of these formats, which are described below.

Nexus format
============

Standard haploid data
---------------------

You can represent your data in a "standard" 0/1 format.
Here's example of a pair of populations from which we've sampled 4 genomes (2
diploid individuals) from two different populations (indicated by the last part
of the taxon label)::

    #NEXUS
    
    BEGIN TAXA;
        DIMENSIONS NTAX=8;
        TAXLABELS
            RMB-5953-a-BabuyanClaro
            RMB-5953-b-BabuyanClaro
            RMB-5954-a-BabuyanClaro
            RMB-5954-b-BabuyanClaro
            RMB-6052-a-Calayan
            RMB-6052-b-Calayan
            RMB-6054-a-Calayan
            RMB-6054-b-Calayan
        ;
    END;

    BEGIN CHARACTERS;
        DIMENSIONS NCHAR=273658;
        FORMAT DATATYPE=STANDARD SYMBOLS="01" MISSING=? GAP=-;
        MATRIX
            RMB-5953-a-BabuyanClaro     001010...
            RMB-5953-b-BabuyanClaro     001010...
            RMB-5954-a-BabuyanClaro     101010...
            RMB-5954-b-BabuyanClaro     001011...
            RMB-6052-a-Calayan          101110...
            RMB-6052-b-Calayan          101110...
            RMB-6054-a-Calayan          001011...
            RMB-6054-b-Calayan          101010...
        ;
    END;

Note, we don't need separate TAXA and CHARACTER blocks like above.
Instead, we can specify a DATA block::

    #NEXUS
    
    BEGIN DATA;
        DIMENSIONS NTAX=8 NCHAR=273658;
        FORMAT DATATYPE=STANDARD SYMBOLS="01" MISSING=? GAP=-;
        MATRIX
            RMB-5953-a-BabuyanClaro     001010...
            RMB-5953-b-BabuyanClaro     001010...
            RMB-5954-a-BabuyanClaro     101010...
            RMB-5954-b-BabuyanClaro     001011...
            RMB-6052-a-Calayan          101110...
            RMB-6052-b-Calayan          101110...
            RMB-6054-a-Calayan          001011...
            RMB-6054-b-Calayan          101010...
        ;
    END;

Both examples above would be equivalent to |eco|, but the
`Nexus Class Library <http://ncl.sourceforge.net/>`_
used by |eco| will report a message about an implicit TAXA block if you use the
latter format.
Either way, in your :ref:`ecoevolity config file <configfile>`,
you need to tell |eco| that the
:ref:`states, or genotypes, are haploid <genotypesarediploid>`
by declaring::

        genotypes_are_diploid: false


Standard diploid data
---------------------

Above, each cell in our matrix represented which state was present
for the character in a particular haploid genome.
We can also represent the same data where each cell represents
the genotype of a diploid individual::

    #NEXUS
    
    BEGIN DATA;
        DIMENSIONS NTAX=4 NCHAR=273658;
        FORMAT DATATYPE=STANDARD SYMBOLS="012" MISSING=? GAP=-;
        MATRIX
            RMB-5953-BabuyanClaro     002020...
            RMB-5954-BabuyanClaro     102021...
            RMB-6052-Calayan          202220...
            RMB-6054-Calayan          101021...
        ;
    END;

Now, "0" represents that the individual has two copies with the 0 state, "2"
represents two copies of the 1 state, and "1" represents a heterozygote.
Again, in your :ref:`ecoevolity config file <configfile>`,
you need to tell |eco| that the
:ref:`states, or genotypes, are diploid <genotypesarediploid>`
by declaring::

        genotypes_are_diploid: true


Nucleotide data
---------------

If you have nucleotide data, the easiest thing is provide the nucleotide
characters to |eco| as is, and let it recode them as biallelic.
Here's an example where we are providing nucleotides as haploid (each cell is a
haploid genotype)::

    #NEXUS
    
    BEGIN DATA;
        DIMENSIONS NTAX=8 NCHAR=273658;
        FORMAT DATATYPE=DNA MISSING=? GAP=-;
        MATRIX
            RMB-5953-a-BabuyanClaro     ACGTAG...
            RMB-5953-b-BabuyanClaro     ACGTAG...
            RMB-5954-a-BabuyanClaro     GCGTAG...
            RMB-5954-b-BabuyanClaro     ACGTAA...
            RMB-6052-a-Calayan          GCGCAG...
            RMB-6052-b-Calayan          GCGCAG...
            RMB-6054-a-Calayan          ACGTAA...
            RMB-6054-b-Calayan          GCGTAG...
        ;
    END;

This is sometimes referred to as "phased" data.
Again, if we are providing a matrix where each cell represents a haploid genotype,
we need to tell |eco| this is so via the
:ref:`config file <configfile>`::

        genotypes_are_diploid: false

We can also represent the same data as "unphased", where each cell represents a
diploid genotype::

    #NEXUS
    
    BEGIN DATA;
        DIMENSIONS NTAX=4 NCHAR=273658;
        FORMAT DATATYPE=DNA MISSING=? GAP=-;
        MATRIX
            RMB-5953-BabuyanClaro     ACGTAG...
            RMB-5954-BabuyanClaro     RCGTAR...
            RMB-6052-Calayan          GCGCAG...
            RMB-6054-Calayan          RCGTAR...
        ;
    END;

We need to indicate this in the
:ref:`config file <configfile>`::
accordingly::

        genotypes_are_diploid: true


Population labels
-----------------

In our nexus character matrix, we need to indicate which population each
row corresponds to.
We can do this with either using a prefix or suffix in the row (or taxon)
labels.
For example, in this example::

    #NEXUS
    
    BEGIN DATA;
        DIMENSIONS NTAX=4 NCHAR=273658;
        FORMAT DATATYPE=DNA MISSING=? GAP=-;
        MATRIX
            RMB-5953-BabuyanClaro     ACGTAG...
            RMB-5954-BabuyanClaro     RCGTAR...
            RMB-6052-Calayan          GCGCAG...
            RMB-6054-Calayan          RCGTAR...
        ;
    END;

we are using the suffixes to indicate that the first two samples came
from a population we are calling ``BabuyanClaro``, and
the last two samples came from a population we are calling
``Calayan``.
In our :ref:`|eco| config file <configfile>`
we have to indicate this with::

        population_name_delimiter: "-"
        population_name_is_prefix: false

This tells |eco| to look for the last bit of each row label that is
separated by a "-" to figure out the population label.

Every nexus file must one or two population labels.
If |eco| finds two population labels, it will model
the comparison as two diverged populations and try to
estimate the time that they diverged:

.. image:: /_static/div-model-singleton.svg
   :align: center
   :width: 50%
   :alt: divergence comparison cartoon

If |eco| finds one population label, it will
model the comparison as a population that underwent
a population-size change and try to estimate the time
that the change occurred:

.. image:: /_static/demog-model-singleton.svg
   :align: center
   :width: 50%
   :alt: demog comparison cartoon

.. note::

    If you like to use underscores as a population label
    delimiter, just watch out for a
    :ref:`gotcha related to how the nexus format treats underscores <underscoregotcha>`

YAML format
===========

The |yaml|-formatted data format is a lot more efficient (much smaller file
sizes).
Instead of a full alignment,
it only contains a list of allele count patterns, followed by a list of the
weight of each pattern (i.e., how many times the allele pattern occurs in the
alignment).
For a very small example let's convert the following nexus alignment
and convert it to the |yaml|_ format that |eco| accepts as input::

    #NEXUS
    
    BEGIN DATA;
        DIMENSIONS NTAX=8 NCHAR=6;
        FORMAT DATATYPE=DNA MISSING=? GAP=-;
        MATRIX
            lizard-953-a-speciesA  ACGTAG
            lizard-953-b-speciesA  ACGTAG
            lizard-954-a-speciesA  GCGTAG
            lizard-954-b-speciesA  ACGTAA
            lizard-152-a-speciesB  GCGCAG
            lizard-152-b-speciesB  GCGCAG
            lizard-154-a-speciesB  ACGTAA
            lizard-154-b-speciesB  GCGTAG
        ;
    END;

To convert these data to |yaml|_ format, we will assume the first nucleotide
(from the top) in each column is state "0", and the second nucleotide (if any)
is state "1".
Doing so gives us the following allele-count patterns in |yaml|_ format::

    ---
    markers_are_dominant: false
    population_labels:
        - speciesA
        - speciesB
    allele_count_patterns:
        - [[1,4], [3,4]]
        - [[0,4], [0,4]]
        - [[0,4], [2,4]]
        - [[1,4], [1,4]]
    pattern_weights:
        - 1
        - 3
        - 1
        - 1
