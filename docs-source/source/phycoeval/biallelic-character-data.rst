|phyco_logo_long|

.. _phycodata:

########
The Data
########

To use |phyco|, you will need **biallelic** genetic characters aligned across
the species for which you wish to infer a phylogeny.
|Phyco| assumes these characters are **unlinked**, but, depending on your
data, you might be better off violating this assumption of the model
(more about this :ref:`below <phyco_unlinked>`).
**Biallelic** means that each character has two possible states.
But, most genomic data are composed of nucleotides, which have four states.
This can be accommodated by
:ref:`coding the nucleotides as biallelic <phyco_codingnucs>`.

.. _phyco_codingnucs:

****************************
Coding nucleotide characters
****************************

|Phyco| will read your nucleotide characters, and recode them as biallelic.
If |phyco| finds a polyallelic character with more than two states across your
samples, it will report an error and stop running.
However, you can tell |phyco| to recode these polyallelic characters as either
being the first state (0) or a different state (1).
Your other option is to exclude these sites.

.. note::

    If you are analyzing nucleotide data as biallelic characters, we strongly
    recommend that you do not try and estimate the frequencies of the two
    states.
    Instead, :ref:`fix the frequencies to be equal <phycofreq_1>`.
    This is because there are many ways to recode the 4 nucleotide states to two
    states.
    Thus, if you try to estimate frequencies of the two states, your results can be
    sensitive to how you decided to code your nucleotides as binary.

.. _phyco_unlinked:

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
:ref:`What should we do with such data <phyco_linkage>`?


.. _phyco_linkage:

What to do with linked characters?
==================================

What if you have loci comprising sequences of linked nucleotides?
Based on simulations
:cite:`Oaks2018ecoevolity,Oaks2018paic,Oaks2019codemog,Oaks2021phycoeval`,
we recommend that you analyze all of your
characters (including the constant ones) and violate the assumption of unlinked
characters.
|Eco| and |phyco| perform better when you use all of the sites (including
linked and constant ones) compared to reducing the data to only one variable
character per locus.
So, your best bet is to include all of the sites from loci.


.. _phyco_nexusfile:

********************
Formatting your data
********************

Currently, |phyco| only accepts two input formats for genetic characters:
Nexus and |yaml|_.
We will describe both and then mention some ways to get your data into these
formats.

Nexus format
============

Standard haploid data
---------------------

You can represent your data in a "standard" 0/1 format.
Here's an example of three species from which we've sampled 4 genomes each (2
diploid individuals)::

    #NEXUS
    
    BEGIN TAXA;
        DIMENSIONS NTAX=12;
        TAXLABELS
            lizard-953-a-speciesA
            lizard-953-b-speciesA
            lizard-954-a-speciesA
            lizard-954-b-speciesA
            lizard-152-a-speciesB
            lizard-152-b-speciesB
            lizard-154-a-speciesB
            lizard-154-b-speciesB
            lizard-331-a-speciesC
            lizard-331-b-speciesC
            lizard-338-a-speciesC
            lizard-338-b-speciesC
        ;
    END;

    BEGIN CHARACTERS;
        DIMENSIONS NCHAR=273658;
        FORMAT DATATYPE=STANDARD SYMBOLS="01" MISSING=? GAP=-;
        MATRIX
            lizard-953-a-speciesA  001010...
            lizard-953-b-speciesA  001010...
            lizard-954-a-speciesA  101010...
            lizard-954-b-speciesA  001011...
            lizard-152-a-speciesB  101110...
            lizard-152-b-speciesB  101110...
            lizard-154-a-speciesB  001011...
            lizard-154-b-speciesB  101010...
            lizard-331-a-speciesC  001010...
            lizard-331-b-speciesC  001010...
            lizard-338-a-speciesC  011010...
            lizard-338-b-speciesC  001011...
        ;
    END;

Note, we don't need separate TAXA and CHARACTER blocks like above.
Instead, we can specify a DATA block::

    #NEXUS
    
    BEGIN DATA;
        DIMENSIONS NTAX=12 NCHAR=273658;
        FORMAT DATATYPE=STANDARD SYMBOLS="01" MISSING=? GAP=-;
        MATRIX
            lizard-953-a-speciesA  001010...
            lizard-953-b-speciesA  001010...
            lizard-954-a-speciesA  101010...
            lizard-954-b-speciesA  001011...
            lizard-152-a-speciesB  101110...
            lizard-152-b-speciesB  101110...
            lizard-154-a-speciesB  001011...
            lizard-154-b-speciesB  101010...
            lizard-331-a-speciesC  001010...
            lizard-331-b-speciesC  001010...
            lizard-338-a-speciesC  011010...
            lizard-338-b-speciesC  001011...
        ;
    END;

Both examples above would be equivalent for |phyco|, but the
`Nexus Class Library <http://ncl.sourceforge.net/>`_
used by |phyco| will report a message about an implicit TAXA block if you use the
latter format.
Either way, in your :ref:`phycoeval config file <phycoconfigfile>`,
you need to tell |phyco| that the
:ref:`states, or genotypes, are haploid <phycogenotypesarediploid>`
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
        DIMENSIONS NTAX=6 NCHAR=273658;
        FORMAT DATATYPE=STANDARD SYMBOLS="012" MISSING=? GAP=-;
        MATRIX
            lizard-953-speciesA  002020...
            lizard-954-speciesA  102021...
            lizard-152-speciesB  202220...
            lizard-154-speciesB  102021...
            lizard-331-speciesC  002020...
            lizard-338-speciesC  012021...
        ;
    END;

Now, "0" represents that the individual has two copies with the 0 state, "2"
represents two copies of the 1 state, and "1" represents a heterozygote.
Again, in your :ref:`phycoeval config file <phycoconfigfile>`,
you need to tell |phyco| that the
:ref:`states, or genotypes, are diploid <phycogenotypesarediploid>`
by declaring::

        genotypes_are_diploid: true


Nucleotide data
---------------

If you have nucleotide data, the easiest thing is provide the nucleotide
characters to |phyco| as is, and let it recode them as biallelic.
Here's an example where we are providing nucleotides as haploid (each cell is a
haploid genotype)::

    #NEXUS
    
    BEGIN DATA;
        DIMENSIONS NTAX=12 NCHAR=273658;
        FORMAT DATATYPE=DNA MISSING=? GAP=-;
        MATRIX
            lizard-953-a-speciesA  ACGTAG...
            lizard-953-b-speciesA  ACGTAG...
            lizard-954-a-speciesA  GCGTAG...
            lizard-954-b-speciesA  ACGTAA...
            lizard-152-a-speciesB  GCGCAG...
            lizard-152-b-speciesB  GCGCAG...
            lizard-154-a-speciesB  ACGTAA...
            lizard-154-b-speciesB  GCGTAG...
            lizard-331-a-speciesC  ACGTAG...
            lizard-331-b-speciesC  ACGTAG...
            lizard-338-a-speciesC  ATGTAG...
            lizard-338-b-speciesC  ACGTAA...
        ;
    END;

This is sometimes referred to as "phased" data.
Again, if we are providing a matrix where each cell represents a haploid
genotype, we need to tell |phyco| this is so via the
:ref:`config file <phycoconfigfile>`::

        genotypes_are_diploid: false

We can also represent the same data as "unphased", where each cell represents a
diploid genotype::

    #NEXUS
    
    BEGIN DATA;
        DIMENSIONS NTAX=6 NCHAR=273658;
        FORMAT DATATYPE=DNA MISSING=? GAP=-;
        MATRIX
            lizard-953-a-speciesA  ACGTAG...
            lizard-954-a-speciesA  RCGTAR...
            lizard-152-a-speciesB  GCGCAG...
            lizard-154-a-speciesB  RCGTAR...
            lizard-331-a-speciesC  ACGTAG...
            lizard-338-a-speciesC  AYGTAR...
        ;
    END;

We need to indicate this in the
:ref:`config file <configfile>`::
accordingly::

        genotypes_are_diploid: true


Population labels
-----------------

In our nexus character matrix, we need to indicate which species (or
population) each row corresponds to.
We can do this with either using a prefix or suffix in the row (or taxon)
labels.
For example, in this nexus data file::

    #NEXUS
    
    BEGIN DATA;
        DIMENSIONS NTAX=6 NCHAR=273658;
        FORMAT DATATYPE=DNA MISSING=? GAP=-;
        MATRIX
            lizard-953-a-speciesA  ACGTAG...
            lizard-954-a-speciesA  RCGTAR...
            lizard-152-a-speciesB  GCGCAG...
            lizard-154-a-speciesB  RCGTAR...
            lizard-331-a-speciesC  ACGTAG...
            lizard-338-a-speciesC  AYGTAR...
        ;
    END;

we are using the suffixes to indicate that the first two samples came
from a species we are calling ``speciesA``,
the next to samples came from a species called ``speciesB``,
and
the last two samples came from a speices we are calling
``speciesC``.
In our |phyco| :ref:`config file <phycoconfigfile>`
we have to indicate this with::

        population_name_delimiter: "-"
        population_name_is_prefix: false

This tells |phyco| to look for the last bit of each row label that is
separated by a "-" to figure out the population label.

.. note::

    If you like to use underscores as a population label
    delimiter, just watch out for a
    :ref:`gotcha related to how the nexus format treats underscores <phycounderscoregotcha>`


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
        DIMENSIONS NTAX=12 NCHAR=6;
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
            lizard-331-a-speciesC  ACGTAG
            lizard-331-b-speciesC  ACGTAG
            lizard-338-a-speciesC  ATGTAG
            lizard-338-b-speciesC  ACGTAA
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
        - speciesC
    allele_count_patterns:
        - [[1,4], [3,4], [0,4]]
        - [[0,4], [0,4], [1,4]]
        - [[0,4], [0,4], [0,4]]
        - [[0,4], [2,4], [0,4]]
        - [[1,4], [1,4], [1,4]]
    pattern_weights:
        - 1
        - 1
        - 2
        - 1
        - 1
