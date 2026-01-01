Version 1.1.0
=============

Breaking changes
----------------

-   Implementing smarter handling of missing data.

    Previously, even if only one leaf population was missing data for a
    character, we threw out the whole character (mirroring the behavior of
    beast/SNAPP).  However such a character can be informative about the pop
    sizes and relationships of other populations.

    Now, such a missing character is essentially ignored only for the clade(s)
    from which it is missing.

    This will be needed for introducing migration (phylo networks), where
    parent nodes of reticulating nodes can recceive no allele copies (going
    back in time) from the daughter. So, we need to account for the probability
    of a node having no allele copies.

    Results of analyses using previous versions that relied on the
    ``--relax-missing-sites`` flag will differ from Version 1.1.0 onward.
    In previous versions, this option enabled sites with missing data from one
    or more tip populations to be completely ignored.
    Now, such sites (unless data are missing from ALL tip populations), will be
    used and influence the results.

Changes
-------

-   Updating test suite behavior to avoid running slow tests by default.

-   Updating dev tools. Updating build scripts for dependencies and adding test
    scripts that are more HPC friendly.

-   Adding cladogram ouptut option to ``sumphycoeval``. This makes it easier to
    see shared/multifurcating divergences.

-   Adding tests to confirm that loci (even individual sites) can be disjoint
    in an alignment.
    This should be true, because each site is independent and only ends up
    being a count of total and "red" alleles (i.e., the sequence labels are
    ignored).
    For example, the following two alignments are represented identically as
    biallelic data by ecoevolity (where "pop-1" and "pop-2" are the population
    labels):

        individual-1-locus-1_pop-1  ATT???
        individual-2-locus-1_pop-1  ATT???
        individual-1-locus-1_pop-2  TTT???
        individual-2-locus-1_pop-2  ATA???
        individual-1-locus-2_pop-1  ???GGC
        individual-2-locus-2_pop-1  ???GGC
        individual-1-locus-2_pop-2  ???GGA
        individual-2-locus-2_pop-2  ???TGC

    and

        individual-1_pop-1  ATTGGC
        individual-2_pop-1  ATTGGC
        individual-1_pop-2  TTTGGA
        individual-2_pop-2  ATATGC


Version 1.0.0
=============

Changes
-------

-   Adding ``phycoeval``, ``sumphycoeval``, and ``simphycoeval`` tools and
    associated code and tests. These tools implement fully phylogenetic models
    of shared or multifurcating divergences under a generalized tree
    distribution.

-   Updating documentation.

-   Documenting the Pitman-Yor and Uniform distributions over divergence models
    for ``ecoevolity``.

-   Updating ``nex2yml`` tool to take both ``ecoevolity`` and ``phycoeval``
    config files.


Version 0.3.2
=============

Changes
-------

-   Fixing bug related to NAN likelihoods. Under extreme combinations of
    parameters, NAN likelihoods were being logged, and no error or warning
    occurred. This release fixes that issue. Any analyses prior to this release
    that did not report NAN likelihoods were not affected by this bug.


Version 0.3.1
=============

Changes
-------

-   Fixing bug in simcoevolity that prevented the '--charsets' and
    '--max-one-variable-site-per-locus' options from being used together.
    There was a sanity check from before the '--charsets' option existed that
    needed to be skipped when the '--charsets' option was used.


Version 0.3.0
=============

Changes
-------

-   Adding new simulation option to simcoevolity. This new option allows
    simulating multi-locus datasets (each locus comprised of multiple linked
    sites) that matches the locus lengths and missing data patterns of the
    empirical dataset. The information to simulate such datasets is provided in
    a nexus 'sets' block that delineates the locus boundaries with 'charsets'.


Version 0.2.1
=============

Changes
-------

-   Adding default constructor for MatrixExponentiator to appease Clang.  This
    change does not affect behavior; it simply allows older versions of Clang
    to compile the code.

-   Adding more documentation.


Version 0.2.0
=============

Changes
-------
-   Allowing the MCMC operator 'TimeRootSizeMixer' to be specified as a
    'global' operator (rather than only tree-specific), and making this the
    default behavior. The default behavior before this release was for each
    tree (population pair) to have it's own 'TimeRootSizeMixer' operator.
    Having 'TimeRootSizeMixer' as a global operator should improve mixing when
    pairs share divergence times.


Version 0.1.0
=============

-   Initial release.
