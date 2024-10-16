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
