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
