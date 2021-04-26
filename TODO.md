# Allowing leaves with missing data

Merge changes to BiallelicPatternProbabilityMatrix and likelihood functions
from retic branch into dev branch.
These changes should allow leaves to have no data for characters.
These changes need to be tested after they are merged in; both to make sure
nothing is broken and to add new tests to make sure the missing data are
handled correctly.
To merge them, run this from the dev branch:

    git checkout retic -- src/ecoevolity/matrix.* src/ecoevolity/likelihood.*

Also, the BiallelicData class will need to be updated to allow characters for
which one or more leaves have no data.

# Cached likelihoods
-   Make current BasePopulationTree templated, and then make a class that
    simply specifies PopulationNode
    -   This will allow us to create a 'CachePopNode' that has a member that is
        a vector of vectors of top/bottom probs
        -   Inner vectors would allow for rate categories.
