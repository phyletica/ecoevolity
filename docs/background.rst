.. _background:

##########
Background
##########

How do community-scale processes affect diversification?
Environmental processes that drive speciation across whole communities of
species predict divergence times across taxa that are temporally clustered.
Our goal with |eco|_ is to provide a tool for testing such predictions
in a full-likelihood Bayesian model choice framework
:cite:`Oaks2018ecoevolity`.
:ref:`The cartoon below <div_island_cartoon>` shows an example where two
species of lizards co-diverge when their island is fragmented by rising sea
levels.

.. _div_island_cartoon:

.. figure:: /_static/div-island-cartoon-event-labels.png
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: divergence model cartoon
   
   A cartoon showing three pairs of insular lizard species, two of which
   co-diverged when the island was fragmented.

If we collect genetic data from all six populations, and we want to
evaluate whether the scenario above is a good explanation for the patterns
in our data, there are other possible scenarios (models) that we need
to consider.
For example, perhaps all three pairs of populations diverged at the
same time:

.. _divmodel_111:

.. figure:: /_static/div-model-111-labels.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: divergence model 111

   All three pairs of lizard species co-diverged.

Also, there are three possible ways our three pairs of species could have
diverged at two different times:

.. _divmodel_211:

.. figure:: /_static/div-model-311-labels.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: divergence model 311

.. figure:: /_static/div-model-131-labels.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: divergence model 131

.. figure:: /_static/div-model-113-labels.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: divergence model 113

   Three possible ways our three pairs of lizard species diverged at two
   different times.

Finally, it's also possible that all three pairs diverged independently:

.. _divmodel_123:

.. figure:: /_static/div-model-213-labels.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: divergence model 123

   Three independent divergences.

If we want to test for patterns of shared (non-independent) divergences, this
last scenario of independent divergences seems like a sensible null model.
Overall, we have five possible models that differ in how the three pairs of
species are partitioned into divergence groups (or "events").
In other words, the models differ in the number of divergence-time parameters
and how we assign our species pairs to them.

Now, we need a way to compare how well these five models explain our genetic
data.
Assuming a probabilistic model of how our sequences evolved (more on that in a
bit), we can derive the function for the probability of the data from one of
our pairs of populations given the time of divergence and the effective sizes
of the ancestral and descendant populations (i.e., the likelihood of the
population history).
With the likelihood function in hand, for each of our five possible models, we
can find the combination of parameter values (divergence times and population
sizes) that maximize the likelihood that the model produced our data.
Then, we can compare the ML scores among our models (while in some way
penalizing models with more divergence-time parameters) to select our preferred
model.
.. However, this ML approach does not allow us to say anything about the
.. probability of the models themselves (only the data).
However, as the number of pairs we wish to compare increases, this
ML approach becomes daunting.
The number of possible models we need to consider is
the `Bell number <http://en.wikipedia.org/wiki/Bell_number>`_
:cite:`Bell1934`.
For 5, 10, and 20 pairs, there are 52, 115,975, and 51,724,158,235,372 possible
divergence models, respectively.
A Bayesian model-averaging approach is appealing, because it allows the data to
determine which models are most relevant, and it allows us to make probability
statements about the models given the data.

Methods like msBayes
:cite:`Hickerson2006,Huang2011`
and dpp-msbayes
:cite:`Oaks2014dpp`
attempt to perform this type of model-averaging
using approximate-likelihood Bayesian computation (ABC).
Rather than compute the likelihood, these methods use simulations
to approximate it.
These methods often struggle to detect multiple divergence times across pairs
of populations
:cite:`Oaks2012,Oaks2014reply`,
are very sensitive to prior assumptions
:cite:`Oaks2012,Hickerson2013,Oaks2014reply`,
or have little information to update *a priori* expectations
:cite:`Oaks2014dpp`.
Furthermore, the loss of information inherent to ABC approaches can prevent
them from accurately selecting models
:cite:`Robert2011,Marin2014,Green2015`.

Our goal with |eco|_ is to overcome some of these limitations by accommodating
genomic data (which we'll denote as :math:`\alldata`), and using all the
information in those data.
We do this by computing the full likelihood of the models described above.
|Eco| models biallelic characters from across the genome as having evolved
along their respective gene trees according to continuous time Markov chain of
state change.
It assumes a two-state equivalent of either the Jukes-Cantor (JC69) or general
time-reversible (GTR) of nucleotide substitution
:cite:`JC1969,Tavare1986`.
Furthermore, |eco| assumes each gene tree branched within the populations
according to a coalescent model
:cite:`Kingman1982`.
Given one of the possible divergence models (:math:`\model`) and values of its
parameters (:math:`parameters`; divergence times, effective population sizes,
and state frequencies), |eco| directly computes the probability of the
biallelic characters,
:math:`\pr(\alldata \given \parameters, \model)`.
By doing so, |eco| effectively integrates over all possible gene trees and
character substitution histories along those tries, during the likelihood
calculation.
This frees |eco| from having to estimate the gene trees, and is all thanks to
the work of David Bryant and colleagues :cite:`Bryant2012`.

Using this likelihood function and Markov chain Monte Carlo (MCMC),
|Eco| jointly samples the posterior across all possible models:

.. math::
    :label: jointposterior

    \pr(\parameters, \model[i] \given \alldata) = 
    \frac{
        \pr(\alldata \given \parameters, \model[i])
        \pr(\parameters \given \model[i])
        \pr(\model[i])
    }{
        \pr(\alldata)
        % \sum_{i}\int_{\parameters}
        % \pr(\alldata \given \parameters, \model[i])
        % \pr(\parameters \given \model[i])
        % \diff{\parameters}
        % \pr(\model[i])
    },

Once we have our samples from the joint posterior above, we can marginalize all
the parameters to approximate the posterior probabilities of the divergence
models:

.. math::
    :label: marginalposterior

    \pr(\model[i] \given \alldata) &= 
    \frac{
        \int_{\parameters}
        \pr(\alldata \given \parameters, \model[i])
        \pr(\parameters \given \model[i])
        \diff{\parameters}
        \;\; \pr(\model[i])
    }{
        \pr(\alldata)
        % \sum_{i}
        % \pr(\alldata \given \model[i])
        % \pr(\model[i])
    } \\
    &=
    \frac{
        \pr(\alldata \given \model[i])
        \pr(\model[i])
    }{
        \pr(\alldata)
        % \sum_{i}
        % \pr(\alldata \given \model[i])
        % \pr(\model[i])
    }

We can get Monte Carlo approximations of these probabilities by simply looking
at the proportion of the posterior samples that are from each model.

.. note::

    The term :math:`\pr(\alldata \given \model[i])` Equation
    :eq:`marginalposterior` is absolutely critical.
    This is the *marginal likelihood* of model :math:`i`, which updates
    our prior to give us the posterior probability for model :math:`i`.
    The marginal likelihood is the model's likelihood "averaged" over all
    possible values of the parameters, and the average is weighted by the
    priors on those parameters.
    As a result, the marginal likelihood of each model can be very sensitive to
    the priors we choose for the parameters, even if we have a lot of
    informative data.
    Accordingly, we **highly recommend** you analyze your data multiple times
    using different prior settings (especially the prior on divergence times)
    to assess how sensitive your results are to the priors.


.. _prior_on_divergence_models:

**************************
Prior on divergence models
**************************

To sample from the posterior in Equation :eq:`jointposterior`,
we have to assume a prior on all the possible ways our pairs of species
diverged (divergence models).
|Eco| treats the divergence model (number of divergence events, and the
assignment of the taxa to the events) as a random variable under a Dirichlet
process :cite:`Ferguson1973,Antoniak1974`.
The basic idea of the Dirichlet process is quite simple; we assign species
pairs to divergence events one at a time following a very simple rule.
When assigning the :math:`n^{th}` pair, we assign it to its own
event (i.e., a new divergence event) with probability

.. math::
    :label: dppnewcat

    \frac{\alpha}{\alpha + n -1}

or you assign it to an existing event :math:`x` with probability

.. math::
    :label: dppexistingcat

    \frac{n_x}{\alpha + n -1}

where :math:`n_x` is the number of pairs already assigned to
event :math:`x`.
Let's walk through an example using our three pairs of lizard species.
First, we have to assign our first pair ("A") to a
divergence event with probability 1.0;
let's call this the "blue" divergence event.
Next we assign the second pair ("B") to either a new ("red") divergence
event with probability :math:`\alpha/\alpha + 1` or to the same "blue"
divergence event as the first pair with probability :math:`1/\alpha + 1`.
For this example, let's say it gets assigned to the "blue" event.
Lastly, we assign the third pair ("C") to either a new ("red") divergence
event with probability :math:`\alpha/\alpha + 2` or to the same "blue"
divergence event as the first two pairs with probability :math:`2/\alpha +
2`.

The GIF below illustrates how the these simple rules determine the
prior probability of all five possible models.
Notice toward the end of the animation, as the concentration parameter
increases we place more probability on the divergence models with more
independent divergence events (less shared divergences).

.. _dpp_tree:

.. figure:: /_static/dpp-3-example.gif
    :align: center
    :width: 600 px
    :figwidth: 90 %
    :alt: DPP example

    An example of the Dirichlet process.

Notice that the Dirichlet process prior (DPP) is not motivated by any
biological processes.
Rather, we use it because it is flexible (we can adjust or estimate the
concentration parameter), and mathematically convenient; it allows us to use
Gibbs sampling :cite:`Neal2000` to sample across divergence models.


************************************
Inferring shared demographic changes
************************************

In addition to the divergence models discussed above, with |eco|, you can also
infer shared changes in effective population size.
For example, insular species of lizards might increase their population sizes
when their island coalesces with another island:

.. figure:: /_static/demog-island-cartoon-event-labels.png
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: divergence model cartoon
   
Also, |eco| allows you to infer a mix of divergences and population
expansions/contractions:

.. figure:: /_static/mixed-island-cartoon-event-labels.png
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: divergence model cartoon
   
