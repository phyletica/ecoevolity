.. _selecting_priors:

****************
Selecting Priors
****************

blah

.. _gamma_intro:

An introduction to the gamma probability distribution
=====================================================

.. code-block:: r

    > x.max = qgamma(0.999, shape=1.0, scale=10.0)
    > x = seq(from=0, to=x.max, by=x.max/1000)
    > dens = dgamma(x, shape=1.0, scale=10.0)
    > plot(x, dens, type='l')

Important priors for the |eco|_ model
=====================================

.. _concentration_parameter:

Concentration parameter of the Dirichlet process
------------------------------------------------

We have to choose a gamma-distributed prior for the concentration parameter
(:math:`\alpha`) of the Dirichlet process that controls the assignment of taxa
to divergence events.

An important point about the concentration parameter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is important to note that values of the DPP concentration parameter are
always specific to the number taxa.
For example, above we saw that for 10 taxa, a concentration parameter
of 3.3 corresponded to a prior mean of 5 divergence events.
However, when the number of taxa is 20, a concentration parameter of
3.3 does **NOT** correspond to a prior mean of 5
    
It actually corresponds to a prior mean of about 6.9.
So, you cannot simply choose your "favorite" prior for the concentration
parameter and apply it blindly for all datasets.
When you are analyzing a dataset with a different number of taxa, you need to
reassess your prior on the concentration parameter.


.. _population_size:

Population size
---------------

Another important parameter for which we need to choose a prior is the
effective population size of the ancestral and descendant populations
in the model.

Let's see how much of our prior probability will fall on values greater
than 0.004 if we use an exponential with a mean of 0.002.
Because the shape parameter is 1, this means the scale parameter is simply
0.002 (remember, the mean is simply the product of the shape and scale):


.. code-block:: r

    > pgamma(0.004, shape=1.0, scale=0.002, lower.tail=F)
    [1] 0.1353353

Perhaps this seems like too much prior probability greater than 0.004; let's
try a scale parameter of 0.001:

.. code-block:: r

    > pgamma(0.004, shape=1.0, scale=0.001, lower.tail=F)
    [1] 0.01831564

If this seems to fit our prior expectations, we can take a look at this prior:

.. code-block:: r

    x.max = qgamma(0.999, shape=1.0, scale=0.001)
    x = seq(from=0, to=x.max, by=x.max/1000)
    dens = dgamma(x, shape=1.0, scale=0.001)
    plot(x, dens, type='l')

which will give us something like:


If you feel you have more prior knowledge than is represented by this
exponential (for example, perhaps you expect the effective population size to
be greater than 0.0005), then you can increase the shape parameter accordingly,
until you end up with a distribution that fits your prior uncertainty.

.. _divergence_time:

Divergence time
---------------

We also need to choose a gamma-distributed prior for the divergence
times of the pairs of populations.

