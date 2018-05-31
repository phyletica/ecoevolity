.. _background:

##########
Background
##########

How do community-scale processes affect diversification?
Environmental processes that drive speciation across whole communities of
species predict divergence times across taxa that are temporally clustered.
Our goal with |eco|_ is to provide a tool for testing such predictions
in a full-likelihood Bayesian model choice framework.
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

   All three lizard lineages diverged when the island fragmented.

As in the :ref:`first scenario <div_island_cartoon>` above perhaps there
were two divergences.
There are three possible ways our three species pairs could have
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

   Three possible ways our three pairs of gecko species diverged at two
   different times.

Lastly, it's also possible that all three pairs diverged independently:

.. _divmodel_123:

.. figure:: /_static/div-model-213-labels.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: divergence model 123

   Three independent divergences.

If we want to test for shared (non-independent) divergences, this last scenario
of independent divergences seems like a sensible null model.
To compare the likelihood of these models...

If researchers are interested in comparing the divergence times among a number
of pairs of populations,
we can approach this as a problem of model choice:
How many divergence events, and what assignment of taxa to those events, best
explain the genetic variation within and between the diverged populations of
each pair (Figure~\ref{fig:divCartoon})?

Biogeographers are often interested in understanding how large-scale processes
affect diversification and community assembly.
One way to approach this challenge is to infer the history of diversification
across co-distributed species and test for patterns predicted by historical
processes of interest (e.g., changes in climate fragmenting communities).
For example, if an event split a community of species 260,000 years ago, we
might expect the divergences to be temporally clustered across multiple species
co-distributed across the barrier created by the event (the ominous "black
rectangle" :ref:`below<divergence_model_111>`).
More specifically, let's say we are interested in investigating three species
of lizards that are co-distributed across the putative barrier.
In order to infer the affect of the historical event on diversification, we
want to compare, across the three species, the timing of the divergence between
the populations on opposite sides of the putative barrier.
If the historical event caused divergence, we would expect that each of the
three pairs of lizard populations (or some subset of them) diverged around the
same time, as shown in :ref:`the figure below<divergence_model_111>`.


We can think of this as a particular *divergence model* where all three pairs
of populations share the same divergence-time parameter.
If we give the divergence-time parameter the index "1", we can use the notation
:math:`\divtimesets{1} = 111` to show that this divergence model assigns
population pairs 1, 2, and 3 to divergence-time parameter 1.
However, this is only one possible divergence model, and happens to be the most
constrained.
With three population pairs, there are 4 other possible models of divergence (5
total possible models).
Three of these models have two divergence-time parameters.
We can assign population-pair 1 to a second divergence-time parameter to get
divergence model :math:`\divtimesets{2} = 211`, as shown in
:ref:`the figure below <divergence_model_211>`.


We can also assign population-pair 2 to divergence-time parameter 2 to get
divergence model :math:`\divtimesets{3} = 121`, as :ref:`shown
below<divergence_model_121>`.

   
   A cartoon showing population-pair 1 assigned to divergence-time parameter 2,
   and population-pairs 2 and 3 assigned to divergence-time parameter 1.

And for the last possible divergence model with two divergence-time parameters,
we assign population-pair 3 to divergence-time parameter 2 to get divergence
model :math:`\divtimesets{4} = 112`, as shown in :ref:`the figure
below<divergence_model_112>`.

.. _divergence_model_121:

.. figure:: /_static/div-model-131-labels.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: divergence model 131
   
   A cartoon showing population-pair 2 assigned to divergence-time parameter 2,
   and population-pairs 1 and 3 assigned to divergence-time parameter 1.

Finally, we can add a third divergence-time parameter so that each pair of
populations is assigned to its own divergence-time parameter (divergence model
:math:`\divtimesets{5} = 123`), as shown in :ref:`the last divergence-model
figure<divergence_model_213>`.
This is the most general model of divergence, and has no co-divergence among
taxa.
Biogeographically, we can think of each free divergence-time parameter
as a "divergence event" during which one or more pairs of populations
can diverge.

.. _divergence_model_112:

.. figure:: /_static/div-model-113-labels.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: divergence model 113
   
   A cartoon showing population-pair 3 assigned to divergence-time parameter 2,
   and population-pairs 1 and 2 assigned to divergence-time parameter 1.

Being energetic herpetologists, we go out and sample individuals from each of
the lizard populations, and from those individuals collect DNA sequence data
from one or more orthologous loci per pair of populations.
We know that the sequences of a locus are related by a genealogy,
and that the shape of this genealogy is governed by demographic processes.
We also know that the genetic variation we see in the data accumulated as the
sequences evolved via mutational processes along the genealogy.
We can modify our cartoon of model :math:`\divtimesets{5} = 123` to better represent this,
as I try to do in .

.. _divergence_model_213:

.. figure:: /_static/div-model-213-labels.svg
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: divergence model 213
   
   A cartoon showing the most general model of divergence where all three
   pairs of lizard populations diverge at unique times.


.. _demog_model_cartoon:

.. figure:: /_static/demog-island-cartoon-event-labels.png
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: demographic model cartoon
   
   A cartoon showing three pairs of lizard populations that co-diverge due to
   an event 260,000 years ago.

.. _mixed_model_cartoon:

.. figure:: /_static/mixed-island-cartoon-event-labels.png
   :align: center
   :width: 600 px
   :figwidth: 90 %
   :alt: demographic model cartoon
   
   A cartoon showing three pairs of lizard populations that co-diverge due to
   an event 260,000 years ago.


Before we go any further, let's clarify some terminology that will be
used throughout the |eco|_ documentation:

.. admonition:: Definitions
    :class: definitions

    Taxon
        A pair of populations that diverged in the past. We are interested in
        comparing the timing of this divergence to other pairs of populations.
        I will use taxon interchangeably with *species* and *population pair*.

    Divergence event
        Synonymous with *divergence-time parameter*. It is a parameter of a
        *divergence model* that represents a time point in the past at which
        one or more of the taxa diverged.

    Divergence model
        A particular assignment (set partiton) of taxa to divergence-time
        parameter(s). It can range from all taxa being assigned to a single
        divergence-time parameter (i.e., "simultaneous" divergence) to each
        taxon being assigned to a unique divergence-time parameter (i.e., no
        co-divergence). Sometimes I get sloppy and just use *model*.


Next, let's jump to the ":ref:`bayesian_divergence_model_choice`" section to
see how we can use the information in the sequence data to infer the temporal
distribution of the population divergences across the three lizard species.


.. _bayesian_divergence_model_choice:

********************************
Bayesian divergence-model choice
********************************

blah

.. _prior_on_divergence_models:

**************************
Prior on divergence models
**************************

In addition to placing priors on all of the parameters of the divergence
models, we also have to place a prior on the divergence models themselves.
This can be a bit tricky, because there can be *a lot* of divergence models.
In our example of :math:`\ncomparisons{} = 3` lizard species above, we saw there were
five possible models of divergence (i.e., there were five possible ways to
assign the three species to divergence-time parameters):
There was only one way to assign the species to both one and three divergence
events,
and there were three ways to assign the three species to two divergence events.
More generally, the number of ways to assign :math:`\ncomparisons{}` taxa to
:math:`n` divergence events is the
`Stirling number of the second kind
<http://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind>`_.
Taking this a step further, there can be anywhere from :math:`1` to
:math:`\ncomparisons{}` divergence events, and so to calculate
the total number of possible divergence models, we need to calculate
the
`Stirling number of the second kind
<http://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind>`_
for :math:`1,2,\ldots,\ncomparisons{}` divergence events and sum them all up
(this is the `Bell number <http://en.wikipedia.org/wiki/Bell_number>`_
:cite:`Bell1934`).
For 3, 5, 10, and 20 taxa, there are 5, 52, 115975, and 51724158235372
possible divergence models, respectively.
The number of possible models quickly explodes as we compare more taxa!
So, how do we put a prior on all of them?!

.. _dpp:

Dirichlet-process prior on divergence models
============================================

|eco|_ takes a non-parametric approach to this problem, and treats the
number of divergence events, and the assignment of the taxa to the events, as a
Dirichlet process :cite:`Ferguson1973`.
This assigns prior probabilities directly to the divergence *models* and avoids
the combinatorial problem created when assigning prior probabilities to the 
*number* of events (Figure ).
Also, the "clumpiness" of the Dirichlet process is controlled by a
concentration parameter (:math:`\alpha`), which makes it a very flexible prior
to use for divergence models (I.e., we can control how much co-divergence we
expect across taxa *a priori*).

The basic idea of the Dirichlet process is quite simple; you assign
random variables (divergence times of our population pairs) to categories
(divergence events) one at a time following a very simple rule. When assigning
the :math:`n^{th}` random variable, you assign it to its own category (i.e.,
a new category) with probability

.. math::
    :label: dppnewcat

    \frac{\alpha}{\alpha + n -1}

or you assign it to an existing category :math:`x` with probability

.. math::
    :label: dppexistingcat

    \frac{n_x}{\alpha + n -1}

where :math:`n_x` is the number of random variables already assigned to
category :math:`x`.
OK, that might not sound very simple, but it is if we just walk through
an example using our three lizard species.
First, we have to assign our first lizard specie ("A") to a
divergence event with probability 1.0 (the species had to diverge
sometime!); let's call this the "blue" divergence event.
Next we assign the second species ("B") to either a new ("red") divergence
event with probability :math:`\alpha/\alpha + 1` or to the same "blue"
divergence event as the first species with probability :math:`1/\alpha + 1`.
For this example, let's say it gets assigned to the "blue" event.
Lastly, we assign the third species ("C") to either a new ("red") divergence
event with probability :math:`\alpha/\alpha + 2` or to the same "blue"
divergence event as the first two species with probability :math:`2/\alpha +
2`.

If we draw out all possible assignments as a tree, we get Figure dpp_tree_
below. You can adjust the concentration parameter to get a feel for how it
affects the prior probability of each divergence model. Notice that as the
concentration parameter increases we place more and more probability on the
divergence models with less clustering (less shared divergences), whereas we
place more and more probability on clustered models (shared divergences) as we
decrease the concentration parameter.

.. _dpp_tree:
.. raw:: html

    <div id="dpp_div" name="dpp_div" align="center">
        <form id="dpp_3_form" name="dpp_3_form">
            <label>Concentration parameter: </label>
            <input id="cparam" type="number" label="label" name="concentration_param" min="0.0" max="100000.0" step="any" value="1.5" onkeypress="parse_key_press(event, 3);"></input>
            <input id="update_3_button" type="button" value="Update" onclick="update_dpp_tree(3,'../_static/dpp-3-example-blank.png',1.0);"></input>
            <input type="text" name="StackOverflow1370021" value="Fix IE bug" style="display: none;"></input>
        </form>
        <canvas id="dpp_3_canvas" width="600" height="450" style="border: 1px solid rgb(211, 211, 211); align:center;"></canvas>
        <script type="text/javascript">
            update_dpp_tree(3, "../_static/dpp-3-example-blank.png", "1");
        </script>
    </div>
