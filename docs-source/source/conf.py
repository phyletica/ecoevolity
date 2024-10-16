# -*- coding: utf-8 -*-

import time

# Ecoevolity documentation build configuration file, created by
# sphinx-quickstart on Thu Jun 29 21:16:35 2017.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.imgmath',
    'sphinxcontrib.bibtex',
    'sphinx.ext.githubpages',
    'sphinx_toolbox.collapse']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
source_suffix = ['.rst', '.md']
# source_suffix = '.rst'

# Add Markdown parser; need to insall recommonmark python package
source_parsers = {
        '.md': 'recommonmark.parser.CommonMarkParser',
        }

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'Ecoevolity'
copyright = u'2017, Jamie R. Oaks'
author = u'Jamie R. Oaks'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = u'0.1'
# The full version, including alpha/beta/rc tags.
release = u'0.1.0'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store',
        # 'tutorials/*.rst',
        'snippets/*.rst',
        ]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Adding explicit sidebars as per Alabaster's documentation
html_sidebars = {
        '**': [
            'about.html',
            'navigation.html',
            'relations.html',
            'searchbox.html',
            'donate.html',
        ]
}

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
        'logo': 'site-logo.svg',
        'github_user': 'phyletica',
        'github_repo': 'ecoevolity',
        'github_button': True,
        'github_banner': False,
        # 'description': 'Estimating evolutionary coevality',
        # 'show_powered_by': True,
        'fixed_sidebar': True,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static', '_static/custom.css']


# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'Ecoevolitydoc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
    "preamble": r"""
\usepackage{xspace}
\newcommand{\given}{\ensuremath{\,|\,}\xspace}
\newcommand{\pr}{\ensuremath{p}}
\newcommand{\data}{\ensuremath{D}\xspace}
\newcommand{\model}[1][]{\ensuremath{M_{#1}}\xspace}
\newcommand{\parameters}[1][]{\ensuremath{\Theta_{#1}}\xspace}
\newcommand{\parameter}[1][]{\ensuremath{\theta_{#1}}\xspace}
\newcommand{\diff}[1]{\ensuremath{\mathrm{d}#1}}
\newcommand{\ncomparisons}{\ensuremath{\mathcal{N}\xspace}}
\newcommand{\nevents}[1][]{\ensuremath{k_{#1}\xspace}}
\newcommand{\nloci}[1][]{\ensuremath{m_{#1}\xspace}}
\newcommand{\comparisondata}[1][]{\ensuremath{D_{#1}}\xspace}
\newcommand{\rgmurate}{\ensuremath{u}\xspace}
\newcommand{\grmurate}{\ensuremath{v}\xspace}
\newcommand{\murate}[1][]{\ensuremath{\mu_{#1}}\xspace}
\newcommand{\murates}[1][]{\ensuremath{\boldsymbol{\mu}_{#1}}\xspace}
\newcommand{\gfreq}[1][]{\ensuremath{\pi_{#1}}\xspace}
\newcommand{\gfreqs}[1][]{\ensuremath{\boldsymbol{\pi}_{#1}}\xspace}
\newcommand{\comparisondivtime}[1][]{\ensuremath{t_{#1}}\xspace}
\newcommand{\comparisondivtimes}[1][]{\ensuremath{\mathbf{t}_{#1}}\xspace}
\newcommand{\divtime}[1][]{\ensuremath{\tau_{#1}}\xspace}
\newcommand{\divtimes}[1][]{\ensuremath{\boldsymbol{\tau}_{#1}}\xspace}
\newcommand{\divtimemodel}[1][]{\ensuremath{T_{#1}}\xspace}
\newcommand{\divtimesets}{\ensuremath{\mathcal{T}}\xspace}
\newcommand{\comparisoneventtime}[1][]{\ensuremath{t_{#1}}\xspace}
\newcommand{\comparisoneventtimes}[1][]{\ensuremath{\mathbf{t}_{#1}}\xspace}
\newcommand{\eventtime}[1][]{\ensuremath{\tau_{#1}}\xspace}
\newcommand{\eventtimes}[1][]{\ensuremath{\boldsymbol{\tau}_{#1}}\xspace}
\newcommand{\eventtimemodel}[1][]{\ensuremath{T_{#1}}\xspace}
\newcommand{\eventtimesets}{\ensuremath{\mathcal{T}}\xspace}
\newcommand{\genetree}[1][]{\ensuremath{g_{#1}}\xspace}
\newcommand{\sptree}[1][]{\ensuremath{S_{#1}}\xspace}
\newcommand{\sptrees}[1][]{\ensuremath{\mathbf{S}_{#1}}\xspace}
\newcommand{\descendantpopindex}[1]{\ensuremath{D{#1}}}
\newcommand{\rootpopindex}[1][]{\ensuremath{R{#1}}\xspace}
\newcommand{\epopsize}[1][]{\ensuremath{N_{e}^{#1}}\xspace}
\newcommand{\comparisonpopsizes}[1][]{\ensuremath{\mathbb{N}_{e}{#1}}\xspace}
\newcommand{\collectionpopsizes}[1][]{\ensuremath{\mathbf{N_{e}}_{#1}}\xspace}
\newcommand{\rootrelativepopsize}{\ensuremath{R_{\epopsize[\rootpopindex]}}\xspace}
\newcommand{\concentration}{\ensuremath{\alpha}\xspace}
\newcommand{\basedistribution}{\ensuremath{H}\xspace}

\newcommand{\tree}{\ensuremath{T}\xspace}
\newcommand{\nTips}{\ensuremath{N}\xspace}
\newcommand{\node}[1]{\ensuremath{t_{#1}}\xspace}
\newcommand{\nodes}{\ensuremath{\boldsymbol{\node{}}}\xspace}
\newcommand{\divTimeSymbol}{\tau}
\newcommand{\divTime}[1]{\ensuremath{\divTimeSymbol_{#1}}\xspace}
\newcommand{\divTimes}{\ensuremath{\boldsymbol{\divTime{}}}\xspace}
\newcommand{\divTimeParentOf}[1]{\ensuremath{y({#1})}\xspace}
\newcommand{\probChangeDimension}{\ensuremath{\psi}\xspace}
\newcommand{\probBreakPolytomy}{\ensuremath{\Upsilon}\xspace}
\newcommand{\nWaysToBreakPolytomy}{\ensuremathmath{k_b}\xspace}
\newcommand{\bellNumber}{\ensuremath{B}\xspace}
\newcommand{\stirling}[2]{\ensuremath{S_2}(#1, #2)\xspace}
\newcommand{\treeClassPriorProb}[1]{\ensuremath{\pi_{\tree}(#1)}\xspace}
\newcommand{\maxSlide}{\ensuremath{\delta_{\divTime{}}}\xspace}
\newcommand{\uniformDeviate}{\ensuremath{u}\xspace}
\newcommand{\probLumpOverProbSplit}{\ensuremath{\gamma_S}\xspace}
\newcommand{\probLumpNeighborOverProbSplitNeighbor}{\ensuremath{\phi_S}\xspace}
\newcommand{\probSplitOverProbLump}{\ensuremath{\gamma_M}\xspace}
\newcommand{\probSplitNeighborOverProbLumpNeighbor}{\ensuremath{\phi_M}\xspace}
\newcommand{\propdens}[1]{\ensuremath{g_{#1}}\xspace}
\newcommand{\multipropdens}[1]{\ensuremath{\boldsymbol{g}_{#1}}\xspace}
\newcommand{\nOf}[2][]{\ensuremath{n_{#1}(#2)}\xspace}
\newcommand{\nDivs}{\ensuremath{\nOf[]{\divTime{}}}\xspace}
\newcommand{\nNodes}{\ensuremath{\nOf{\node{}}}\xspace}
\newcommand{\nTrees}{\ensuremath{\nOf{\tree{}}}\xspace}
\newcommand{\nSharedDivs}{\ensuremath{\nOf[s]{\divTime{}}}\xspace}
\newcommand{\nNodesMappedTo}[1]{\ensuremath{\nOf[]{\node{} \mapsto #1}}\xspace}
\newcommand{\nPolyNodesMovingTo}[1]{\ensuremath{\nOf[p]{\node{} \Rightarrow #1}}\xspace}
\newcommand{\nWaysToSplitAllPolytomies}{\ensuremath{\Phi}\xspace}
\newcommand{\probDivTimePartition}{\ensuremath{\Xi}\xspace}
\newcommand{\modelState}{\ensuremath{\Theta}\xspace}
\newcommand{\multiplier}{\ensuremath{m}\xspace}
\newcommand{\proposed}{\ensuremath{^{\prime}}\xspace}
\newcommand{\tuningparameter}{\ensuremath{\lambda}\xspace}
\newcommand{\uniformdeviate}{\ensuremath{u}\xspace}
\newcommand{\observedallelecount}[1][]{\ensuremath{n_{#1}}\xspace}
\newcommand{\observedredallelecount}[1][]{\ensuremath{r_{#1}}\xspace}
\newcommand{\nodeallelecount}[2]{\ensuremath{n_{#1}^{#2}}}
\newcommand{\noderedallelecount}[2]{\ensuremath{r_{#1}^{#2}}}
\newcommand{\allelecount}[1][]{\ensuremath{\nodeallelecount{#1}{}}\xspace}
\newcommand{\redallelecount}[1][]{\ensuremath{\noderedallelecount{#1}{}}\xspace}
\newcommand{\leafallelecounts}[1][]{\ensuremath{\mathbf{n}_{#1}}\xspace}
\newcommand{\leafredallelecounts}[1][]{\ensuremath{\mathbf{r}_{#1}}\xspace}
\newcommand{\alldata}[1][]{\ensuremath{\mathbf{D}}\xspace}
\newcommand{\allepopsizes}{\ensuremath{\boldsymbol{N_{e}}}\xspace}
\newcommand{\alphaOfDivTimeBetaPrior}{\ensuremath{\alpha_\divTimeSymbol}}
""",
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'Ecoevolity.tex', u'Ecoevolity Documentation',
     u'Jamie R. Oaks', 'manual'),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'ecoevolity', u'Ecoevolity Documentation',
     [author], 1)
]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'Ecoevolity', u'Ecoevolity Documentation',
     author, 'Ecoevolity', 'One line description of project.',
     'Miscellaneous'),
]



# -- Options for Epub output ----------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']

imgmath_image_format = "svg"
imgmath_latex_preamble = latex_elements["preamble"]

rst_epilog = """
.. |jro| replace:: Jamie Oaks
.. _jro: http://phyletica.org

.. |phyleticalab| replace:: Phyletica Lab 
.. _phyleticalab: http://phyletica.org

.. |eco| replace:: ecoevolity
.. _eco: http://phyletica.org/ecoevolity
.. |Eco| replace:: Ecoevolity
.. _Eco: http://phyletica.org/ecoevolity
.. |eco_gh| replace:: ecoevolity
.. _eco_gh: https://github.com/phyletica/ecoevolity
.. |eco_url| replace:: http://phyletica.org/ecoevolity
.. |eco_gh_url| replace:: https://github.com/phyletica/ecoevolity
.. |eco_copyright| replace:: **Copyright 2015-{this_year} Jamie R. Oaks**

.. |pyco| replace:: pycoevolity 
.. _pyco: https://github.com/phyletica/pycoevolity
.. |Pyco| replace:: Pycoevolity
.. _Pyco: https://github.com/phyletica/pycoevolity

.. |simco| replace:: ``simcoevolity``
.. |sumco| replace:: ``sumcoevolity``
.. |dpprobs| replace:: DPprobs
.. |cdpprobs| replace:: ``dpprobs``
.. |ceco| replace:: ``ecoevolity``
.. |pyco-sumchains| replace:: ``pyco-sumchains``
.. |pyco-sumtimes| replace:: ``pyco-sumtimes``
.. |pyco-sumevents| replace:: ``pyco-sumtimes``

.. |phyco| replace:: phycoeval
.. _phyco: http://phyletica.org/ecoevolity
.. |Phyco| replace:: Phycoeval
.. _Phyco: http://phyletica.org/ecoevolity
.. |phyco_gh| replace:: phycoeval
.. _phyco_gh: https://github.com/phyletica/ecoevolity
.. |phyco_url| replace:: http://phyletica.org/ecoevolity
.. |phyco_gh_url| replace:: https://github.com/phyletica/ecoevolity
.. |phyco_copyright| replace:: **Copyright 2021-{this_year} Jamie R. Oaks**
.. |cphyco| replace:: ``phycoeval`` 
.. |simphyco| replace:: simphycoeval
.. |Simphyco| replace:: Simphycoeval
.. |csimphyco| replace:: ``simphycoeval``
.. |sumphyco| replace:: sumphycoeval
.. |Sumphyco| replace:: Sumphycoeval
.. |csumphyco| replace:: ``sumphycoeval``

.. |eco_logo_long| image:: /_static/ecoevolity-logo.svg
                    :width: 95%
                    :alt: Ecoevolity
                    :target: index.html

.. |phyco_logo_long| image:: /_static/phycoeval-logo-long.svg
                    :width: 95%
                    :alt: Phycoeval
                    :target: index.html

.. |Tracer| replace:: Tracer
.. _Tracer: https://github.com/beast-dev/tracer/releases
.. |Figtree| replace:: Figtree
.. _Figtree: https://github.com/rambaut/figtree/releases

.. |git| replace:: Git
.. _git: http://git-scm.com/

.. |yaml| replace:: YAML 
.. _yaml: http://yaml.org/
.. |yamllint| replace:: http://www.yamllint.com/

.. |icytree| replace:: IcyTree 
.. _icytree: https://icytree.org

.. |gpl3| replace:: http://www.gnu.org/licenses/gpl-3.0-standalone.html

.. |cmake| replace:: CMake
.. _cmake: https://cmake.org/
""".format(this_year = time.strftime('%Y'))
