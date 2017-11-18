# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
SPHINXPROJ    = Ecoevolity
SOURCEDIR     = .
# Updating path to build html in neighboring directory, within which there is a
# clone of the gh-pages branch of *this* repo as the 'html' subdirectory.
# Keeping the gh-pages branch outside so that it is easy to compile the
# documentation html there and commit and push it.
# Idea from https://daler.github.io/sphinxdoc-test/includeme.html
# BUILDDIR      = _build
BUILDDIR      = ../../ecoevolity-docs

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
