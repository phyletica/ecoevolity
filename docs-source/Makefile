# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
# Calling python to build can be very useful, for example, if we want to build
# the documentation within a particular Python environment, because
# sphinx-build will use the "globally" installed Python.
# SPHINXBUILD   = sphinx-build
SPHINXBUILD   = python -msphinx
SOURCEDIR     = source
BUILDDIR      = build
PUBLISHDIR	  = ../docs

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: publish help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

publish: | clean html
	@echo ""
	@echo "Copying files to published directory..."
	@echo ""
	@mkdir -p "$(PUBLISHDIR)" \
		&& rsync -ap --progress --delete $(BUILDDIR)/html/. "$(PUBLISHDIR)"
	@echo ""
	@echo "Committing and pushing all changes to source and published html files..."
	@echo ""
	@git add "$(SOURCEDIR)" "$(PUBLISHDIR)" \
		&& git commit \
			-m "Updating documentation at $$(date -u +"%Y-%m-%d %H:%M:%S")" \
			-m "This commit was automated by the Sphinx documentation make file." \
			-m "All changes to source and built/published html files are committed." \
		&& git push
