# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d build/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) source

.PHONY: help clean html dirhtml pickle json htmlhelp qthelp latex changes linkcheck doctest rcheck extract

source/conf.py: ../DESCRIPTION.in source/conf.py.in ./substVersion
	./substVersion

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  dirhtml   to make HTML files named index.html in directories"
	@echo "  pickle    to make pickle files"
	@echo "  json      to make JSON files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  qthelp    to make HTML files and a qthelp project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  pdf       to make LaTeX files and create a PDF"
	@echo "  changes   to make an overview of all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  rcheck    extract and run R code"

clean:
	-rm -rf build/*
	-rm -rf source/ext/__pycache__
	-rm source/*.R source/*.Rout
	-rm source/*.omx # from checkpointing in examples
	-rm source/conf.py
	-rm source/cache/*.pdf source/cache/*.png

html-initial: rcheck
	-rm -rf build/doctrees  # latex interferes with html somehow
	# https://bitbucket.org/arjones6/sphinx-numfig/issue/7/make-html-needs-to-be-run-before-make
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) build/html

html: html-initial
	cd build/html && perl -pi -e 's,http://openmx\.ssri\.psu\.edu/svn/trunk/demo/,_static/demo/,g' *.html
	cd build/html && perl -pi -e 's,\.R">_static/demo/,.R">,g' *.html
	@echo
	@echo "Build finished. The HTML pages are in build/html."

dirhtml: rcheck
	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) build/dirhtml
	@echo
	@echo "Build finished. The HTML pages are in build/dirhtml."

pickle: rcheck
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) build/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

json: rcheck
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) build/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp: rcheck
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) build/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in build/htmlhelp."

qthelp: rcheck
	$(SPHINXBUILD) -b qthelp $(ALLSPHINXOPTS) build/qthelp
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in build/qthelp, like this:"
	@echo "# qcollectiongenerator build/qthelp/OpenMx.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile build/qthelp/OpenMx.qhc"

latex: rcheck
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) build/latex
	@echo
	@echo "Build finished; the LaTeX files are in build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

pdf: rcheck
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) build/latex
	cd build/latex && make all-pdf

changes:
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) build/changes
	@echo
	@echo "The overview file is in build/changes."

rcheck: source/conf.py
	./rcheck

extract: source/conf.py
	./rcheck extract

linkcheck:
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) build/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in build/linkcheck/output.txt."

# doctest:
# 	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) build/doctest
# 	@echo "Testing of doctests in the sources finished, look at the " \
# 	      "results in build/doctest/output.txt."
