# Makefile for building Sphinx documentation

SPHINXBUILD   = sphinx-build
SPHINXPROJ    = CBMPy
SOURCEDIR     = source
BUILDDIR      = _build

# Targets
.PHONY: help html latex dirhtml singlehtml htmlhelp qthelp epub text clean

help:
	@echo "Please use \`make <target>' where <target> is one of:"
	@echo "  html       to make standalone HTML files"
	@echo "  latex      to make LaTeX files (can be compiled)"
	@echo "  dirhtml    to make HTML files named index.html in directories"
	@echo "  singlehtml to make a single large HTML file"
	@echo "  epub       to make an epub"
	@echo "  text       to make text files"
	@echo "  clean      to clean build directory"

html:
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

latex:
	$(SPHINXBUILD) -b latex $(SOURCEDIR) $(BUILDDIR)/latex
	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
	@echo "Run \`make all' in that directory to run these through xelatex."

dirhtml:
	$(SPHINXBUILD) -b dirhtml $(SOURCEDIR) $(BUILDDIR)/dirhtml
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."

singlehtml:
	$(SPHINXBUILD) -b singlehtml $(SOURCEDIR) $(BUILDDIR)/singlehtml
	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

htmlhelp:
	$(SPHINXBUILD) -b htmlhelp $(SOURCEDIR) $(BUILDDIR)/htmlhelp
	@echo "Build finished; the HTML help pages are in $(BUILDDIR)/htmlhelp."

epub:
	$(SPHINXBUILD) -b epub $(SOURCEDIR) $(BUILDDIR)/epub
	@echo "Build finished. The epub file is in $(BUILDDIR)/epub."

text:
	$(SPHINXBUILD) -b text $(SOURCEDIR) $(BUILDDIR)/text
	@echo "Build finished. The text files are in $(BUILDDIR)/text."

clean:
	rm -rf $(BUILDDIR)/*
