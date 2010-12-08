# R package Makefile template
# Vince Buffalo <vsbuffaloAAAAA@gmail.com> (with poly-A tail removed)
#
# Replace 'package' with the name of your package (which should be a
# directory). My layout for package development is:
# package_build
#  - TODO
#  - Readme.md (for Github)
#  - package (actual R package directory)
#  - Makefile (this)
#  - test-all.R (but this is just a wrapper around package/tests
#

DIR=qrqc
PKG=$(DIR)_*.tar.gz

build: clean
	R CMD build $(DIR)

install:
	R CMD install $(PKG)

all: build install

clean: clean-pdf
	rm -f qrqc_*.tar.gz

clean-pdf:
	rm -f $(DIR).pdf

pdf: clean-pdf
	R CMD Rd2pdf $(DIR)