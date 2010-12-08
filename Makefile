# R package Makefile template
# Vince Buffalo <vsbuffaloAAAAA@gmail.com> (with poly-A tail removed)
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