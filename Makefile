# R package Makefile template
# Vince Buffalo <vsbuffaloAAAAA@gmail.com> (with poly-A tail removed)
#

DIR=qrqc
PKG=$(DIR)_*.tar.gz

build: clean check
	R CMD BUILD $(DIR)

check: clean
	R CMD CHECK $(DIR)

install:
	R CMD INSTALL $(PKG)

all: build install

clean: clean-pdf
	rm -f qrqc_*.tar.gz
	rm -rf qrqc/inst/doc/auto
	rm -f qrqc/src/*.o qrqc/src/*.so
	rm -rf .DS_Store
	find . | grep .DS_Store | xargs rm

clean-pdf:
	rm -f $(DIR).pdf

pdf: clean-pdf
	R CMD Rd2pdf $(DIR)