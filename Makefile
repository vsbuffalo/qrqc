DIR=qrqc
PKG=qrqc_0.9.tar.gz

build:
	R CMD build $(DIR)

install:
	R CMD install $(PKG)

all: build install

clean-pdf:
	rm $(DIR).pdf

pdf: clean-pdf
	R CMD Rd2pdf $(DIR)