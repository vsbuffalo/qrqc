DIR=qrqc
PKG=qrqc_0.5.tar.gz

build:
	R CMD build $(DIR)

install:
	R CMD install $(PKG)

all: build install
