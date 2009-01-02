PACKAGE=$(wildcard *.tar.gz)

check:
	R CMD CHECK pkg

build:
	R CMD BUILD pkg

install: build
	R CMD INSTALL ${PACKAGE}

test: install
	export RCMDCHECK=FALSE;\
	cd pkg/tests;\
	R --slave < doRUnit.R

clean:
	rm -rf pkg.Rcheck
	rm -rf ${PACKAGE}
	rm pkg/inst/unitTests/report*
