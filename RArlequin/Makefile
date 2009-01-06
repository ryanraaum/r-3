PACKAGE=RArlequin

check:
	R CMD CHECK pkg

build:
	R CMD INSTALL --build pkg

install: build
	R CMD INSTALL ${PACKAGE}*.tar.gz

test: install
	export RCMDCHECK=FALSE;\
	cd pkg/tests;\
	R --slave < doRUnit.R

clean:
	rm -rf pkg.Rcheck
	rm -rf ${PACKAGE}*.tar.gz
	rm -f pkg/inst/unitTests/report*
