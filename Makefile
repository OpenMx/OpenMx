REXEC = R
RCOMMAND = CMD
RBUILD = build
RINSTALL = INSTALL
BUILDARGS = --dsym  # need for MacOS debug symbols
RCHECK = check
RPDF = Rd2pdf
BUILDPRE = 999.0.0
BUILDNO = $(shell ./inst/tools/findBuildNum.sh)
TARGET = OpenMx_$(BUILDPRE)-$(BUILDNO).tar.gz 
PDFFILE = $(RBUILD)/OpenMx.pdf
DOCTESTGEN = inst/tools/docTestGenerator.sh
DOCTESTFILE = inst/tools/testDocs.R
ifdef CPUS
   # snowfall
   TESTFILE = inst/tools/parallelTestModels.R
else
   TESTFILE = inst/tools/testModels.R
endif
NIGHTLYPPMLFILE = inst/tools/testNightlyPPML.R
RPROFTESTFILE = inst/tools/rprofTestModels.R
FAILTESTFILE = inst/tools/failTestModels.R
MEMORYTESTFILE = inst/tools/memoryTestModels.sh

INSTALLMAKEFLAGS=""
#INSTALLMAKEFLAGS="-j 8"   #much faster compiles

# subdirectories
RSOURCE = R
RDOCUMENTS = man
RDATA = data

# file types
RFILES = $(wildcard R/*.R)

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo ""	
	@echo "BUILDS"
	@echo ""
	@echo "  build         create an OpenMx binary for unix systems (no cross-compilation)"
	@echo "  build32       create an OpenMx binary for i386 systems"
	@echo "  build64       create an OpenMx binary for x86_64 systems"
	@echo "  buildppc      create an OpenMx binary for ppc systems"
	@echo "  srcbuild      create an OpenMx source release (alias for 'internal-build')"
	@echo "  winbuild      create an OpenMx binary on windows systems (no cross-compilation)"
	@echo "  winbuild-biarch  create an OpenMx binary for [32|64] bit windows systems"
	@echo "  cran-check    build OpenMx without NPSOL and run CRAN check"
	@echo "  cran-build    build OpenMx without NPSOL"
	@echo "  cran-winbuild build OpenMx without NPSOL and build a binary package"
	@echo ""		
	@echo "INSTALL"
	@echo ""	
	@echo "  install       build and install OpenMx on this machine"
	@echo "                (On unix, CFLAGS can be overridden in /etc/R/Makeconf)"
	@echo ""
	@echo "DOCUMENTATION"
	@echo ""	
	@echo "  pdf           create a pdf file (in build) of the OpenMx R documentation"
	@echo "  html          create Sphinx documentation (in docs/build/html) in html format"
	@echo ""
	@echo "TESTING"
	@echo ""	
	@echo "  test          run the test suite"
	@echo "  torture       run the test suite with gctorture(TRUE)"
	@echo "  check         run the R package checking system on the OpenMx package"		
	@echo "  nightly       run the nightly test suite"			
	@echo "  nightlyPPML   run the nightly test suite with PPML"			
	@echo "  testdocs      test the examples in the Sphinx documentation"	
	@echo "  failtest      run the failing test suite"
	@echo "  memorytest    run the test suite under the Valgrind memory debugger"
	@echo "  rproftest     run the test suite under the Rprof R profiler"
	@echo ""
	@echo "CLEANING"
	@echo ""	
	@echo "  clean      remove all files from the build directory"
	@echo "  veryclean  remove all files from the build directory and all *~ files"

r-libs-user-dir:
	./inst/tools/mk-r-libs-user-dir

internal-build: build/$(TARGET)

dev-doc:
	./util/rox

build/$(TARGET): $(RFILES) src/omxSymbolTable.h src/omxSymbolTable.cpp
	@if grep Rprintf src/*.cpp; then echo "*** Rprintf is not thread-safe. Use mxLog or mxLogBig."; exit 1; fi
	@if [ `grep setFinalReturns src/*.cpp | wc -l` -gt 3 ]; then echo "*** setFinalReturns is deprecated. Use populateAttrFun or addOutput."; exit 1; fi
	rm -f inst/no-npsol
	cp DESCRIPTION DESCRIPTION.bak
	sed '/Version:/d' DESCRIPTION.bak > DESCRIPTION
	echo "Version: "$(BUILDPRE)"-"$(BUILDNO) >> DESCRIPTION	
	echo '#define HAS_NPSOL 1' > src/npsolswitch.h
	echo 'NPSOL_LIBS=$(NPSOL)' > src/Makevars.win
	cat src/Makevars.win.in >> src/Makevars.win
	cp .Rbuildignore-npsol .Rbuildignore
	mkdir -p build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) build ..
	mv DESCRIPTION.bak DESCRIPTION
	rm -f man/genericFitDependencies.Rd man/imxAddDependency.Rd man/MxAlgebraFunction.Rd \
		man/omxCheckCloseEnough.Rd \
		man/mxExpectationBA81.Rd

# Developers only. This rule is for testing builds without NPSOL. 
cran: $(RFILES) src/omxSymbolTable.h src/omxSymbolTable.cpp clean
	touch inst/no-npsol
	cp DESCRIPTION DESCRIPTION.bak
	sed '/Version:/d' DESCRIPTION.bak > DESCRIPTION
	echo "Version: "$(BUILDPRE)"-"$(BUILDNO) >> DESCRIPTION	
	echo '#define HAS_NPSOL 0' > src/npsolswitch.h
	echo 'NPSOL_LIBS=' > src/Makevars.win
	cat src/Makevars.win.in >> src/Makevars.win
	cp .Rbuildignore-cran .Rbuildignore
	./util/rox
	mkdir -p build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) build ..
	mv DESCRIPTION.bak DESCRIPTION

cran-build: cran
	rm -f man/genericFitDependencies.Rd man/imxAddDependency.Rd man/MxAlgebraFunction.Rd \
		man/omxCheckCloseEnough.Rd \
		man/mxExpectationBA81.Rd
	ls -lh $(RBUILD)/OpenMx_*.tar.gz

cran-winbuild: cran
	rm -f man/genericFitDependencies.Rd man/imxAddDependency.Rd man/MxAlgebraFunction.Rd \
		man/omxCheckCloseEnough.Rd \
		man/mxExpectationBA81.Rd
	cd $(RBUILD) && R CMD INSTALL --build OpenMx_*.tar.gz

cran-check: cran
	cd $(RBUILD) && R CMD check --as-cran OpenMx_*.tar.gz
	rm -f man/genericFitDependencies.Rd man/imxAddDependency.Rd man/MxAlgebraFunction.Rd \
		man/omxCheckCloseEnough.Rd \
		man/mxExpectationBA81.Rd

pdf:
	rm -rf $(PDFFILE); $(REXEC) $(RCOMMAND) $(RPDF) --title="OpenMx Reference Manual" --output=$(PDFFILE) .
	cd docs; make latex; cd build/latex; make all-pdf

src/omxSymbolTable.h: data/omxSymbolTable.tab inst/tools/genSymbolTableHeader.R
	$(REXEC) --slave --vanilla -f inst/tools/genSymbolTableHeader.R  > src/omxSymbolTable.h

src/omxSymbolTable.cpp: data/omxSymbolTable.tab inst/tools/genSymbolTableSource.R
	$(REXEC) --slave --vanilla -f inst/tools/genSymbolTableSource.R  > src/omxSymbolTable.cpp

html: internal-build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) --html --build $(TARGET)
	rm -f build/$(TARGET)
	cd $(RBUILD); tar -zxf *gz
	mv build/OpenMx/html build/html
	mv build/OpenMx/demo build/demo
	cp build/html/* docs/source/static/Rdoc
	cp build/demo/* docs/source/static/demo
	cd docs; make clean; make html

common-build: clean internal-build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) $(BUILDARGS) --build $(TARGET)

common-build32: clean internal-build
	cd $(RBUILD); $(REXEC) --arch i386 $(RCOMMAND) $(RINSTALL) $(BUILDARGS) --build $(TARGET)

common-build64: clean internal-build
	cd $(RBUILD); $(REXEC) --arch x86_64 $(RCOMMAND) $(RINSTALL) $(BUILDARGS) --build $(TARGET)

common-buildppc: clean internal-build
	cd $(RBUILD); $(REXEC) --arch ppc $(RCOMMAND) $(RINSTALL) $(BUILDARGS) --build $(TARGET)

post-build:
	rm -f $(RBUILD)/$(TARGET)
	cd $(RBUILD); gunzip *;\
	tar --delete --file=`ls` OpenMx/npsol;\
	gzip *.tar


build32: common-build32 post-build

build64: common-build64 post-build

buildppc: common-buildppc post-build

build: common-build post-build

srcbuild: clean internal-build

winbuild: common-build

winbuild-biarch:
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) --force-biarch --build $(TARGET)

install: clean internal-build
	cd $(RBUILD); MAKEFLAGS=$(INSTALLMAKEFLAGS) $(REXEC) $(RCOMMAND) $(RINSTALL) $(BUILDARGS) $(TARGET) 

check: internal-build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RCHECK) $(TARGET)

rproftest:
	$(REXEC) --vanilla --slave < $(RPROFTESTFILE)

testdocs:
	$(DOCTESTGEN)
	$(REXEC) --vanilla --slave < $(DOCTESTFILE)

test:
	$(REXEC) --vanilla --slave -f $(TESTFILE)

test-lisrel:
	$(REXEC) --vanilla --slave -f $(TESTFILE) --args lisrel

test-csolnp:
	IMX_OPT_ENGINE=CSOLNP $(REXEC) --vanilla --slave -f $(TESTFILE) --args csolnp

torture:
	$(REXEC) -d "gdb --batch --command util/gdb-where" --vanilla --slave -f $(TESTFILE) --args gctorture

nightly:
	$(REXEC) --vanilla --slave -f $(TESTFILE) --args nightly

nightlyPPML:
	$(REXEC) --vanilla --slave < $(NIGHTLYPPMLFILE)	

failtest:
	$(REXEC) --vanilla --slave < $(FAILTESTFILE)

memorytest:
	$(MEMORYTESTFILE)

clean:
	rm -rf $(RBUILD)/*

veryclean: clean
	find . -name "*~" -exec rm -rf '{}' \;
