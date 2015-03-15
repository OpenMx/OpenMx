OPENMP = yes
export OPENMP
REXEC = R
export REXEC

# --dsym is need for MacOS debug symbols
# --force-biarch is for Windows 64/32 fat binary packages
BUILDARGS = --force-biarch --dsym

VERSION = $(shell ./inst/tools/findBuildNum.sh)

TARGET = OpenMx_$(VERSION).tar.gz 
PDFFILE = build/OpenMx.pdf
DOCTESTGEN = inst/tools/docTestGenerator.sh
DOCTESTFILE = inst/tools/testDocs.R
ifdef CPUS
   # snowfall
   TESTFILE = inst/tools/parallelTestModels.R
else
   TESTFILE = inst/tools/testModels.R
endif
RPROFTESTFILE = inst/tools/rprofTestModels.R
FAILTESTFILE = inst/tools/failTestModels.R
MEMORYTESTFILE = inst/tools/memoryTestModels.sh
GDBWRAP = $(shell if which gdb >/dev/null; then echo '-d gdb --debugger-args="--nx --batch --return-child-result --command util/gdb-where"'; fi)

#INSTALLMAKEFLAGS="--debug=b"   #debug dependencies
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
	@echo "  build-simple  create an OpenMx binary for unix systems without OpenMP"
	@echo "  srcbuild      create an OpenMx source release"
	@echo "  cran-build    build OpenMx without NPSOL"
	@echo "  cran-check    build OpenMx without NPSOL and run CRAN check"
	@echo ""		
	@echo "INSTALL"
	@echo ""	
	@echo "  install       install OpenMx with NPSOL"
	@echo "  cran-install  install OpenMx without NPSOL"
	@echo ""
	@echo "DOCUMENTATION"
	@echo ""	
	@echo "  pdf           create a pdf file (in build) of the OpenMx R documentation"
	@echo "  html          create Sphinx documentation (in docs/build/html) in html format"
	@echo "  doc.tar.bz2   create doc tarball suitable for our website"
	@echo ""
	@echo "TESTING"
	@echo ""	
	@echo "  test          run the test suite"
	@echo "  test-failing  run the failing test collection"
	@echo "  torture       run the test suite with gctorture(TRUE)"
	@echo "  check         run the R package checking system on the OpenMx package"		
	@echo "  nightly       run the nightly test suite"			
	@echo "  testdocs      test the examples in the Sphinx documentation"	
	@echo "  failtest      run the failing test suite"
	@echo "  memorytest    run the test suite under the Valgrind memory debugger"
	@echo "  rproftest     run the test suite under the Rprof R profiler"
	@echo ""
	@echo "CLEANING"
	@echo ""	
	@echo "  clean      remove all files from the build directory"
	@echo "  veryclean  remove all files from the build directory and all *~ files"
	@echo "  autodep    regenerate src/autodep"

r-libs-user-dir:
	./inst/tools/mk-r-libs-user-dir

code-style: $(RFILES)
	@echo "Checking code style"
	@if grep Rf_unprotect src/*.cpp; then echo "*** Rf_unprotect is error prone. Use ScopedProtect instead."; exit 1; fi
	@if grep UNPROTECT src/*.cpp; then echo "*** UNPROTECT is error prone. Use ScopedProtect instead."; exit 1; fi
	@if grep Rprintf src/*.cpp; then echo "*** Rprintf is not thread-safe. Use mxLog or mxLogBig."; exit 1; fi
	@if [ `grep strncmp src/*.cpp | wc -l` -gt 0 ]; then echo "*** Use strEQ instead of strncmp."; exit 1; fi
	@if [ `grep setFinalReturns src/*.cpp | wc -l` -gt 2 ]; then echo "*** setFinalReturns is deprecated. Use populateAttrFun or addOutput."; exit 1; fi
	@if grep --color=always --exclude '*.rda' --exclude '.R*' -r "@" demo models; then echo '*** Access of @ slots must be done using $$'; fi

build-prep:
	-[ -d build ] && rm -r ./build
	mkdir build
	git archive --format=tar HEAD | (cd build; tar -xf -)

cran-build: build-prep
	cd build && ./util/prep cran && $(REXEC) CMD build .

build: build-prep
	cd build && ./util/prep npsol && $(REXEC) CMD INSTALL $(BUILDARGS) --build .

build-simple: build-prep
	cd build && ./util/prep npsol && OPENMP=no $(REXEC) CMD INSTALL $(BUILDARGS) --build .

srcbuild: build-prep
	cd build && ./util/prep npsol && $(REXEC) CMD build .

cran-check: cran-build
	$(REXEC) CMD check build/OpenMx_*.tar.gz | tee cran-check.log
	wc -l OpenMx.Rcheck/00check.log
	@if [ $$(wc -l OpenMx.Rcheck/00check.log | cut -d ' ' -f 1) -gt 105 ]; then echo "CRAN check problems have grown; see cran-check.log" ; false; fi

pdf:
	./util/prep npsol
	rm -f $(PDFFILE); $(REXEC) CMD Rd2pdf --title="OpenMx Reference Manual" --output=$(PDFFILE) .
	cd docs; make pdf

html:
	./util/prep npsol
	cd build && R CMD INSTALL --html --no-libs --no-test-load --build ..
	cd build && tar -zxf *gz
	mv build/OpenMx/html/* docs/source/static/Rdoc
	mv build/OpenMx/demo/* docs/source/static/demo
	cd docs && make html
	cd docs/build/html && perl -pi -e 's,http://openmx\.psyc\.virginia.edu/svn/trunk/demo/,_static/demo/,g' *.html
	cd docs/build/html && perl -pi -e 's,\.R">_static/demo/,.R">,g' *.html

doc.tar.bz2: html pdf
	-rm -r build/$(VERSION)
	mkdir -p build/$(VERSION)
	mv docs/build/html/* build/$(VERSION)
	mv docs/build/latex/OpenMx.pdf build/$(VERSION)/OpenMxUserGuide.pdf
	mv build/OpenMx.pdf build/$(VERSION)
	cd build && tar jcf ../doc.tar.bz2 $(VERSION)
	-rm -r build/$(VERSION)

install: code-style
	./util/prep npsol
	MAKEFLAGS="$(INSTALLMAKEFLAGS)" $(REXEC) CMD INSTALL $(BUILDARGS) .

cran-install: code-style
	./util/prep cran
	MAKEFLAGS="$(INSTALLMAKEFLAGS)" $(REXEC) CMD INSTALL $(BUILDARGS) .

rproftest:
	$(REXEC) --vanilla --slave < $(RPROFTESTFILE)

testdocs:
	$(DOCTESTGEN)
	$(REXEC) --vanilla --slave < $(DOCTESTFILE)

test:
	$(REXEC) $(GDBWRAP) --vanilla --slave -f $(TESTFILE)

test-failing:
	$(REXEC) $(GDBWRAP) --vanilla --slave -f $(TESTFILE) --args failing

test-lisrel:
	$(REXEC) $(GDBWRAP) --vanilla --slave -f $(TESTFILE) --args lisrel

test-csolnp:
	IMX_OPT_ENGINE=CSOLNP $(REXEC) $(GDBWRAP) --vanilla --slave -f $(TESTFILE) --args csolnp

torture:
	$(REXEC) $(GDBWRAP) --vanilla --slave -f $(TESTFILE) --args gctorture

nightly:
	$(REXEC) $(GDBWRAP)  --vanilla --slave -f $(TESTFILE) --args nightly

failtest:
	$(REXEC) --vanilla --slave < $(FAILTESTFILE)

memorytest:
	$(MEMORYTESTFILE)

autodep:
	@echo "WARNING: These dependencies are not exact because they don't consider #defined CPP macros."
	cd src && gcc -MM *.cpp | perl -pe 's,[^\s]+/RcppEigen/[^\s]+, ,g' | perl -pe 's,^\s*\\\n,,' | perl -pe 's,:,: Makevars.in,' > autodep

clean:
	mkdir -p build
	-rm build/OpenMx_*.tar.gz
	-rm src/*.o
	-rm src/*.so
	-rm src/*.dll
	-rm DESCRIPTION
	-rm runtimes.csv

veryclean: clean
	-find . -name "*~" -exec rm -f '{}' \;
	-rm src/omxSymbolTable.*
	-rm src/libnpsol.a
	-grep -l 'Generated by roxygen2' man/*.Rd | xargs rm -f
