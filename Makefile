ifeq ($(OPENMP),)
  OPENMP = yes
endif
export OPENMP

ifeq ($(REXEC),)
  REXEC = R
endif
export REXEC

ifeq ($(NOT_CRAN),)
  NOT_CRAN = true
endif
export NOT_CRAN

# --dsym is need for MacOS debug symbols
# --force-biarch is for Windows 64/32 fat binary packages
BUILDARGS = --force-biarch --dsym

VERSION = $(shell ./inst/tools/findBuildNum.sh)

TARGET = OpenMx_$(VERSION).tar.gz
PDFFILE = staging/OpenMx.pdf
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

# Possibly causing grief on the UVa compute cluster
#GDBWRAP = $(shell if which gdb >/dev/null; then echo '-d gdb --debugger-args="--nx --batch --return-child-result --command util/gdb-where"'; fi)
GDBWRAP=

#MAKEFLAGS="--debug=b"   #debug dependencies
#MAKEFLAGS="-j 8"   #much faster compiles

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
	@echo "  packages-help print some hints for generating PACKAGES files"
	@echo ""
	@echo "INSTALL"
	@echo ""
	@echo "  r-libs-user-dir create R_LIBS_USER to contain installed packages"
	@echo "  cran-install    install OpenMx without NPSOL"
	@echo "  install         install OpenMx with NPSOL"
	@echo "  gpu-install     install OpenMx with CUDA (requires CUDA toolkit)"
	@echo ""
	@echo "DOCUMENTATION"
	@echo ""
	@echo "  pdf           create a pdf file (in staging) of the OpenMx R documentation"
	@echo "  html          create Sphinx documentation (in docs/build/html) in html format"
	@echo "  doc.tar.bz2   create doc tarball suitable for our website"
	@echo ""
	@echo "TESTING"
	@echo ""
	@echo "  test               run the test suite"
	@echo "  cran-check         build OpenMx without NPSOL and run CRAN check"
	@echo ""
	@echo "  test-failing  run the failing test collection"
	@echo "  torture       run the test suite with gctorture(TRUE)"
	@echo "  nightly       run the nightly test suite"
	@echo "  testdocs      test the examples in the Sphinx documentation"
	@echo "  failtest      run the failing test suite"
	@echo "  memorytest    run the test suite under the Valgrind memory debugger"
	@echo "  rproftest     run the test suite under the Rprof R profiler"
	@echo ""
	@echo "CLEANING"
	@echo ""
	@echo "  clean      remove all files from the staging directory"
	@echo "  veryclean  remove all files from the staging directory and all *~ files"
	@echo "  autodep    regenerate src/autodep"
	@echo ""
	@echo "For extra compiler diagnostics, touch ./.devmode"

r-libs-user-dir:
	./inst/tools/mk-r-libs-user-dir

code-style: $(RFILES)
	@echo "Checking code style"
	@if [ `grep R_CheckUserInterrupt src/*.cpp | wc -l` -gt 1 ]; then echo "*** omxGlobal::interrupted instead of R_CheckUserInterrupt."; exit 1; fi
	@if grep Rf_unprotect src/*.cpp; then echo "*** Rf_unprotect is error prone. Use ProtectedSEXP or Rcpp instead."; exit 1; fi
	@if grep UNPROTECT src/*.cpp; then echo "*** UNPROTECT is error prone. Use ProtectedSEXP or Rcpp instead."; exit 1; fi
	@if [ `grep Rf_error src/*.cpp | wc -l` -gt 10 ]; then echo "*** Use mxThrow instead of Rf_error."; exit 1; fi
	@if grep Rprintf src/*.cpp; then echo "*** Rprintf is not thread-safe. Use mxLog or mxLogBig."; exit 1; fi
	@if [ `grep strncmp src/*.cpp | wc -l` -gt 0 ]; then echo "*** Use strEQ instead of strncmp."; exit 1; fi
	@if [ `grep globalenv R/*.R | wc -l` -gt 6 ]; then echo "*** globalenv() interferes with testthat and is not allowed"; exit 1; fi
	@if [ `grep GetRNGstate src/*.cpp | wc -l` -gt 1 ]; then echo "*** Use BorrowRNGState instead of GetRNGstate."; exit 1; fi
	@if grep --color=always --exclude '*.rda' --exclude '*.RData' --exclude '.R*' --exclude '*.pdf' --exclude MatrixErrorDetection.R -r "@" demo inst/models; then echo '*** Access of @ slots must be done using $$'; fi

staging-clean:
	-[ -d staging ] && rm -r ./staging
	mkdir staging

staging-prep: staging-clean
	@if [ $$(git status --short --untracked-files=no 2> /dev/null | wc -l) != 0 ]; then \
	  echo '***'; echo "*** UNCOMMITTED CHANGES IGNORED ***"; \
	  echo '***'; echo "*** Use 'git diff' to see what is uncommitted"; \
          echo '***'; fi
	git archive --format=tar HEAD | (cd staging; tar -xf -)

cran-build: staging-prep
	+cd staging && sh ./util/prep cran build && $(REXEC) CMD build --resave-data .

build: staging-prep
	+cd staging && sh ./util/prep npsol build && $(REXEC) CMD INSTALL $(BUILDARGS) --build .

build-simple: staging-prep
	+cd staging && sh ./util/prep npsol build && OPENMP=no $(REXEC) CMD INSTALL $(BUILDARGS) --build .

packages-help:
	@echo 'To generate a PACKAGES file, use:'
	@echo '  echo "library(tools); write_PACKAGES('"'.', type='source'"')" | R --vanilla'
	@echo '  echo "library(tools); write_PACKAGES('"'.', type='mac.binary'"', latestOnly=FALSE)" | R --vanilla # for OS/X'

srcbuild: staging-prep packages-help
	+cd staging && sh ./util/prep npsol build && $(REXEC) CMD build --resave-data .

cran-check: cran-build
	+cd staging && _R_CHECK_FORCE_SUGGESTS_=false $(REXEC) CMD check OpenMx_*.tar.gz | tee cran-check.log
	wc -l staging/OpenMx.Rcheck/00check.log
	@if [ $$(wc -l staging/OpenMx.Rcheck/00check.log | cut -d ' ' -f 1) -gt 82 ]; then echo "CRAN check problems have grown; see cran-check.log" ; false; fi

roxygen:
	sh ./util/rox

docprep: roxygen staging-clean

pdf: docprep
	sh ./util/prep cran install
	rm -f $(PDFFILE); $(REXEC) CMD Rd2pdf --title="OpenMx Reference Manual" --output=$(PDFFILE) .
	cd docs; make pdf

html: docprep
	sh ./util/prep cran install
	cd staging && R CMD INSTALL --library=/tmp --html --no-libs --no-test-load --build ..
	cd staging && tar -zxf *gz
	mv staging/OpenMx/html/* docs/source/static/Rdoc
	mv staging/OpenMx/demo/* docs/source/static/demo
	cd docs && make html

doc.tar.bz2:
	cd docs && make clean
	$(MAKE) -j1 html pdf
	-rm -r staging/$(VERSION)
	mkdir -p staging/$(VERSION)
	mv docs/build/html/* staging/$(VERSION)
	mv docs/build/latex/OpenMx.pdf staging/$(VERSION)/OpenMxUserGuide.pdf
	mv staging/OpenMx.pdf staging/$(VERSION)
	cd staging && tar jcf ../doc.tar.bz2 $(VERSION)
	git checkout DESCRIPTION

install: code-style
	sh ./util/prep npsol install
	+$(REXEC) CMD INSTALL --no-test-load --with-keep.source $(BUILDARGS) . ;\
	git checkout DESCRIPTION

cran-install: code-style
	sh ./util/prep cran install
	+$(REXEC) CMD INSTALL --no-test-load --with-keep.source $(BUILDARGS) . ;\
	git checkout DESCRIPTION

gpu-install: code-style
	sh ./util/prep npsol install include-cuda
	+GPU_BUILD=yes $(REXEC) CMD INSTALL --no-test-load --with-keep.source $(BUILDARGS) . ;\
	git checkout DESCRIPTION

rproftest:
	$(REXEC) --vanilla --slave < $(RPROFTESTFILE)

testdocs:
	$(DOCTESTGEN)
	$(REXEC) --vanilla --slave < $(DOCTESTFILE)

test:
	sh ./tools/test && $(REXEC) $(GDBWRAP) --vanilla --slave -f $(TESTFILE)

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

nightly-parallel: nightly-csolnp nightly-slsqp  # use with -j3

nightly-csolnp:
	IMX_OPT_ENGINE=CSOLNP $(REXEC) $(GDBWRAP)  --vanilla --slave -f $(TESTFILE) --args nightly

nightly-npsol:
	IMX_OPT_ENGINE=NPSOL $(REXEC) $(GDBWRAP)  --vanilla --slave -f $(TESTFILE) --args nightly

nightly-slsqp:
	IMX_OPT_ENGINE=SLSQP $(REXEC) $(GDBWRAP)  --vanilla --slave -f $(TESTFILE) --args nightly

failtest:
	$(REXEC) --vanilla --slave < $(FAILTESTFILE)

memorytest:
	$(MEMORYTESTFILE)

autodep:
	@echo "WARNING: These dependencies are not exact because they don't consider #defined CPP macros."
	cd src && gcc \
    -MM *.cpp *.c | perl -pe 's,\S*/(R|Rcpp|BH|RcppEigen|rpf|StanHeaders)/include/\S*,,g' | perl -pe 's,^\s*\\\n,,'  |perl -pe 's,:,: Makevars,'  > autodep

clean:
	cd docs && make clean
	mkdir -p staging
	-rm staging/OpenMx_*.tar.gz
	-rm src/*.o
	-rm src/*.so
	-rm src/*.dll
	-rm runtimes.csv
	-rm src/omxSymbolTable.*
	-rm -r inst/debug
	-rm R/*.R.bak

veryclean: clean
	-rm DESCRIPTION
	-find . -name "*~" -exec rm -f '{}' \;
	-rm src/libnpsol.a
	-grep -l 'Generated by roxygen2' man/*.Rd | xargs rm -f
	-rm -r revdep/*.Rcheck
	-rm revdep/*.tar.gz
