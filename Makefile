REXEC = R
RCOMMAND = CMD
RBUILD = build
RINSTALL = INSTALL
RCHECK = check
RPDF = Rd2dvi
TARGET = OpenMx_1.1.0-1783.tar.gz
PDFFILE = $(RBUILD)/OpenMx.pdf
DOCTESTGEN = inst/tools/docTestGenerator.sh
DOCTESTFILE = inst/tools/testDocs.R
ifdef CPUS
   TESTFILE = inst/tools/parallelTestModels.R
else
   TESTFILE = inst/tools/testModels.R
   CPUS = 1
endif
NIGHTLYFILE = inst/tools/testNightly.R
RPROFTESTFILE = inst/tools/rprofTestModels.R
FAILTESTFILE = inst/tools/failTestModels.R
MEMORYTESTFILE = inst/tools/memoryTestModels.R

# subdirectories
RSOURCE = R
RDOCUMENTS = man
RDATA = data

# file types
RDFILES =$(wildcard man/*.Rd)
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
	@echo "  winbuild32    create an OpenMx binary on 32-bit windows systems"
	@echo "  winbuild64    create an OpenMx binary on 64-bit windows systems"
	@echo ""		
	@echo "INSTALL"
	@echo ""	
	@echo "  install       build and install OpenMx on this machine"
	@echo ""
	@echo "DOCUMENTATION"
	@echo ""	
	@echo "  pdf           create a pdf file (in build) of the OpenMx R documentation"
	@echo "  html          create Sphinx documentation (in docs/build/html) in html format"
	@echo ""
	@echo "TESTING"
	@echo ""	
	@echo "  test (CPUS=n) run the test suite"	
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

internal-build: build/$(TARGET)

build/$(TARGET): $(RFILES) $(RDFILES)
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RBUILD) ..

pdf:
	rm -rf $(PDFFILE); $(REXEC) $(RCOMMAND) $(RPDF) --pdf --title="OpenMx Reference Manual" --output=$(PDFFILE) .
	cd docs; make latex; cd build/latex; make all-pdf

html: internal-build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) --html --build $(TARGET)
	rm -f build/$(TARGET)
	cd $(RBUILD); gnutar -zxf *.tgz
	mv build/OpenMx/html build/html
	mv build/OpenMx/demo build/demo
	cp build/html/* docs/source/static/Rdoc
	cp build/demo/* docs/source/static/demo
	cd docs; make clean; make html

common-build: clean internal-build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) --build $(TARGET)

common-build32: clean internal-build
	cd $(RBUILD); $(REXEC) --arch i386 $(RCOMMAND) $(RINSTALL) --build $(TARGET)

common-build64: clean internal-build
	cd $(RBUILD); $(REXEC) --arch x86_64 $(RCOMMAND) $(RINSTALL) --build $(TARGET)

common-buildppc: clean internal-build
	cd $(RBUILD); $(REXEC) --arch ppc $(RCOMMAND) $(RINSTALL) --build $(TARGET)

post-build:
	rm -f $(RBUILD)/$(TARGET)
	cd $(RBUILD); gunzip *;\
	gnutar --delete --file=`ls` OpenMx/npsol;\
	gzip *.tar


build32: common-build32 post-build

build64: common-build64 post-build

buildppc: common-buildppc post-build

build: common-build post-build

srcbuild: clean internal-build

winbuild: common-build

winbuild32: clean internal-build
	cd $(RBUILD); $(REXEC) --arch i386 $(RCOMMAND) $(RINSTALL) --build $(TARGET)

winbuild64: clean internal-build
	cd $(RBUILD); $(REXEC) --arch x64 $(RCOMMAND) $(RINSTALL) --build $(TARGET)

install: clean internal-build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) $(TARGET)

check: internal-build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RCHECK) $(TARGET)

rproftest:
	$(REXEC) --vanilla --slave < $(RPROFTESTFILE)

testdocs:
	$(DOCTESTGEN)
	$(REXEC) --vanilla --slave < $(DOCTESTFILE)

test:
	$(REXEC) --vanilla --slave --cpus=$(CPUS) < $(TESTFILE)
	
nightly:
	$(REXEC) --vanilla --slave < $(NIGHTLYFILE)	

failtest:
	$(REXEC) --vanilla --slave < $(FAILTESTFILE)

memorytest:
	$(REXEC) -d "valgrind --tool=memcheck --leak-check=full --suppressions=inst/tools/OpenMx.supp --quiet" --vanilla --slave < $(MEMORYTESTFILE)

clean:
	rm -rf $(RBUILD)/*
	rm -rf models/passing/temp-files/*
	rm -rf models/failing/temp-files/*

veryclean: clean
	find . -name "*~" -exec rm -rf '{}' \;
