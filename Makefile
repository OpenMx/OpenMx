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
ROXDOC = man/genericFitDependencies.Rd man/imxAddDependency.Rd man/MxAlgebraFunction.Rd \
	man/omxCheckCloseEnough.Rd \
	man/mxExpectationBA81.Rd man/imxPPML.Rd man/imxPPML.Test.Battery.Rd \
	man/imxCheckVariables.Rd \
	man/imxCheckMatrices.Rd \
	man/imxCheckVariables.Rd \
	man/imxConstraintRelations.Rd \
	man/imxConvertIdentifier.Rd \
	man/imxConvertLabel.Rd \
	man/imxConvertSubstitution.Rd \
	man/imxExtractNames.Rd \
	man/imxExtractReferences.Rd \
	man/imxGenerateNamespace.Rd \
	man/imxIdentifier.Rd \
	man/imxIsDefinitionVariable.Rd \
	man/imxLocateIndex.Rd \
	man/imxLocateLabel.Rd \
	man/imxPreprocessModel.Rd \
	man/imxReverseIdentifier.Rd \
	man/imxSeparatorChar.Rd \
	man/imxVerifyName.Rd \
	man/imxVerifyReference.Rd \
	man/imxCheckMatrices.Rd \
	man/imxCheckVariables.Rd \
	man/imxConstraintRelations.Rd \
	man/imxConvertIdentifier.Rd \
	man/imxConvertLabel.Rd \
	man/imxConvertSubstitution.Rd \
	man/imxCreateMatrix.Rd \
	man/imxDataTypes.Rd \
	man/imxDeparse.Rd \
	man/imxDependentModels.Rd \
	man/imxDiff.Rd \
	man/imxDmvnorm.Rd \
	man/imxEvalByName.Rd \
	man/imxExtractMethod.Rd \
	man/imxExtractNames.Rd \
	man/imxExtractReferences.Rd \
	man/imxFilterDefinitionVariables.Rd \
	man/imxFlattenModel.Rd \
	man/imxFreezeModel.Rd \
	man/imxGenSwift.Rd \
	man/imxGenerateLabels.Rd \
	man/imxGenerateNamespace.Rd \
	man/imxGenericModelBuilder.Rd \
	man/imxHasNPSOL.Rd \
	man/imxIdentifier.Rd \
	man/imxIndependentModels.Rd \
	man/imxInitModel.Rd \
	man/imxIsDefinitionVariable.Rd \
	man/imxIsPath.Rd \
	man/imxLocateFunction.Rd \
	man/imxLocateIndex.Rd \
	man/imxLocateLabel.Rd \
	man/imxLookupSymbolTable.Rd \
	man/imxModelBuilder.Rd \
	man/imxModelTypes.Rd \
	man/imxMpiWrap.Rd \
	man/imxOriginalMx.Rd \
	man/imxPPML.Rd \
	man/imxPPML.Test.Battery.Rd \
	man/imxPreprocessModel.Rd \
	man/imxReplaceMethod.Rd \
	man/imxReplaceModels.Rd \
	man/imxReservedNames.Rd \
	man/imxReverseIdentifier.Rd \
	man/imxSameType.Rd \
	man/imxSeparatorChar.Rd \
	man/imxSfClient.Rd \
	man/imxSimpleRAMPredicate.Rd \
	man/imxSquareMatrix.Rd \
	man/imxSymmetricMatrix.Rd \
	man/imxTypeName.Rd \
	man/imxVerifyMatrix.Rd \
	man/imxVerifyModel.Rd \
	man/imxVerifyName.Rd \
	man/imxVerifyReference.Rd \
	man/MxFlatModel.Rd \
	man/MxLISRELModel.Rd \
	man/MxRAMModel.Rd \
	man/MxCharOrList-class.Rd \
	man/MxCharOrNumber-class.Rd \
	man/MxListOrNull-class.Rd \
	man/MxOptionalChar-class.Rd \
	man/MxOptionalCharOrNumber-class.Rd \
	man/MxOptionalLogical-class.Rd \
	man/MxOptionalMatrix-class.Rd \
	man/MxOptionalNumeric-class.Rd \
	man/mxComputeOnce.Rd \
	man/mxComputeIterate.Rd \
	man/mxComputeSequence.Rd \
	man/MxLISRELModel-class.Rd \
	man/MxRAMModel-class.Rd \
	man/mxComputeHessianQuality.Rd \
	man/mxComputeStandardError.Rd \
	man/mxComputeGradientDescent.Rd \
	man/mxComputeNewtonRaphson.Rd \
	man/mxComputeNumericDeriv.Rd \
	man/mxFitFunctionMultigroup.Rd \
	man/mxComputeNothing.Rd \
	man/mxComputeReportDeriv.Rd

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
	@echo "  install       install OpenMx with NPSOL"
	@echo "  cran-install  install OpenMx without NPSOL"
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
	@echo "  autodep    regenerate src/autodep"

r-libs-user-dir:
	./inst/tools/mk-r-libs-user-dir

internal-build: build/$(TARGET)

dev-doc:
	./util/rox

code-style: $(RFILES) src/omxSymbolTable.h src/omxSymbolTable.cpp
	@echo "Checking code style"
	@if grep Rprintf src/*.cpp; then echo "*** Rprintf is not thread-safe. Use mxLog or mxLogBig."; exit 1; fi
	@if [ `grep setFinalReturns src/*.cpp | wc -l` -gt 3 ]; then echo "*** setFinalReturns is deprecated. Use populateAttrFun or addOutput."; exit 1; fi

npsol-prep: code-style
	rm -f inst/no-npsol
	cp DESCRIPTION DESCRIPTION.bak
	sed '/Version:/d' DESCRIPTION.bak > DESCRIPTION
	echo "Version: "$(BUILDPRE)"-"$(BUILDNO) >> DESCRIPTION	
	echo 'NPSOL_DIR= "..\inst"' > src/Makevars.win
	echo 'NPSOL=-lnpsol$$(WIN) -L$$(NPSOL_DIR)' >> src/Makevars.win
	echo 'NPSOL_LIBS=$$(NPSOL)' >> src/Makevars.win
	cat src/Makevars.win.in >> src/Makevars.win
	echo '#define HAS_NPSOL 1' > src/npsolswitch.h
	cp .Rbuildignore-npsol .Rbuildignore
	-./util/rox

no-npsol-prep: code-style
	touch inst/no-npsol
	cp DESCRIPTION DESCRIPTION.bak
	sed '/Version:/d' DESCRIPTION.bak > DESCRIPTION
	echo "Version: "$(BUILDPRE)"-"$(BUILDNO) >> DESCRIPTION	
	echo '#define HAS_NPSOL 0' > src/npsolswitch.h
	echo 'NPSOL_LIBS=' > src/Makevars.win
	cat src/Makevars.win.in >> src/Makevars.win
	cp .Rbuildignore-cran .Rbuildignore
	-./util/rox

build/$(TARGET): npsol-prep
	mkdir -p build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) build ..
	mv DESCRIPTION.bak DESCRIPTION
	rm -f $(ROXDOC)

cran-build: no-npsol-prep
	$(REXEC) $(RCOMMAND) build .
	rm -f $(ROXDOC)

cran-winbuild: no-npsol-prep
	rm -f $(ROXDOC)
	cd $(RBUILD) && R CMD INSTALL --build OpenMx_*.tar.gz

cran-check: cran-build
	cd .. && R CMD check OpenMx_*.tar.gz
	rm -f $(ROXDOC)

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

install: npsol-prep
	MAKEFLAGS=$(INSTALLMAKEFLAGS) $(REXEC) $(RCOMMAND) $(RINSTALL) $(BUILDARGS) .

cran-install: no-npsol-prep
	MAKEFLAGS=$(INSTALLMAKEFLAGS) $(REXEC) $(RCOMMAND) $(RINSTALL) $(BUILDARGS) .

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
	$(REXEC) -d "gdb --nx --batch --return-child-result --command util/gdb-where" --vanilla --slave -f $(TESTFILE) --args gctorture

nightly:
	$(REXEC) --vanilla --slave -f $(TESTFILE) --args nightly

nightlyPPML:
	$(REXEC) --vanilla --slave < $(NIGHTLYPPMLFILE)	

failtest:
	$(REXEC) --vanilla --slave < $(FAILTESTFILE)

memorytest:
	$(MEMORYTESTFILE)

autodep:
	@echo "WARNING: These dependencies are not exact because they don't consider #defined CPP macros."
	cd src && gcc -MM *.cpp *.c | perl -pe 's,\bEigen[^\s]*\s, ,g' | perl -pe 's,^\s*\\\n,,' > autodep

clean:
	-rm $(RBUILD)/OpenMx_*.tar.gz
	-rm src/*.o
	-rm src/*.so

veryclean: clean
	find . -name "*~" -exec rm -rf '{}' \;
	-rm -f $(ROXDOC)
	-rm src/omxSymbolTable.*
	-rm src/libnpsol.a
