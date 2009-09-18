REXEC = R
RCOMMAND = CMD
RBUILD = build
RINSTALL = INSTALL
RCHECK = check
RPDF = Rd2dvi
TARGET = OpenMx_0.1.4-827.tar.gz
PDFFILE = $(RBUILD)/OpenMx.pdf
TESTFILE = inst/tools/testModels.R
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
	@echo "  build      create an OpenMx binary (in build) that can be exported"
	@echo "  install    build and install OpenMx on this machine"
	@echo "  pdf        create a pdf file (in build) of the OpenMx R documentation"
	@echo "  html       create Sphinx documentation (in docs/build/html) in html format"
	@echo "  test       run the OpenMx test suite"
	@echo "  failtest   run the OpenMx failing test suite"
	@echo "  memorytest run the OpenMx test suite under the Valgrind memory debugger"
	@echo "  rproftest  run the OpenMx test suite under the Rprof R profiler"
	@echo "  check      run the R package checking system on the OpenMx package"
	@echo "  clean      remove all files from the build directory"
	@echo "  veryclean  remove all files from the build directory and all *~ files"

internal-build: build/$(TARGET)

build/$(TARGET): $(RFILES) $(RDFILES)
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RBUILD) ..

pdf:
	rm -rf $(PDFFILE); $(REXEC) $(RCOMMAND) $(RPDF) --pdf --title="OpenMx Reference Manual" --output=$(PDFFILE) .
	cd docs; make latex; cd build/latex; make all-pdf

html: build
	rm -f build/$(TARGET)
	cd $(RBUILD); gunzip *.gz; tar -xf *.tar
	mv build/OpenMx/html build/html
	cp build/html/* docs/source/static/Rdoc
	cd docs; make clean; make html

build: clean internal-build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) --build $(TARGET)

install: clean internal-build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) $(TARGET)

check: internal-build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RCHECK) $(TARGET)

rproftest:
	$(REXEC) --vanilla --slave < $(RPROFTESTFILE)

test:
	$(REXEC) --vanilla --slave < $(TESTFILE)

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
