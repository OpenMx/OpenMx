REXEC = R
RCOMMAND = CMD
RBUILD = build
RINSTALL = INSTALL
RCHECK = check
RPDF = Rd2dvi
TARGET = OpenMx_0.1-1.tar.gz
DOCFILE = $(RBUILD)/OpenMx.pdf
TESTFILE = inst/tools/testModels.R

# subdirectories
RSOURCE = R
RDOCUMENTS = man
RDATA = data

# file types
RDFILES = *.Rd
RFILES = *.R


nothing:
	@echo \
	'Please type make [build | install | doc | check | clean | veryclean]'

internal-build: build/$(TARGET)

build/$(TARGET): $(RSOURCE)/$(RFILES) $(RDOCUMENTS)/$(RDFILES)
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RBUILD) ..

doc:
	rm -rf $(DOCFILE); $(REXEC) $(RCOMMAND) $(RPDF) --pdf --title="OpenMx Reference Manual" --output=$(DOCFILE) .

$(RSOURCE)/$(RFILES):

$(RDOCUMENTS)/$(RDFILES):

build: clean internal-build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) --build $(TARGET)

install: clean internal-build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) $(TARGET)

check: internal-build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RCHECK) $(TARGET)

test:
	$(REXEC) --vanilla --slave < $(TESTFILE)

clean:
	rm -rf $(RBUILD)/*

veryclean: clean
	find . -name "*~" -exec rm -rf '{}' \;
