REXEC = R
RCOMMAND = CMD
RBUILD = build
RINSTALL = INSTALL
RCHECK = check
TARGET = OpenMx_0.1-1.tar.gz

# subdirectories
RSOURCE = R
RDOCUMENTS = man
RDATA = data

# file types
RDFILES = *.Rd
RFILES = *.R


nothing:
	@echo \
	'Please type make [build | install | check | clean | veryclean]'

build: build/$(TARGET)

build/$(TARGET): $(RSOURCE)/$(RFILES) $(RDOCUMENTS)/$(RDFILES)
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RBUILD) ..

$(RSOURCE)/$(RFILES):

$(RDOCUMENTS)/$(RDFILES):

install: build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) $(TARGET)

check: build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RCHECK) $(TARGET)

clean:
	rm -rf $(RBUILD)/*

veryclean: clean
	find . -name "*~" -exec rm -rf '{}' \;
