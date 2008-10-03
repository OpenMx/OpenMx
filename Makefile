REXEC = R
RCOMMAND = CMD
RBUILD = build
RINSTALL = install
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
	'Please type make [build | install | check | doc | clean]'

build: build/$(TARGET)

build/$(TARGET): $(RSOURCE)/$(RFILES) $(RDOCUMENTS)/$(RDFILES)
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RBUILD) ..

$(RSOURCE)/$(RFILES):

$(RDOCUMENTS)/$(RDFILES):

install: build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) $(TARGET)

check: build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RCHECK) $(TARGET)

doc: install
	rm -rf $(RDOCUMENTS)/$(RDFILES)
	cd $(RSOURCE); $(REXEC) --vanilla < ../support/document.R
	cd $(RDATA); $(REXEC) --vanilla < ../support/document.data.R
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RBUILD) ..
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) $(TARGET)

clean:
	rm -rf $(RBUILD)/*
