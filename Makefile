RDFILES = *.Rd
REXEC = R
RSOURCE = R
RCOMMAND = CMD
RBUILD = build
RINSTALL = install
RCHECK = check
TARGET = OpenMx_0.1-0.tar.gz

nothing:
	@echo \
	'Please type make [build | install | check | doc | clean]'

build: build/$(TARGET)

build/$(TARGET): R/*.R man/*.Rd
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RBUILD) ..

R/*.R:

man/*.Rd:

install: build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) $(TARGET)

check: build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RCHECK) $(TARGET)

doc: install
	cd $(RSOURCE); $(REXEC) --vanilla < ../build/document.R
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) $(TARGET)

clean:
	rm -rf $(RBUILD)/$(TARGET)
