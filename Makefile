RDFILES = *.Rd
REXEC = R
RCOMMAND = CMD
RBUILD = build
RINSTALL = install
RCHECK = check

nothing:
	@echo \
	'Please type "make build", "make install", "make check", or "make clean"'

build: build/*

build/*: R/*
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RBUILD) ..

R/*:

install: build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) OpenMx*

check: build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RCHECK) OpenMx*

clean:
	rm -rf $(RBUILD)/OpenMx*
