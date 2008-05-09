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

build/*:
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RBUILD) ..

install: build
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RINSTALL) OpenMx*

check:
	rm -rf base
	ln -s . base
	$(REXEC) $(RCOMMAND) $(RCHECK) base
	rm -rf base

clean:
	rm -rf base
	rm -rf base.Rcheck
	rm -rf $(RBUILD)/OpenMx*
