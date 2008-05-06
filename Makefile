RDFILES = *.Rd
REXEC = R
RCOMMAND = CMD
RBUILD = build
RTEST = check

nothing:
	@echo \
	'Please type "make build", "make test", or "make clean"'


build: clean
	cd $(RBUILD); $(REXEC) $(RCOMMAND) $(RBUILD) ..

test: clean
	ln -s . base
	$(REXEC) $(RCOMMAND) $(RTEST) base
	rm base

clean:
	rm -rf base.Rcheck
	rm -rf $(RBUILD)/OpenMx*
