## This is a new, repo Makefile for viraltraj 2022 May 29 (Sun)

current: target
-include target.mk

# -include makestuff/perl.def

# vim_session:
#	 bash -cl "vmt"

######################################################################

autopipeR = defined
pardirs += R

######################################################################

Sources += $(wildcard *.R)

baseparam.Rout: baseparam.R
	$(wrapR)
	
simulation_growth_sensitivity_k.Rout: simulation_growth_sensitivity_k.R baseparam.rda R/viraltraj.rda R/simgrowth.rda
	$(pipeR)

summarize_growth_sensitivity_k.Rout: summarize_growth_sensitivity_k.R simulation_growth_sensitivity_k.rds baseparam.rda R/viraltraj.rda R/util.rda
	$(pipeR)
	
simulation_growth_sensitivity_corr.Rout: simulation_growth_sensitivity_corr.R baseparam.rda R/viraltraj.rda R/simgrowth.rda
	$(pipeR)

summarize_growth_sensitivity_corr.Rout: summarize_growth_sensitivity_corr.R simulation_growth_sensitivity_corr.rda baseparam.rda R/viraltraj.rda R/util.rda
	$(pipeR)

######################################################################
### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

Makefile: makestuff/00.stamp
makestuff/%.stamp:
	- $(RM) makestuff/*.stamp
	(cd makestuff && $(MAKE) pull) || git clone $(msrepo)/makestuff
	touch $@

-include makestuff/os.mk

-include makestuff/texi.mk
-include makestuff/pipeR.mk
-include makestuff/pdfpages.mk

-include makestuff/git.mk
-include makestuff/visual.mk
