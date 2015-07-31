SUBDIRS := pkg/SBML2ADOLC/. pkg/SBML2FORTRAN/. pkg/LIMEX4_3A/. pkg/NLSCON/. ext/LimexDL/Model_ODE/. ext/. $(wildcard ext/*/.)  
# e.g. SUBDIRS == "... ext/. ext/foo/. ext/bar/."
TARGETS := all clean mclean # whatever else, but must not contain '/'

# ext/foo/.all ext/bar/.all ext/foo/.clean ext/bar/.clean
SUBDIRS_TARGETS := \
    $(foreach t,$(TARGETS),$(addsuffix $t,$(SUBDIRS)))

.PHONY : $(TARGETS) $(SUBDIRS_TARGETS)

# static pattern rule, expands into:
# all clean : % : foo/.% bar/.%
$(TARGETS) : % : $(addsuffix %,$(SUBDIRS))
	@echo 'Done "$*" target'

# here, for foo/.all:
#   $(@D) is foo
#   $(@F) is .all, with leading period
#   $(@F:.%=%) is just all
$(SUBDIRS_TARGETS) :
	$(MAKE) -C $(@D) $(@F:.%=%)

