TARGETS := $(wildcard */) 
MTARGETS := $(foreach t,$(TARGETS),$(addsuffix Makefile,$t))

all : prepare $(MTARGETS)

prepare : 
	/bin/rm -f $(MTARGETS)

$(MTARGETS) : $(@F)
	cd $(@D) && /usr/bin/ruby extconf.rb

clean :
	# /bin/rm -f $(MTARGETS)

mclean :
	/bin/rm -f $(MTARGETS)
	@for a in $(MTARGETS); do cp ./mkfilestub $$a; done
