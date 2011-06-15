SUBDIRS = EISMINT exact-isothermal

data::
		list='$(SUBDIRS)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done
