SUBDIRS = EISMINT exact-isothermal

data::
		list='$(SUBDIRS)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done

check-data::
		list='$(SUBDIRS)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done

dist-gold-standard::
		list='$(SUBDIRS)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done