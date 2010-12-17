SUBDIRS = eis EISMINT-3 glint
tests = eis EISMINT-3

data::
		list='$(tests)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done

dist-gold-standard::
		list='$(tests)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done

