tests = EISMINT-3
SUBDIRS = src EISMINT-3

data::
		list='$(tests)'; for dir in $$list; do \
		(cd $$dir && $(MAKE) $@); \
		done

check-data::
		list='$(tests)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done

dist-gold-standard::
		list='$(tests)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done

