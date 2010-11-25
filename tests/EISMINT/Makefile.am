tests = benchmarks EISMINT-1 EISMINT-2 isos lin-non_lin
SUBDIRS = src $(tests)
if BUILD_EXT_EROSION
 SUBDIRS += erosion
endif

graphics::
		list='$(tests)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done

check-data::
		list='$(tests)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done
