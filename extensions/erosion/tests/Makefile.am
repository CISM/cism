TEST_CONFIGS = erosion0.config erosion1.config erosion2.config erosion3.config erosion4.config \
               er-transport0.config er-transport1.config er-transport2.config er-transport3.config er-transport4.config fenscan.config

EXTRA_DIST = trans_velo.config.in trans_velo.py $(TEST_CONFIGS)

data = $(subst config,nc,$(TEST_CONFIGS))

MODEL_BINARY=$(top_builddir)/extensions/erosion/src/simple_erosion

fenscan.nc:	fenscan.config
		$(PYTHON) $(top_srcdir)/scripts/glide_launch.py $(LAUNCH_OPTIONS) -m $(top_builddir)/extensions/erosion/src/eis_erosion @SGE_OPTS@ $<


-include $(top_builddir)/make.rules
