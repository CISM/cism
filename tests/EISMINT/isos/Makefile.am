graphics = profile.eps #spots.eps 
data = isos.llfa.nc isos.llra.nc isos.elfa.nc isos.elra.nc
BUILT_SOURCES = 

EXTRA_DIST = isos.elfa.config  isos.elra.config  isos.llfa.config  isos.llra.config

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(data) $(graphics) *.jid $(BUILT_SOURCES)

MODEL_BINARY=$(top_builddir)/tests/EISMINT/src/simple_glide

#spots.eps:	isos.llfa.nc isos.llra.nc isos.elfa.nc isos.elra.nc
#		$(plotGCSpot) --ij 15 15  -vtopg $^ -o $@
#		$(plotGCSpot) --ij 15 15 --ij 15 22 --ij 15 27 -vtopg $^ -o $@

profile.eps:	isos.llfa.nc isos.llra.nc isos.elfa.nc isos.elra.nc
		$(plotGCProfile) -t 18  -v topg -p $(top_srcdir)/PyGC/data/eismint.prof $^ -o $@

include $(top_srcdir)/extra_rules.am
