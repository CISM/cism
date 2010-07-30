#if BUILD_GRAPHICS
#  graphics = 
#endif
data = b1.nc  b2.nc  b3.nc  b4.nc  b5.nc  b6.nc 
BUILT_SOURCES = 

EXTRA_DIST = b1.config  b2.config  b3.config  b4.config  b5.config  b6.config

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(data) $(graphics) *.jid $(BUILT_SOURCES)

MODEL_BINARY=$(top_builddir)/tests/EISMINT/src/simple_glide

include $(top_srcdir)/extra_rules.am

