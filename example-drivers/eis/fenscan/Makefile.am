#if BUILD_GRAPHICS
#  graphics = 
#endif
data = fenscan.nc fenscan-gthf.nc
BUILT_SOURCES = 

EXTRA_DIST = fenscan.config fenscan-gthf.config europe.20km.nc mb2.data temp-exp.model specmap.data

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(data) $(graphics) $(BUILT_SOURCES)

MODEL_BINARY=$(top_builddir)/example-drivers/eis/src/eis_glide

-include $(top_builddir)/make.rules

