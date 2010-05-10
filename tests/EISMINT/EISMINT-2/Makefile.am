#if BUILD_GRAPHICS
#  graphics = e2.a.ps e2.b.ps e2.c.ps e2.d.ps e2.f.ps e2.g.ps e2.h.ps e2.i.ps e2.j.ps e2.k.ps e2.l.ps e2.a-comp.ps e2.b-comp.ps e2.c-comp.ps e2.d-comp.ps e2.g-comp.ps e2.h-comp.ps
#endif
data = e2.a.nc e2.b.nc e2.c.nc e2.d.nc e2.f.nc e2.g.nc e2.h.nc e2.i.nc e2.j.nc e2.k.nc e2.l.nc
BUILT_SOURCES = 

EXTRA_DIST = e2.a.config  e2.b.config  e2.c.config  e2.d.config  e2.f.config  e2.g.config  e2.h.config  e2.i.config  e2.j.config  e2.k.config  e2.l.config mound.nc trough.nc

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(data) $(graphics) $(BUILT_SOURCES)

MODEL_BINARY=$(top_builddir)/tests/EISMINT/src/simple_glide

e2.b.nc:        e2.a.nc
e2.c.nc:        e2.a.nc
e2.d.nc:        e2.a.nc
e2.j.nc:        e2.i.nc
e2.l.nc:        e2.k.nc

e2.a-comp.ps:	e2.a.nc
	@GLIMMER_PREFIX@/bin/plotEISMINT2stats.py -e A $^ $@
e2.b-comp.ps:	e2.b.nc e2.a.nc
	@GLIMMER_PREFIX@/bin/plotEISMINT2stats.py -e B $^ $@
e2.c-comp.ps:	e2.c.nc e2.a.nc
	@GLIMMER_PREFIX@/bin/plotEISMINT2stats.py -e C $^ $@
e2.d-comp.ps:	e2.d.nc e2.a.nc
	@GLIMMER_PREFIX@/bin/plotEISMINT2stats.py -e D $^ $@
e2.g-comp.ps:	e2.g.nc
	@GLIMMER_PREFIX@/bin/plotEISMINT2stats.py -e G $^ $@
e2.h-comp.ps:	e2.h.nc
	@GLIMMER_PREFIX@/bin/plotEISMINT2stats.py -e H $^ $@

-include $(top_builddir)/make.rules

.nc.ps:
	@GLIMMER_PREFIX@/bin/plotEISMINT.py -T -1  -p ../PyCF/eismint.prof $< $@
