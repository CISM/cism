graphics = e2-comp.eps $(subst .nc,_btemp.eps,$(data)) $(subst .nc,_bvel.eps,$(data))
data = l-5a.nc nl-5a.nc l-10a.nc nl-10a.nc l-20a.nc nl-20a.nc l-30a.nc nl-30a.nc
BUILT_SOURCES = 

EXTRA_DIST = l-10a.config  l-20a.config  l-5a.config  nl-10a.config  nl-20a.config  nl-5a.config l-30a.config  nl-30a.config

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(data) $(graphics) *.jid $(BUILT_SOURCES)

MODEL_BINARY=$(top_builddir)/tests/EISMINT/src/simple_glide

e2-comp.eps:	$(data)
		$(plotEISMINT2stats) -e H $^ -o $@

#e2-comp-t.ps:	$(data)
#	@GLIMMER_PREFIX@/bin/plotCFstats.py -e -m -p ../PyCF/eismint.prof $^ $@


%_btemp.eps:	%.nc
		$(plotGCvar) --pmt -T -1 -cthk -vtemp -l -1 --leg v $^ -o $@

%_bvel.eps:	%.nc
		$(plotGCvar) -T -1 -cthk -vbvel -l -1 --leg v $^ -o $@

include $(top_srcdir)/extra_rules.am
