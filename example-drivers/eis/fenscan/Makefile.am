graphics = europe.20km.eps fenscan_is-17.eps fenscan-gthf_is-17.eps fenscan_bvel-17.eps fenscan-gthf_bvel-17.eps profile-20.eps profile-comp.eps

data = fenscan.nc fenscan-gthf.nc
BUILT_SOURCES = 

EXTRA_DIST = fenscan.config fenscan-gthf.config europe.20km.nc mb2.data temp-exp.model specmap.data prof_a.data

MOSTLYCLEANFILES = $(graphics)
CLEANFILES = $(data) $(graphics) $(BUILT_SOURCES)

MODEL_BINARY=$(top_builddir)/example-drivers/eis/src/eis_glide

-include $(top_builddir)/make.rules


europe.20km.eps:	europe.20km.nc
			$(plotGCvar) -v topg --not_projected -p prof_a.data --legend=v $^ -o $@
%_is-17.eps:		%.nc
			$(plotGCvar) -v is -t -17 -c thk --land -i is --legend=v $^ -o $@
%_bvel-17.eps:		%.nc
			$(plotGCvar) -v bvel -t -17 -c thk --land --vectors bvel --legend=v $^ -o $@

profile-20.eps:		fenscan.hot.nc
			$(plotGCProfile) -p prof_a.data --not_projected --legend -v temp -v bmlt  -t -20 $^ -o $@

profile-comp.eps:	fenscan.hot.nc fenscan-gthf.hot.nc
			$(plotGCProfile) -p prof_a.data --not_projected -v is -v bheatflx -v bmlt -t -20 --size 3 3.5 $^ -o $@
