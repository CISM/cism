EXTRA_DIST = eis_varlist.tex glide_varlist.tex glint_varlist.tex intro.tex \
             varlist.tex drivers.tex  example.tex glide.tex installation.tex \
             ug.tex varlist.tex.in glint.tex glint_api.tex mbals.tex \
             glum_api.tex glide_api.tex

#BUILT_SOURCES = glide_varlist.tex eis_varlist.tex glint_varlist.tex
CLEANFILES = $(BUILT_SOURCES)

SUBDIRS = figs

glide_varlist.tex:	@GLIMMER_PREFIX@/share/glimmer/glide_vars.def varlist.tex.in
			$(PYTHON) @GLIMMER_PREFIX@/bin/generate_ncvars.py $^
eis_varlist.tex:	@GLIMMER_PREFIX@/share/glimmer/eis_vars.def varlist.tex.in
			$(PYTHON) @GLIMMER_PREFIX@/bin/generate_ncvars.py $^
glint_varlist.tex:	@GLIMMER_PREFIX@/share/glimmer/glint_vars.def varlist.tex.in
			$(PYTHON) @GLIMMER_PREFIX@/bin/generate_ncvars.py $^
