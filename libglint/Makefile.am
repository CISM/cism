lib_LTLIBRARIES = libglint.la

variable_DATA = glint_vars.def glint_mbal_vars.def
variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak 
BUILT_SOURCES = glint_io.F90 glint_mbal_io.F90
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak $(wildcard *.mod)

libglint_la_SOURCES = glimmer_daily_pdd.F90 glimmer_pdd.F90 \
                      glint_climate.F90 glint_initialise.F90 glint_io.F90 \
                      glint_timestep.F90 glint_type.F90 glint_main.F90 \
                      glint_interp.F90 glint_global_grid.F90 \
                      glint_precip_param.F90 glint_mbal.F90 glint_constants.F90\
                      glint_mbal_coupling.F90 glint_mbal_io.F90 \
                      glint_smb.F90 glint_mpinterp.F90 glint_commandline.F90 \
                      glint_global_interp.F90

AM_FCFLAGS = $(F90_MODFLAG)../libglimmer $(F90_MODFLAG)../libglide $(F90_MODFLAG)../libglimmer-solve $(NETCDF_FCFLAGS)

#some special rules for automatically generating ncdf src code
glint_io.F90:	glint_vars.def ../libglimmer/ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/scripts/generate_ncvars.py $^
glint_mbal_io.F90:	glint_mbal_vars.def ../libglimmer/ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/scripts/generate_ncvars.py $^


f90_dep.mak:	$(srcdir)/*.F90 $(BUILT_SOURCES)
		$(PYTHON) $(top_srcdir)/scripts/f90_dependencies.py -l -o $@ $^

-include f90_dep.mak

glide_mask.inc : $(top_srcdir)/libglide/glide_mask.inc