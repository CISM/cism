lib_LTLIBRARIES = libglide.la
bin_PROGRAMS  =  
check_PROGRAMS = test_setup

variable_DATA = glide_vars.def glide_lithot_vars.def time_vars.def
variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak $(RESTART_INC_INS)
BUILT_SOURCES = glide_io.F90  glide_lithot_io.F90
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak

libglide_la_SOURCES = glide.F90 glide_setup.F90 glide_types.F90 glide_temp.F90 \
                      glide_thck.F90 glide_velo.F90 glide_mask.F90 \
                      glide_stop.F90 glide_io.F90 glide_lithot_io.F90 glide_nc_custom.F90 \
                      isostasy.F90 isostasy_el.F90 isostasy_setup.F90 \
                      isostasy_types.F90 glide_lithot.F90 glide_lithot3d.F90 \
                      glide_lithot1d.F90 glide_profile.F90 glide_diagnostics.F90 \
                      glimmer_routing.F90 glimmer_ncio.F90 glimmer_ncparams.F90  \
                      glimmer_restart_gcm.F90

test_setup_SOURCES = test_setup.F90
test_setup_LDADD = libglide.la ../libglimmer/libglimmer.la ../libglimmer/libglimmer-IO.la $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS)

AM_FCFLAGS = $(F90_MODFLAG)../libglimmer $(F90_MODFLAG)../libglimmer-solve $(NETCDF_FCFLAGS)
AM_LDFLAGS=-all-static

#AM_FCFLAGS = $(NETCDF_FCFLAGS)
#LDADD = $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS)
#AM_CFLAGS=$(CFORTRAN_CFLAGS)

#some special rules for automatically generating ncdf src code
%_io.F90:	$(srcdir)/%_vars.def $(top_srcdir)/libglimmer/ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/scripts/generate_ncvars.py $^


f90_dep.mak:	$(srcdir)/*.F90 $(BUILT_SOURCES)
		$(PYTHON) $(top_srcdir)/scripts/f90_dependencies.py -l -o $@ $^

-include f90_dep.mak
