lib_LTLIBRARIES = libglimmer.la libglimmer-IO.la
bin_PROGRAMS  =  nc2config
check_PROGRAMS = test_config test_integrate test_commandline test_ts test_writestats test_filenames
TESTS = test_filenames

variable_DATA = ncdf_template.F90.in
variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak 
CLEANFILES = glimmer_vers.F90 f90_dep.mak *.mod
libglimmer_la_SOURCES = glimmer_vers.F90 glimmer_commandline.F90 \
                       glimmer_global.F90 glimmer_utils.F90\
                       glimmer_log.F90 glimmer_paramets.F90 glimmer_physcon.F90 kelvin.F90 \
                       profile.F90  glimmer_searchcircle.F90 \
                       glimmer_coordinates.F90 glimmer_deriv.F90 \
                       glimmer_integrate.F90 \
                       glimmer_map_CFproj.F90 glimmer_map_proj4.F90 \
                       glimmer_map_types.F90 glimmer_map_init.F90 \
                       glimmer_map_trans.F90 glimmer_scales.F90
nodist_libglimmer_la_SOURCES = glimmer_vers.F90
libglimmer_la_LDFLAGS=-all-static

libglimmer_IO_la_SOURCES = glimmer_anomcouple.F90 glimmer_filenames.F90 glimmer_ncdf.F90 \
                          glimmer_ts.F90 glimmer_config.F90 \
                          glimmer_writestats.F90 ncdf_utils.F90 writestats.c writestats.h cfortran.h 
libglimmer_IO_la_LDFLAGS=-all-static

nc2config_SOURCES = nc2config.F90

test_config_SOURCES = test_config.F90
test_integrate_SOURCES = test_integrate.F90
test_commandline_SOURCES = test_commandline.F90
test_ts_SOURCES = test_ts.F90
test_writestats_SOURCES = test_writestats.F90
test_filenames_SOURCES = test_filenames.F90

AM_FCFLAGS = $(NETCDF_FCFLAGS)
LDADD = libglimmer-IO.la libglimmer.la $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS)
AM_CFLAGS=$(CFORTRAN_CFLAGS)

# handling subversion versioning
if HAVE_SVN
svnversion.sed::
		$(srcdir)/createSVNVersion.sh $(top_srcdir)
glimmer_vers.F90:	svnversion.sed glimmer_vers.F90.in
		sed -f svnversion.sed glimmer_vers.F90.in > glimmer_vers.F90
else
glimmer_vers.F90:	glimmer_vers.F90.in
		cp $< $@
endif

-include f90_dep.mak
-include $(top_builddir)/make.rules
