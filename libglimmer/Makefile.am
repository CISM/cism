if EXT_BLAS
 BLAS_DIR =
else
 BLAS_DIR = blas
endif

SUBDIRS = $(BLAS_DIR) solve

lib_LIBRARIES = libglimmer.a
bin_PROGRAMS  =  
check_PROGRAMS = test_config test_integrate test_commandline


EXTRA_DIST = f90_dep.mak 
CLEANFILES = $(BUILT_SOURCES) glimmer_vers.F90 f90_dep.mak
libglimmer_a_SOURCES = glimmer_vers.F90 glimmer_commandline.F90 \
                       glimmer_global.F90 glimmer_utils.F90\
                       glimmer_log.F90 glimmer_paramets.F90 glimmer_physcon.F90 kelvin.F90 \
                       profile.F90  \
                       glimmer_coordinates.F90 \
                       glimmer_integrate.F90 \
                       glimmer_map_CFproj.F90 glimmer_map_proj4.F90 \
                       glimmer_map_types.F90 glimmer_map_init.F90 \
                       glimmer_map_trans.F90
nodist_libglimmer_a_SOURCES = glimmer_vers.F90

test_config_SOURCES = test_config.F90 glimmer_config.F90 glimmer_log.F90
test_integrate_SOURCES = test_integrate.F90 glimmer_global.F90 \
                         glimmer_physcon.F90 glimmer_integrate.F90
test_commandline_SOURCES = test_commandline.F90 glimmer_commandline.F90

AM_FCFLAGS = $(NETCDF_CPPFLAGS) $(HDF5_CPPFLAGS)
LDADD = $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS)

# handling subversion versioning
if HAVE_SVN
svnversion.sed::
		$(srcdir)/createSVNVersion.sh
glimmer_vers.F90:	svnversion.sed glimmer_vers.F90.in
		sed -f svnversion.sed glimmer_vers.F90.in > glimmer_vers.F90
else
glimmer_vers.F90:	glimmer_vers.F90.in
		cp $< $@
endif


f90_dep.mak:	$(srcdir)/*.F90 $(BUILT_SOURCES) glimmer_vers.F90
		$(PYTHON) $(top_srcdir)/scripts/f90_dependencies.py -o $@ $^

-include f90_dep.mak
