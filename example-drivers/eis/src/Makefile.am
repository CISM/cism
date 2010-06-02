lib_LTLIBRARIES = libglimmer-eis.la
bin_PROGRAMS  =  eis_glide

variable_DATA = eis_vars.def
variabledir = $(pkgdatadir)
dist_variable_DATA = $(variable_DATA)

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak
BUILT_SOURCES = eis_io.F90
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak $(wildcard *.mod)

libglimmer_eis_la_SOURCES = eis_forcing.F90 eis_ela.F90 eis_temp.F90 eis_slc.F90 eis_cony.F90 \
                       eis_io.F90 eis_types.F90

eis_glide_SOURCES = eis_glide.F90 
eis_glide_LDADD = libglimmer-eis.la $(top_builddir)/libglide/libglide.la $(top_builddir)/libglimmer-solve/libglimmer-solve.la \
                    $(top_builddir)/libglimmer/libglimmer-IO.la $(top_builddir)/libglimmer/libglimmer.la \
                    $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS) $(BLAS_LIBS) $(LAPACK_LIBS) $(PARDISO_LIBS) $(UMFPACK_LIBS)

AM_FCFLAGS = $(NETCDF_FCFLAGS) $(F90_MODFLAG)$(top_builddir)/libglimmer-solve $(F90_MODFLAG)$(top_builddir)/libglide $(F90_MODFLAG)$(top_builddir)/libglimmer



#some special rules for automatically generating ncdf src code
eis_io.F90:	eis_vars.def $(top_srcdir)/libglimmer/ncdf_template.F90.in
		$(PYTHON) $(top_srcdir)/scripts/generate_ncvars.py $^

f90_dep.mak:	$(srcdir)/*.F90 $(BUILT_SOURCES)
		$(PYTHON) $(top_srcdir)/scripts/f90_dependencies.py -l -o $@ $^

-include f90_dep.mak

# this fixes the build order: f90_dep.mak has .lo for eis_glide.o
# but since eis_glide isn't a library, there is no eis_glide.lo
# target, there is only a eis_glide.o target which never depends on
# the other targets, so it builds early and is missing it's .mod files
eis_glide.o : eis_glide.lo

glide_mask.inc : $(top_srcdir)/libglide/glide_mask.inc