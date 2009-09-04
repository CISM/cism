bin_PROGRAMS  =  glint_example glex_ebm

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak eismint3_vars.def
BUILT_SOURCES = 
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak



glint_example_SOURCES = glint_example.F90 glint_example_clim.F90
glex_ebm_SOURCES = glex_ebm.F90 glex_ebm_clim.F90

LDADD = $(top_builddir)/libglint/libglint.la  $(top_builddir)/libglide/libglide.la \
        $(top_builddir)/libglimmer-solve/libglimmer-solve.la \
        $(top_builddir)/libglimmer/libglimmer-IO.la $(top_builddir)/libglimmer/libglimmer.la \
        $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS) $(BLAS_LIBS)

AM_FCFLAGS = $(NETCDF_CPPFLAGS) $(HDF5_CPPFLAGS) -I$(top_builddir)/libglimmer-solve -I$(top_builddir)/libglide \
             -I$(top_builddir)/libglimmer -I$(top_builddir)/libglint



f90_dep.mak:	$(srcdir)/*.F90 $(BUILT_SOURCES)
		$(PYTHON) $(top_srcdir)/scripts/f90_dependencies.py -o $@ $^

-include f90_dep.mak
