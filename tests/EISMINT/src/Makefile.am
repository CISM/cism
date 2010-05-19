bin_PROGRAMS  =  simple_glide
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak $(wildcard *.mod)
simple_glide_SOURCES = simple_forcing.F90 simple_glide.F90 
simple_glide_LDADD = $(top_builddir)/libglide/libglide.la $(top_builddir)/libglimmer-solve/libglimmer-solve.la \
                    $(top_builddir)/libglimmer/libglimmer-IO.la $(top_builddir)/libglimmer/libglimmer.la \
                    $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS) $(BLAS_LIBS)

AM_FCFLAGS = $(NETCDF_FCFLAGS) $(F90_MODFLAG)$(top_builddir)/libglimmer-solve $(F90_MODFLAG)$(top_builddir)/libglide $(F90_MODFLAG)$(top_builddir)/libglimmer

f90_dep.mak:	$(simple_glide_SOURCES)
		$(PYTHON) $(top_srcdir)/scripts/f90_dependencies.py -o $@ $^

-include f90_dep.mak
