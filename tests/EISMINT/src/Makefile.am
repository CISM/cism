bin_PROGRAMS  =  simple_glide
lib_LTLIBRARIES = libglimmer-simple.la

libglimmer_simple_la_SOURCES = simple_forcing.F90
libglimmer_simple_la_LDFLAGS=-all-static

simple_glide_SOURCES = simple_glide.F90
simple_glide_LDADD = libglimmer-simple.la $(top_builddir)/libglide/libglide.la $(top_builddir)/libglimmer-solve/libglimmer-solve.la \
                    $(top_builddir)/libglimmer/libglimmer-IO.la $(top_builddir)/libglimmer/libglimmer.la \
                    $(NETCDF_LDFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LIBS) $(BLAS_LIBS)

AM_FCFLAGS = $(NETCDF_FCFLAGS) $(F90_MODFLAG)$(top_builddir)/libglimmer-solve $(F90_MODFLAG)$(top_builddir)/libglide $(F90_MODFLAG)$(top_builddir)/libglimmer


CLEANFILES = f90_dep.mak *.mod
-include f90_dep.mak
-include $(top_builddir)/make.rules
