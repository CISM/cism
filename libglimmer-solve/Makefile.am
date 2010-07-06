if EXT_BLAS
 BLAS_DIR =
 BLAS_LIB = 
else
 BLAS_DIR = blas
 BLAS_LIB = blas/libblas.la
endif

SUBDIRS = SLAP $(BLAS_DIR)

lib_LTLIBRARIES = libglimmer-solve.la

libglimmer_solve_la_SOURCES = glimmer_sparse.F90
libglimmer_solve_la_LIBADD = SLAP/libslap.la $(BLAS_LIB)
libglimmer_solve_la_LDFLAGS=-all-static
nodist_include_HEADERS := $(wildcard *.mod)

AM_FCFLAGS = $(F90_MODFLAG)../libglimmer
CLEANFILES = glimmer_vers.F90 f90_dep.mak *.mod

-include f90_dep.mak
-include $(top_builddir)/make.rules
