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
nodist_include_HEADERS := $(wildcard *.mod)

AM_FCFLAGS = -I../libglimmer

f90_dep.mak:	$(srcdir)/*.F90 $(BUILT_SOURCES) glimmer_vers.F90
		$(PYTHON) $(top_srcdir)/scripts/f90_dependencies.py -l -o $@ $^

-include f90_dep.mak
