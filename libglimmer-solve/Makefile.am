SUBDIRS = SLAP

lib_LTLIBRARIES = libglimmer-solve.la

libglimmer_solve_la_SOURCES = glimmer_sparse.F90
libglimmer_solve_la_LIBADD = SLAP/libslap.la
nodist_include_HEADERS := $(wildcard *.mod)

AM_FCFLAGS = -I../libglimmer

f90_dep.mak:	$(srcdir)/*.F90 $(BUILT_SOURCES) glimmer_vers.F90
		$(PYTHON) $(top_srcdir)/scripts/f90_dependencies.py -l -o $@ $^

-include f90_dep.mak
