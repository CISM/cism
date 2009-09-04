SUBDIRS = SLAP

lib_LTLIBRARIES = libglimmer-solve.la

libglimmer_solve_la_SOURCES = glimmer_sparse.F90
libglimmer_solve_la_LIBADD = SLAP/libslap.la
nodist_include_HEADERS := $(wildcard *.mod)

AM_FCFLAGS = -I../libglimmer
