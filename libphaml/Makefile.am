FC = $(MPIFC)

if USE_PHAML_GRAPHICS
   PHAML_GRAPHICS = phaml_graphics
endif

bin_PROGRAMS = phaml_slave simple_phaml $(PHAML_GRAPHICS)


lib_LTLIBRARIES = libphaml_user_mod.la libphaml_support.la  \
     libphaml_example.la libphaml_pde.la

nodist_include_HEADERS := $(wildcard *.mod)
EXTRA_DIST = f90_dep.mak 
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak $(wildcard *.mod)

#the usermod and support libraries are needed for all modules
libphaml_user_mod_la_SOURCES =  phaml_user_mod.F90 
libphaml_user_mod_la_LIBADD = $(ac_cv_phaml_prefix)/lib/libphaml.a

libphaml_support_la_SOURCES = phaml_support.F90 
libphaml_support_la_LIBADD = libphaml_user_mod.la

#---------------------------------------------------------------------
#new phaml modules to be used by glimmer
libphaml_example_la_SOURCES = phaml_example.F90 
libphaml_example_la_LIBADD = libphaml_user_mod.la libphaml_support.la 

#---------------------------------------------------------------------

#the pde library
libphaml_pde_la_SOURCES = phaml_pde.F90 phaml_example_pde.F90
libphaml_pde_la_LIBADD = libphaml_user_mod.la \
    $(top_builddir)/libglide/libglide.la \
    $(top_builddir)/libglimmer-solve/libglimmer-solve.la \
    $(top_builddir)/libglimmer/libglimmer-IO.la \
    $(top_builddir)/libglimmer/libglimmer.la

#this slave program must be compiled for phaml to work
nodist_phaml_slave_SOURCES = $(ac_cv_phaml_prefix)/src/slave.f90 
phaml_slave_LDADD = $(ac_cv_phaml_prefix)/lib/libphaml.a \
    $(MPILIBS) libphaml_pde.la
     

#this builds the graphics slave
nodist_phaml_graphics_SOURCES = $(ac_cv_phaml_prefix)/src/graphmain.f90
phaml_graphics_LDADD = $(ac_cv_phaml_prefix)/lib/libphaml.a	\
        $(MPILIBS) libphaml_pde.la

##############################################################
#a simple driver program
simple_phaml_SOURCES = simple_phaml.F90
simple_phaml_LDADD =  $(ac_cv_phaml_prefix)/lib/libphaml.a \
    $(top_builddir)/libglide/libglide.la \
    $(top_builddir)/libglimmer-solve/libglimmer-solve.la \
    $(top_builddir)/libglimmer/libglimmer-IO.la \
    $(top_builddir)/libglimmer/libglimmer.la \
    $(NETCDF_LDFLAGS) $(NETCDF_LIBS) $(MPILIBS) \
    libphaml_user_mod.la libphaml_example.la libphaml_pde.la


AM_FCFLAGS = $(NETCDF_FCFLAGS) $(F90_MODFLAG)$(top_builddir)/libglimmer-solve \
    $(F90_MODFLAG)$(top_builddir)/libglide $(F90_MODFLAG)$(top_builddir)/libglimmer


f90_dep.mak:	$(srcdir)/*.F90 $(BUILT_SOURCES)
		$(PYTHON) $(top_srcdir)/scripts/f90_dependencies.py -l -o $@ $^

-include f90_dep.mak

