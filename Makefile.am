SUBDIRS = m4macros doc scripts libglimmer libglimmer-solve libglide libglint \
          utils tests example-drivers extensions
if BUILD_PYTHON
  SUBDIRS += PyGC
endif
if USE_PHAML
 SUBDIRS += libphaml
endif

EXTRA_DIST = make.rules.in glimmer-cism-config.in
bin_SCRIPTS=glimmer-cism-config
BUILT_SOURCES = make.rules

datadir = tests example-drivers

data::
	list='$(datadir)'; for dir in $$list; do \
          (cd $$dir && $(MAKE) $@); \
        done

make.rules:	Makefile make.rules.in
		sed -e 's|@AM_SGE_OPTS[@]|$(AM_SGE_OPTS)|g' '$(srcdir)/$@.in' >$@
