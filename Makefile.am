if USE_PHAML
 PHAML_DIR = libphaml
endif
SUBDIRS = m4macros doc scripts libglimmer libglimmer-solve libglide libglint $(libphaml) utils tests example-drivers

EXTRA_DIST = make.rules glimmer-cism-config.in
bin_SCRIPTS=glimmer-cism-config

datadir = tests example-drivers

data::
	list='$(datadir)'; for dir in $$list; do \
          (cd $$dir && $(MAKE) $@); \
        done
