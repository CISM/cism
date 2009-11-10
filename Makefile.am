SUBDIRS = m4macros doc scripts libglimmer libglimmer-solve libglide libglint utils tests example-drivers
EXTRA_DIST = make.rules glimmer-cism-config.in

datadir = tests example-drivers

data::
	list='$(datadir)'; for dir in $$list; do \
          (cd $$dir && $(MAKE) $@); \
        done
