SUBDIRS = m4macros scripts libglimmer libglimmer-solve libglide libglint \
          utils extensions tests example-drivers doc
if BUILD_PYTHON
  SUBDIRS += PyGC
endif

EXTRA_DIST = glimmer-cism-config.in
bin_SCRIPTS=glimmer-cism-config

datadir = tests example-drivers

data::
	list='$(datadir)'; for dir in $$list; do \
          (cd $$dir && $(MAKE) $@); \
        done

check-data::
	list='$(datadir)'; for dir in $$list; do \
          (cd $$dir && $(MAKE) $@); \
        done

dist-gold-standard::
	@find $(top_builddir)/tests/*/*/*.nc -size +100 -print > test_filenames  
	@tar -cvzf gold-standard_$(distdir).tgz -T test_filenames   
	@rm test_filenames 
	@echo "Output of tests compressed to" gold-standard_$(distdir).tgz