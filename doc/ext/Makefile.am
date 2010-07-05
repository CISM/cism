SUBDIRS=
if BUILD_EXT_EROSION
SUBDIRS += erosion
endif

BUILT_SOURCES = ext.tex

editext = sed -e 's|@build_extension_erosion[@]|$(build_extension_erosion)|g' 

ext.tex:	ext.tex.in Makefile
		$(editext) '$(srcdir)/$@.in' >$@
