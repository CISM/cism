if DO_DOC_PDF
 pdf_doc = glimmer.pdf
 dirs = dg ug num tut common 
endif
if DO_DOC_PS
 ps_doc = glimmer.ps
 dirs = dg ug num tut common 
endif
doc_DATA = $(ps_doc) $(pdf_doc)
docdir = $(prefix)/doc
SUBDIRS = $(dirs)

EXTRA_DIST = dg.tex dg_pdf.tex dg_ps.tex glimmer.tex glimmer_pdf.tex \
             glimmer_ps.tex num.tex num_pdf.tex num_ps.tex ug.tex ug_pdf.tex \
             ug_ps.tex tut_pdf.tex tut_ps.tex tut.tex glimmer.bib

CLEANFILES = *~ *.aux *.dvi *.log *.out *.toc *.bbl *blg *.html \
             glimmer_html.4ct  glimmer_html.4tc  glimmer_html.css  glimmer_html.idv  \
             glimmer_html.lg  glimmer_html.tmp  glimmer_html.xref

editdoxy = sed -e 's|@top_srcdir[@]|$(top_srcdir)|g' \
               -e 's|@GLIMMER_VERSION[@]|$(GLIMMER_VERSION)|g' \
               -e 's|@PACKAGE_NAME[@]|$(PACKAGE_NAME)|g'
include $(top_srcdir)/m4macros/amdoxygen.am

glimmer-cism.cfg:	$(srcdir)/glimmer-cism.cfg.in Makefile
			$(editdoxy) '$(srcdir)/$@.in' >$@

www::	glimmer.tex
	install -d www
	htlatex glimmer_html "html,3,info" "" " -d./www/ -m 644 "
	$(BIBTEX) glimmer_html
	htlatex glimmer_html "html,3,info" "" " -d./www/ -m 644 "
	find $(SUBDIRS) -name '*.png' -exec install -m 644 -D {} ./www/{} \;

.tex.dvi::
	$(LATEX) $<
	$(BIBTEX) $*
	$(LATEX) $<
	$(LATEX) $<

%.ps:	%_ps.dvi
	dvips -o $@ $<

%.pdf:	%_pdf.dvi
	$(DVIPDFT) -o $@ $<
