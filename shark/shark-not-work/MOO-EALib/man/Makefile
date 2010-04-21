SRCS  = 
TRGTS = $(SRCS:%.tex=%.ps)
DVIPS = dvips

.PHONY:	.tex .dvi .ps

all: 	$(TRGTS) 
	@printf "\n  # There are already precompiled manuals for library MOO-EALib"
	@printf "\n  # but there are no sources available.\n"	

%.dvi: %.tex
	latex $< > /dev/null
	./rlatex $< > /dev/null

%.ps: %.dvi
	$(DVIPS) -q -f $< -o $@ 

clean:
	-rm -f *.dvi *.aux *.toc *.log *.bbl *.blg *.glo *.idx *.ilg *.ind

distclean: clean
	-rm -f *.ps






