SRCS  = Array-mini.tex \
	Reference_Array.tex

TRGTS = $(SRCS:%.tex=%.ps)
DVIPS = dvips

.PHONY:	.tex .dvi .ps

all: $(TRGTS) 

Reference_Array.dvi: Reference_Array.tex
	latex Reference_Array > /dev/null
	-bibtex Reference_Array > /dev/null
	latex Reference_Array > /dev/null
	makeindex -q Reference_Array > /dev/null
	latex Reference_Array > /dev/null

Array-mini.dvi: Array-mini.tex

%.dvi: %.tex
	./rlatex $< > /dev/null

%.ps: %.dvi
	$(DVIPS) -q -f $< -o $@ 

clean:
	-rm -f *.dvi *.aux *.toc *.log *.bbl *.blg *.glo *.idx *.ilg *.ind

distclean: clean
	-rm -f *.ps

dos2unix:
	@printf "### Removing carriage returns from all source files..."
	@for i in `find . -type f \( -name "*.tex" -o -name "*.txt" -o -name "*.sty" -o -name "*.cls" -o -name "rlatex" -o -name "Makefile.in" \) -print`; do mv $$i tmp ; tr -d "\r" < tmp > $$i ; rm -f tmp ; done
	@echo "done."



