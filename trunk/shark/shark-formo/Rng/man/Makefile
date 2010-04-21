SRCS  = Reference_Rng.tex

TRGTS = $(SRCS:%.tex=%.ps)
DVIPS = dvips

.PHONY:	.tex .dvi .ps

all: $(TRGTS) 

Reference_Rng.dvi: Reference_Rng.tex
	latex Reference_Rng > /dev/null
	-bibtex Reference_Rng > /dev/null
	latex Reference_Rng > /dev/null
	makeindex -q Reference_Rng > /dev/null
	latex Reference_Rng > /dev/null

#Reference_Rng.dvi: Reference_Rng.tex

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




