SRCS  = Reference_EALib.tex \
	QuickReference_EALib.tex \
	Schnellreferenz_EALib.tex \
	Referenz_EALib.tex \
	EALib_Structure_Introduction.tex \
	EALib.tex

TRGTS = $(SRCS:%.tex=%.ps)
DVIPS = dvips

.PHONY:	.tex .dvi .ps

all: $(TRGTS) 

Reference_EALib.dvi: Reference_EALib.tex
	latex Reference_EALib > /dev/null
	-bibtex Reference_EALib > /dev/null
	latex Reference_EALib > /dev/null
	makeindex -q  Reference_EALib > /dev/null
	latex Reference_EALib > /dev/null

EALib.dvi: EALib.tex
	latex $< > /dev/null
	bibtex $(<:%.tex=%) > /dev/null
	./rlatex $< > /dev/null

QuickReference_EALib.dvi: QuickReference_EALib.tex
Schnellreferenz_EALib.dvi: Schnellreferenz_EALib.tex
Referenz_EALib.dvi: Referenz_EALib.tex
EALib_Structure_Introduction.dvi: EALib_Structure_Introduction.tex

%.dvi: %.tex
	latex $< > /dev/null		
	./rlatex $< > /dev/null

%.ps: %.dvi
	$(DVIPS) -q -f $< -o $@ 

clean:
	-rm -f *.dvi *.aux *.toc *.log *.bbl *.blg *.glo *.idx *.ilg *.ind

distclean: clean
	-rm -f *.ps






