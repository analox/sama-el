REM To execute this batch file, you need the free MiKTeX package,
REM available at http://www.miktex.org/latest.html.

latex Reference_EALib.tex
bibtex Reference_EALib
latex Reference_EALib.tex
makeindex Reference_EALib
latex Reference_EALib.tex
dvips Reference_EALib
latex QuickReference_EALib.tex
bibtex QuickReference_EALib
latex QuickReference_EALib.tex
makeindex QuickReference_EALib
latex QuickReference_EALib.tex
dvips QuickReference_EALib
latex Schnellreferenz_EALib.tex
bibtex Schnellreferenz_EALib
latex Schnellreferenz_EALib.tex
makeindex Schnellreferenz_EALib
latex Schnellreferenz_EALib.tex
dvips Schnellreferenz_EALib
latex Referenz_EALib.tex
bibtex Referenz_EALib
latex Referenz_EALib.tex
makeindex Referenz_EALib
latex Referenz_EALib.tex
dvips Referenz_EALib
latex EALib_Structure_Introduction.tex
latex EALib_Structure_Introduction.tex
dvips EALib_Structure_Introduction
latex EALib.tex
bibtex EALib
latex EALib.tex
latex EALib.tex
dvips EALib
