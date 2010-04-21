REM To execute this batch file, you need the free MiKTeX package,
REM available at http://www.miktex.org/latest.html.

latex Reference_Rng.tex
bibtex Reference_Rng
latex Reference_Rng.tex
makeindex Reference_Rng
latex Reference_Rng.tex
dvips Reference_Rng
