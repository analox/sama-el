REM To execute this batch file, you need the free MiKTeX package,
REM available at http://www.miktex.org/latest.html.

latex Array-mini.tex
latex Array-mini.tex
dvips Array-mini
latex Reference_Array.tex
bibtex Reference_Array
latex Reference_Array.tex
makeindex Reference_Array
latex Reference_Array.tex
dvips Reference_Array
