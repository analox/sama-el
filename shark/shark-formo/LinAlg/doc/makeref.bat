doxygen LinAlg.cfg
CALL vcvars32
cl prepdoc.cpp
prepdoc
cd latex
rem latex refman.tex
rem makeindex refman.idx
rem latex refman.tex
rem latex refman.tex
rem dvips refman
rem ren refman.ps LinAlg_Reference.ps
rem cd ..


