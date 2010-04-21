doxygen Rng.cfg
CALL vcvars32
cl prepdoc.cpp
prepdoc
rem cd latex
rem latex refman.tex
rem makeindex refman.idx
rem latex refman.tex
rem latex refman.tex
rem dvips refman
rem ren refman.ps Rng_Reference.ps
rem cd ..


