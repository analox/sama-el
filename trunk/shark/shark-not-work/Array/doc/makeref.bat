CALL vcvars32
cl prepdoc1.cpp
cl prepdoc2.cpp
prepdoc1
doxygen Array.cfg
prepdoc2
rem del ..\..\include\Array\ArrayOp.h
rem copy arraytmpfile ..\..\include\Array\ArrayOp.h
rem del arraytmpfile
rem cd latex
rem latex refman.tex
rem makeindex refman.idx
rem latex refman.tex
rem latex refman.tex
rem dvips refman
rem ren refman.ps Array_Reference.ps
rem cd ..


