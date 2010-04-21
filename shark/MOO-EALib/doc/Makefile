PACKAGE = MOO-EALib
SHELL 	= /bin/sh
CXX 	= g++
DVIPS 	= dvips

all: 	doxydoc \
	prepdoc \
	prepare 

#\
#	refman \
#	rename

clean:
	-rm -fR *.o latex/*.dvi latex/*.aux latex/*.toc latex/*.log *.aux *.idx *.log *.toc *.SUP *.dsw *.ncb *.pdb *.obj *.opt *.plg *~ win_tmp

distclean: clean
	-rm -fR prepdoc prepdoc.exe $(PACKAGE)_Reference.html html latex

doxydoc: $(PACKAGE).cfg
	doxygen $(PACKAGE).cfg 

prepdoc.o: prepdoc.cpp
	@echo " ## Compiling documentation preparing program... "
	$(CXX) -c prepdoc.cpp 

prepdoc: prepdoc.o
	@echo " ## Linking documentation preparing program... "
	$(CXX) -o prepdoc prepdoc.o 

prepare: prepdoc
	@echo " ## Preparing documentation... "
	@./prepdoc 

refman: latex/Makefile
	@echo " ## Creating documentation (postscript-version)... "
	cd latex ; \
	make ; \
	make ps 

rename: latex/refman.ps
	@mv latex/refman.ps latex/$(PACKAGE)_Reference.ps 









