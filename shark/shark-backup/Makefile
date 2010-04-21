############### get machine type
ARCH_MACH	= i686


# top level makefile for Shark
# runpath for Linux: -Wl,-rpath,${LDRUNPATH}

############### package name and version number
#
PACKAGE		= Shark
VERSION		= 1.4.4


############### variables set by configure script
#

FIND		= /usr/bin/find
MAKE		= /usr/bin/make
PERL		= /usr/bin/perl
ARCH_MACH       = i686
HOST		= Linux
CFLAGS		= -D__LINUX__ -fpic
DEFS		= -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" 
LD_SHARED	= -shared
X_LIBS_R	= 
X_LIBS		=  -L/usr/X11R6/lib
LIBEXT		= .so

PWD		= $(shell pwd)

ifeq ($(INSTALLDIR),)
 INSTALLDIR = $(HOME)/$(PACKAGE)
endif

LDLIBS 			= 
FLAG_PAGODE		= 
FLAG_MOEALIB		= __is_existent__
FLAG_SPREADCAT		= 
FLAG_SIMANN		= 
FLAG_SELFTEST		= 
SELFTEST_DIR		= notestsdir

ifdef LDLIBS
	USE_MOTIF	=	yes
else
	USE_MOTIF	=	no
endif

############### list of libraries to work on
#

ifeq ($(USE_MOTIF), yes)
	LIBRARIES = 		Array \
				Rng \
				Metric \
				LinAlg \
				Mixture \
				EALib \
				TestData \
				ReClaM 
	ifdef FLAG_SPREADCAT
		LIBRARIES +=	SpreadCAT
	endif
	ifdef FLAG_PAGODE
		LIBRARIES +=	Pagode
	endif
else
	LIBRARIES = 		Array \
				Rng \
				Metric \
				LinAlg \
				Mixture \
				EALib \
				TestData \
				ReClaM 
endif

ifdef FLAG_SIMANN
	LIBRARIES +=	SimAnn
endif

ifdef FLAG_MOEALIB
		LIBRARIES +=	MOO-EALib
endif

LIBRARIES_D:= $(LIBRARIES:%=%-debug)

REFERENCES =			Array 	\
				Rng   	\
				LinAlg	\
				ReClaM \
				EALib \
				Metric \
				Mixture
ifdef FLAG_SIMANN
	        REFERENCES +=	SimAnn
endif

ifdef FLAG_MOEALIB
	        REFERENCES +=	MOO-EALib
endif





############### general make commands
#

.PHONY:
		everything all debug examples examples-debug man doc \
		clean tar dist distclean depend independ install     \
		install-ini install-clean rights configure dos2unix  \
		dos2unix_all unix2dos admin_diff

all:		$(LIBRARIES)

debug:		$(LIBRARIES_D)

everything:	$(LIBRARIES) $(LIBRARIES_D) examples examples-debug man doc

configure: \
	Array-conf	\
	EALib-conf	\
	LinAlg-conf	\
	Metric-conf	\
	Mixture-conf	\
	MOO-EALib-conf	\
	Pagode-conf	\
	ReClaM-conf	\
	Rng-conf	\
	SimAnn-conf	\
	SpreadCAT-conf	\
	TestData-conf

Array::\
	Array-conf
	@echo
	@echo "### Building Array library"
	@( cd Array; $(MAKE) -s )
	@echo "### Array library built."

Array-debug::\
	Array-conf
	@echo
	@echo "### Building Array library (debug)"
	@( cd Array; $(MAKE) -s debug )
	@echo "### Array library (debug) built."

Array-conf::\
	Array/Makefile

Array/Makefile::\
	Array/Makefile.in Array/configure.in
	@printf "\n### Configuring Array... "
	@( cd Array; autoconf; ./configure 1>/dev/null 2>/dev/null )
	@echo "Done."

EALib::\
	EALib-conf
	@echo
	@echo "### Building EALib library"
	@( cd EALib; $(MAKE) -s )
	@echo "### EALib library built."

EALib-debug::\
	EALib-conf
	@echo
	@echo "### Building EALib library (debug)"
	@( cd EALib; $(MAKE) -s debug )
	@echo "### EALib library (debug) built."

EALib-conf::\
	EALib/Makefile

EALib/Makefile::\
	EALib/Makefile.in EALib/configure.in
	@printf "\n### Configuring EALib... "
	@( cd EALib; autoconf; ./configure 1>/dev/null 2>/dev/null )
	@echo "Done."

LinAlg::\
	LinAlg-conf
	@echo
	@echo "### Building LinAlg library"
	@( cd LinAlg; $(MAKE) -s )
	@echo "### LinAlg library built."

LinAlg-debug::\
	LinAlg-conf
	@echo
	@echo "### Building LinAlg library (debug)"
	@( cd LinAlg; $(MAKE) -s debug )
	@echo "### LinAlg library (debug) built."

LinAlg-conf::\
	LinAlg/Makefile

LinAlg/Makefile::\
	LinAlg/Makefile.in LinAlg/configure.in
	@printf "\n### Configuring LinAlg... "
	@( cd LinAlg; autoconf; ./configure 1>/dev/null 2>/dev/null )
	@echo "Done."

Metric::\
	Metric-conf
	@echo
	@echo "### Building Metric library"
	@( cd Metric; $(MAKE) -s )
	@echo "### Metric library built."

Metric-debug::\
	Metric-conf
	@echo
	@echo "### Building Metric library (debug)"
	@( cd Metric; $(MAKE) -s debug )
	@echo "### Metric library (debug) built."

Metric-conf::\
	Metric/Makefile

Metric/Makefile::\
	Metric/Makefile.in Metric/configure.in
	@printf "\n### Configuring Metric... "
	@( cd Metric; autoconf; ./configure 1>/dev/null 2>/dev/null )
	@echo "Done."

Mixture::\
	Mixture-conf
	@echo
	@echo "### Building Mixture library"
	@( cd Mixture; $(MAKE) -s )
	@echo "### Mixture library built."

Mixture-debug::\
	Mixture-conf
	@echo
	@echo "### Building Mixture library (debug)"
	@( cd Mixture; $(MAKE) -s debug )
	@echo "### Mixture library (debug) built."

Mixture-conf::\
	Mixture/Makefile

Mixture/Makefile::\
	Mixture/Makefile.in Mixture/configure.in
	@printf "\n### Configuring Mixture... "
	@( cd Mixture; autoconf; ./configure 1>/dev/null 2>/dev/null )
	@echo "Done."

ReClaM::\
	ReClaM-conf
	@echo
	@echo "### Building ReClaM library"
	@( cd ReClaM; $(MAKE) -s )
	@echo "### ReClaM library built."

ReClaM-debug::\
	ReClaM-conf
	@echo
	@echo "### Building ReClaM library (debug)"
	@( cd ReClaM; $(MAKE) -s debug )
	@echo "### ReClaM library (debug) built."

ReClaM-conf::\
	ReClaM/Makefile

ReClaM/Makefile::\
	ReClaM/Makefile.in ReClaM/configure.in
	@printf "\n### Configuring ReClaM... "
	@( cd ReClaM; autoconf; ./configure 1>/dev/null 2> /dev/null )
	@echo "Done."

Rng::\
	Rng-conf
	@echo
	@echo "### Building Rng library"
	@( cd Rng; $(MAKE) -s )
	@echo "### Rng library built."

Rng-debug::\
	Rng-conf
	@echo
	@echo "### Building Rng library (debug)"
	@( cd Rng; $(MAKE) -s debug )
	@echo "### Rng library (debug) built."

Rng-conf::\
	Rng/Makefile

Rng/Makefile::\
	Rng/Makefile.in Rng/configure.in
	@printf "\n### Configuring Rng... "
	@( cd Rng; autoconf; ./configure 1>/dev/null 2> /dev/null )
	@echo "Done."

TestData::\
	TestData-conf
	@echo
	@echo "### Building TestData libraries"
	@( cd TestData; $(MAKE) -s )
	@echo "### TestData libraries built."

TestData-debug::\
	TestData-conf
	@echo
	@echo "### Building TestData libraries (debug)"
	@( cd TestData; $(MAKE) -s debug )
	@echo "### TestData libraries (debug) built."

TestData-conf::\
	TestData/Makefile

TestData/Makefile::\
	TestData/Makefile.in TestData/configure.in
	@printf "\n### Configuring TestData... "
	@( cd TestData; autoconf; ./configure 1>/dev/null 2>/dev/null )
	@echo "Done."

examples:
	@-for i in $(LIBRARIES); do \
	    echo "### Building examples for" $$i; \
	    ( cd $$i; $(MAKE) -s examples ); \
	    echo "### Done."; \
	    echo ; \
	done

examples-debug:
	@-for i in $(LIBRARIES); do \
	    echo ; \
	    echo "### Building examples for" $$i "(debug)"; \
	    ( cd $$i; $(MAKE) -s examples-debug ); \
	    echo "### Done."; \
	done

man: rights
	@-for i in $(LIBRARIES); do \
	    echo ""; \
	    echo "### Building documentation for library" `echo $$i | cut -f1 -d/`; \
	    ( cd $$i; $(MAKE) -s man ); \
	    echo "### Done."; \
	done

doc:
	@-make selftest-man
	@-for i in $(REFERENCES); do \
	    echo ""; \
	    echo "### Building reference for library" `echo $$i | cut -f1 -d/`; \
	    ( cd $$i; $(MAKE) -s doc ); \
	    echo "### Done."; \
	done


depend: configure
	@echo
	@echo "### Establishing dependencies"
	@-for i in $(LIBRARIES); do \
	    ( cd $$i; $(MAKE) -s depend ); \
	done
	@echo "### Done."

independ: configure
	@echo
	@echo "### Removing dependencies"
	@-for i in $(LIBRARIES); do \
		( cd $$i; $(MAKE) -s independ ); \
	done
	@echo "### Done."

clean:
	@echo
	@-for i in $(LIBRARIES); do \
	    echo "### Cleaning directory" $$i; \
	    ( cd $$i; $(MAKE) -s clean ); \
	done
	@rm -fR win_tmp *~ config.* *.ncb *.opt *.log

distclean: clean
	@-for i in $(LIBRARIES); do \
		echo; \
	    echo "### Cleaning distribution" $$i; \
	    ( cd $$i; $(MAKE) -s distclean ); \
	    echo "### Done."; \
	done
	@rm -fR lib* man


changeversion:
	@echo Input the new version number:
	@(read ver;\
	liste=`egrep -R -e '^[ 	]{0,}[Vv][eE][rR]([sS][iI][oO][nN]){0,1}[ 	]{0,}=[ 	]{0,}([0-9]\.){0,}[0-9]{1,}([# 	]{0,}.{0,}$$)' *| cut -d":" -f1`; \
	for i in $$liste; do \
		echo "Processing the file $$i"; \
		sed -e "s/\(^[ 	]\{0,\}[vV][eE][rR]\([sS][iI][oO][nN]\)\{0,1\}[ 	]\{0,\}=[ 	]\{0,\}\)\([0-9]\{1,\}\.\)\{0,\}[0-9]\{1,\}\([# 	]\{0,\}.\{0,\}\$\)/\1$$ver\4/g" $$i >$$i-new;\
	done ;\
	for i in $$liste; do \
		mv -f $$i-new $$i;\
	done ;\
	)


ARCDIR		:= $(PACKAGE)-$(VERSION)

dist:   tar Pagode-package SimAnn-package SpreadCAT-package
	@printf " ## Archiving and compressing Shark... "
	@for i in `cat Admin/files_used`; do tar rf $(ARCDIR).tar $$i; done;
	@gzip --best $(ARCDIR).tar; mv $(ARCDIR).tar.gz $(ARCDIR).tgz
	@rm -f configure.in_others
	@rm -f Mixture/configure.in_others
	@echo "### Done."

zipdist:   tar SimAnn-zippackage
	@printf " ## Archiving and compressing Shark... "
	@for i in `cat Admin/files_used`; do zip -9 $(ARCDIR).zip $$i > \dev\null; done; 
	@rm -f configure.in_others
	@rm -f Mixture/configure.in_others
	@echo "### Done."

tar:	\
	unix2dos \
	independ

# Show differences of the files of the current
# distribution and the distribution stored in
# directory 'DIFF_DIR'.
# Only the distribution files as given in
# './Admin/files_used' are compared.
#DIFF_DIR = ~/source/Linux/modules/
DIFF_DIR = ~/SONN/sf-repository/Shark

admin_diff:
	@for i in `sed '/.*Metric/d' Admin/files_used | \
	sed '/.*\.dsp/d'                              | \
	sed '/.*\.dsw/d'                              | \
	sed '/.*Pagode/d'                             | \
	sed '/.*SimAnn/d'                             | \
	sed '/.*SpreadCAT/d'                          | \
	sed '/.*TestData/d'                           | \
	sed '/.*ZNconfig.h/d'                         | \
	sed '/.*ZNminmax.h/d'                         | \
	sed '/\.\/CHANGES\.txt/d'	              | \
	sed '/\.\/COMPATIBILITY\.txt/d'               | \
	sed '/\.\/COPYING\.txt/d'                     | \
	sed '/\.\/INDEX\.txt/d'                       | \
	sed '/\.\/INSTALL\.txt/d'                     | \
	sed '/\.\/cleanUp/d'                          | \
	sed '/\.\/configure\.in/d'                    | \
	sed '/\.\/Makefile\.in/d'                     | \
	sed '/\.\/make_/d'                            | \
	sed '/\.\/makedocs\.bat/d'                    | \
	sed '/files_used/d'                           | \
	sed '/\.\/shark_index\.html/d'                | \
	sed 's/\.\///'				      | \
	cat`;                                           \
	do echo diff ./$$i $(DIFF_DIR)/$$i;             \
        diff -s ./$$i $(DIFF_DIR)/$$i;                  \
        echo --------------------------------; done;

install: install-clean
	@echo
	@echo "### Installing Shark"
	@mkdir -p lib_$(ARCH_MACH)
	@mkdir -p lib_$(ARCH_MACH)/debug
	@mkdir -p man
	@printf " ## Library links... "
	@( for i in `$(FIND) $$PWD -type f -name "lib*_$(ARCH_MACH)${LIBEXT}"  -a \
	! -type l -print | grep -v debug`; do for j in `basename $$i | 	\
	sed 's/_$(ARCH_MACH)//'`; do ln -s $$i lib_$(ARCH_MACH)/$$j; done; \
	done )
	@( for i in `$(FIND) $$PWD -type f -name "lib*_$(ARCH_MACH)${LIBEXT}" -a \
	! -type l -print | grep -v debug`; do ln -s $$i lib_$(ARCH_MACH); \
	done )
	@echo "Done."
	@printf " ## Library (debug) links... "
	@( for i in `$(FIND) $$PWD -type f -name "lib*_$(ARCH_MACH)${LIBEXT}" -a \
	! -type l -print | grep debug`; do for j in `basename $$i | \
	sed 's/_$(ARCH_MACH)//'`; do ln -sf $$i lib_$(ARCH_MACH)/debug/$$j; \
	done; done )
	@( for i in `$(FIND) $$PWD -type f -name "lib*_$(ARCH_MACH)${LIBEXT}" -a \
	! -type l -print | grep debug`; do ln -sf $$i lib_$(ARCH_MACH)/debug; \
	done )
	@echo "Done."
	@printf " ## Man(ual) links... "
	@( for i in `$(FIND) $$PWD \( -type f -name "*.ps" -o \
	-name "*.ps.gz" \) -a ! -type l -print | grep -v Figures | \
	grep -v Images`; do ln -sf $$i man; done )
	@echo "Done."
	@echo "### Done."
	@echo

install-ini: install-clean
	@echo
	@echo "### Installing Shark (INI-version)"
	@mkdir -p lib_$(ARCH_MACH)
	@mkdir -p lib_$(ARCH_MACH)/debug
	@mkdir -p man
	@printf " ## Library links... "
	@( for i in `$(FIND) $$PWD -type f -name "lib*_$(ARCH_MACH)${LIBEXT}" -a \
	! -type l -print | grep -v debug | sed 's/_vs_[0-9]*/_vs/'`; \
	do for j in `basename $$i | sed 's/_$(ARCH_MACH)//'`; \
	do ln -s $$i lib_$(ARCH_MACH)/$$j; done; done )
	@( for i in `$(FIND) $$PWD -type f -name "lib*_$(ARCH_MACH)${LIBEXT}" -a \
	! -type l -print | grep -v debug | sed 's/_vs_[0-9]*/_vs/'`; \
	do ln -s $$i lib_$(ARCH_MACH); done )
	@printf " ## Library (debug) links... "
	@( for i in `$(FIND) $$PWD -type f -name "lib*_$(ARCH_MACH)${LIBEXT}" -a \
	! -type l -print | grep debug | sed 's/_vs_[0-9]*/_vs/'`; \
	do for j in `basename $$i | sed 's/_$(ARCH_MACH)//'`; \
	do ln -s $$i lib_$(ARCH_MACH)/debug/$$j; done; done )
	@( for i in `$(FIND) $$PWD -type f -name "lib*_$(ARCH_MACH)${LIBEXT}" -a \
	! -type l -print | grep debug | sed 's/_vs_[0-9]*/_vs/'`; \
	do ln -s $$i lib_$(ARCH_MACH)/debug; done )	
	@printf " ## Man(ual) links... "
	@( for i in `$(FIND) $$PWD \( -type f -name "*.ps" -o \
	-name "*.ps.gz" \) -a ! -type l -print | grep -v Figures | \
	grep -v Images`; do ln -s $$i man; done )
	@echo "Done."
	@echo "### Done."
	@echo


install-clean:
	@-rm -f lib_$(ARCH_MACH)/*${LIBEXT} lib_$(ARCH_MACH)/debug/* man/*

rights:
	@printf "\n### Setting rights... "
	@$(FIND) -type d -exec chmod 775 {} \;
	@$(FIND) -type f -exec chmod a+rX,ug+w {} \;
	@$(FIND) `$(FIND) \( -type d -name "exa*" -o -name "bin" \) \
	-o \( -type f -name "rlatex" -o -name "cleanUp" -o -name "configure" \)` -type f -exec file {} \; | \
	grep -i "\( script\)\|\(elf.*executable\)\|\(executable\)" | cut -f1 -d: | \
	xargs chmod 775
	@echo " Done!"

dos2unix:
	@printf "### Windows => Unix/Linux - Removing carriage returns from sourcefiles... "
	@for i in `$(FIND) . -type f \( \
	-name "*.h" -o 			\
	-name "*.cpp" -o		\
	-name "Makefile.in" -o  	\
	-name "configure.in" -o  	\
	-name "*.tex" -o		\
	-name "*.bib" -o		\
	-name "*.sty" -o		\
	-name "rlatex" -o		\
	-name "*.txt" 			\) -print`; \
	do mv $$i tmp ; tr -d "\r" < tmp > $$i ; rm -f tmp ; done
	@echo "Done."

dos2unix_all:
	@printf "### Windows => Unix/Linux - Removing carriage returns from all files... "
	@for i in `$(FIND) . -type f -name "*" -print`; do mv $$i tmp ; \
	tr -d "\r" < tmp > $$i ; rm -f tmp ; done
	@echo "Done."

unix2dos:
	@printf "### Transforming all Visual C++ project and workspace files to Windows format... "
	@for i in `$(FIND) . -type f \(\
	-name "*.dsp" -o            \
	-name "*.dsw"               \) -print`; \
	do $(PERL) -pi -e 's/\n/\r\n/;' $$i; done
	@echo "Done."


###############  package make commands for the MOO-EALib
#

ifdef FLAG_MOEALIB
MOO-EALib::\
	MOO-EALib-conf
	@echo
	@echo "### Building MOO-EALib library"
	@( cd MOO-EALib; $(MAKE) -s )
	@echo "### MOO-EALib library built."

MOO-EALib-debug::\
	MOO-EALib-conf
	@echo
	@echo "### Building MOO-EALib library (debug)"
	@( cd MOO-EALib; $(MAKE) -s debug )
	@echo "### MOO-EALib library (debug) built."

MOO-EALib-conf::\
	MOO-EALib/Makefile

MOO-EALib/Makefile::\
	MOO-EALib/Makefile.in MOO-EALib/configure.in
	@printf "\n### Configuring MOO-EALib... "
	@( cd MOO-EALib; autoconf; ./configure 1>/dev/null 2>/dev/null )
	@echo "Done."
else
MOO-EALib::
	@echo "### There is no copy of MOO-EALib wihtin your SHARK directory."

MOO-EALib-debug::    MOO-EALib
MOO-EALib-conf::     MOO-EALib
MOO-EALib/Makefile:: MOO-EALib
endif


###############  package make commands for Pagode
#

ifdef FLAG_PAGODE
Pagode::\
	Pagode-conf
	@echo
	@echo "### Building Pagode library"
	@( cd Pagode; $(MAKE) -s )
	@echo "### Pagode library built."

Pagode-debug::\
	Pagode-conf
	@echo
	@echo "### Building Pagode library (debug)"
	@( cd Pagode; $(MAKE) -s debug )
	@echo "### Pagode library (debug) built."

Pagode-conf::\
	Pagode/Makefile

Pagode/Makefile::\
	Pagode/Makefile.in Pagode/configure.in
	@printf "\n### Configuring Pagode... "
	@( cd Pagode; autoconf; ./configure 1>/dev/null 2> /dev/null )
	@echo "Done."

Pagode-package::  
	@printf " ## Archiving and compressing Pagode files... "
	@for i in `cat Admin/files_used_pagode`; do tar rf Pagode-$(VERSION).tar $$i; done; \
	gzip --best Pagode-$(VERSION).tar; mv Pagode-$(VERSION).tar.gz Pagode-$(VERSION).tgz
	@echo "### Done."

Pagode-zippackage::  
	@printf " ## Archiving and compressing Pagode files... "
	@for i in `cat Admin/files_used_pagode`; do zip -9 Pagode-$(VERSION).zip $$i > \dev\null; done;
	@echo "### Done."
else
Pagode::
	@echo "### There is no copy of Pagode wihtin your SHARK directory."

Pagode-debug::    Pagode
Pagode-conf::     Pagode
Pagode/Makefile:: Pagode
Pagode-package::  Pagode
Pagode-zippackage::  Pagode
endif



###############  package make commands for SimAnn
#

ifdef FLAG_SIMANN
SimAnn::\
	SimAnn-conf
	@echo
	@echo "### Building SimAnn library"
	@( cd SimAnn; $(MAKE) -s )
	@echo "### SimAnn library built."

SimAnn-debug::\
	SimAnn-conf
	@echo
	@echo "### Building SimAnn library (debug)"
	@( cd SimAnn; $(MAKE) -s debug )
	@echo "### SimAnn library (debug) built."

SimAnn-conf::\
	SimAnn/Makefile

SimAnn/Makefile::\
	SimAnn/Makefile.in SimAnn/configure.in
	@printf "\n### Configuring SimAnn... "
	@( cd SimAnn; autoconf; ./configure 1>/dev/null 2> /dev/null )
	@echo "Done."

SimAnn-package::  
	@printf " ## Archiving and compressing SimAnn files... "
	@for i in `cat Admin/files_used_simann`; do tar rf SimAnn-$(VERSION).tar $$i; done; \
	gzip --best SimAnn-$(VERSION).tar; mv SimAnn-$(VERSION).tar.gz SimAnn-$(VERSION).tgz
	@echo "### Done."

SimAnn-zippackage::  
	@printf " ## Archiving and compressing SimAnn files... "
	@for i in `cat Admin/files_used_simann`; do zip -9 SimAnn-$(VERSION).zip $$i > \dev\null; done;
	@echo "### Done."
else
SimAnn::
	@echo "### There is no copy of SimAnn wihtin your SHARK directory."

SimAnn-debug::    SimAnn
SimAnn-conf::     SimAnn
SimAnn/Makefile:: SimAnn
SimAnn-package::  SimAnn
SimAnn-zippackage::  SimAnn
endif



###############  package make commands for SpreadCat
#

ifdef FLAG_SPREADCAT
SpreadCAT::\
	SpreadCAT-conf
	@echo
	@echo "### Building SpreadCAT library"
	@( cd SpreadCAT; $(MAKE) -s )
	@echo "### SpreadCAT library built."

SpreadCAT-debug::\
	SpreadCAT-conf
	@echo
	@echo "### Building SpreadCAT library (debug)"
	@( cd SpreadCAT; $(MAKE) -s debug )
	@echo "### SpreadCAT library (debug) built."

SpreadCAT-conf::\
	SpreadCAT/Makefile

SpreadCAT/Makefile::\
	SpreadCAT/Makefile.in SpreadCAT/configure.in
	@printf "\n### Configuring SpreadCAT... "
	@( cd SpreadCAT; autoconf; ./configure 1>/dev/null 2> /dev/null )
	@echo "Done."

SpreadCAT-package::  
	@printf " ## Archiving and compressing SpreadCAT files... "
	@for i in `cat Admin/files_used_spreadcat`; do tar rf SpreadCAT-$(VERSION).tar $$i; done; \
	gzip --best SpreadCAT-$(VERSION).tar; mv SpreadCAT-$(VERSION).tar.gz SpreadCAT-$(VERSION).tgz
	@echo "### Done."

SpreadCAT-zippackage::  
	@printf " ## Archiving and compressing SpreadCAT files... "
	@for i in `cat Admin/files_used_spreadcat`; do zip -9 SpreadCAT-$(VERSION).zip $$i > \dev\null; done;
	@echo "### Done."
else
SpreadCAT::
	@echo "### There is no copy of SpreadCAT wihtin your SHARK directory."

SpreadCAT-debug::    SpreadCAT
SpreadCAT-conf::     SpreadCAT
SpreadCAT/Makefile:: SpreadCAT
SpreadCAT-package::  SpreadCAT
SpreadCAT-zippackage::  SpreadCAT
endif





###############  commands for selftest
#

ifdef FLAG_SELFTEST

AVAILABLE_TESTS =	Array-selftest \
			LinAlg-selftest \
			EALib-selftest \
			ReClaM-selftest \
			Metric-selftest \
			Rng-selftest \
			TestData-selftest

ifdef FLAG_SPREADCAT
	        AVAILABLE_TESTS +=	SpreadCAT-selftest
endif
ifdef FLAG_PAGODE
	        AVAILABLE_TESTS +=	Pagode-selftest
endif
ifdef FLAG_SIMANN
	        AVAILABLE_TESTS +=	SimAnn-selftest
endif
ifdef FLAG_MOEALIB
	        AVAILABLE_TESTS +=	MOO-EALib-selftest
endif


AVAILABLE_TESTS_INSTALL:= $(AVAILABLE_TESTS:%=%-install)

AVAILABLE_TESTS_FULL:= $(AVAILABLE_TESTS:%=%-full)

AVAILABLE_TESTS_CLEAN:= $(AVAILABLE_TESTS:%=%-clean)



selftest-install:  $(AVAILABLE_TESTS_INSTALL)

selftest:   $(AVAILABLE_TESTS)

selftest-full:	$(AVAILABLE_TESTS_FULL)

selftest-clean:  $(AVAILABLE_TESTS_CLEAN)

selftest-man: 
	@-(cd selftest-doc; make -s)

Array-selftest-install:
	@-(cd $(SELFTEST_DIR)/Array; make -s )

Array-selftest:
	@-(cd $(SELFTEST_DIR)/Array; make -s selftest)

Array-selftest-clean:
	@-(cd $(SELFTEST_DIR)/Array; make -s clean)

Array-selftest-full:
	@-(cd $(SELFTEST_DIR)/Array; make -s selftest-full)

LinAlg-selftest-install:
	@-(cd $(SELFTEST_DIR)/LinAlg; make -s )

LinAlg-selftest:
	@-(cd $(SELFTEST_DIR)/LinAlg; make -s selftest)

LinAlg-selftest-clean:
	@-(cd $(SELFTEST_DIR)/LinAlg; make -s clean)

LinAlg-selftest-full:
	@-(cd $(SELFTEST_DIR)/LinAlg; make -s selftest-full)

EALib-selftest-install:
	@-(cd $(SELFTEST_DIR)/EALib; make -s ) 

EALib-selftest:
	@-(cd $(SELFTEST_DIR)/EALib; make -s selftest)

EALib-selftest-clean:
	@-(cd $(SELFTEST_DIR)/EALib; make -s clean) 

EALib-selftest-full:
	@-(cd $(SELFTEST_DIR)/EALib; make -s selftest-full)

ReClaM-selftest-install:
	@-(cd $(SELFTEST_DIR)/ReClaM; make -s )

ReClaM-selftest:
	@-(cd $(SELFTEST_DIR)/ReClaM; make -s selftest)

ReClaM-selftest-clean:
	@-(cd $(SELFTEST_DIR)/ReClaM; make -s clean)

ReClaM-selftest-full:
	@-(cd $(SELFTEST_DIR)/ReClaM; make -s selftest-full)

Rng-selftest-install:
	@-(cd $(SELFTEST_DIR)/Rng; make -s ) 

Rng-selftest:
	@-(cd $(SELFTEST_DIR)/Rng; make -s selftest) 

Rng-selftest-clean:
	@-(cd $(SELFTEST_DIR)/Rng; make -s clean)

Rng-selftest-full:
	@-(cd $(SELFTEST_DIR)/Rng; make -s selftest-full) 

Metric-selftest-install:
	@-(cd $(SELFTEST_DIR)/Metric; make -s ) 

Metric-selftest:
	@-(cd $(SELFTEST_DIR)/Metric; make -s selftest) 

Metric-selftest-clean:
	@-(cd $(SELFTEST_DIR)/Metric; make -s clean)

Metric-selftest-full:
	@-(cd $(SELFTEST_DIR)/Metric; make -s selftest-full) 

TestData-selftest-install:
	@-(cd $(SELFTEST_DIR)/TestData; make -s )

TestData-selftest:
	@-(cd $(SELFTEST_DIR)/TestData; make -s selftest)

TestData-selftest-clean:
	@-(cd $(SELFTEST_DIR)/TestData; make -s clean) 

TestData-selftest-full:
	@-(cd $(SELFTEST_DIR)/TestData; make -s selftest-full)


ifdef FLAG_SPREADCAT
SpreadCAT-selftest-install:
	@-(cd $(SELFTEST_DIR)/SpreadCAT; make -s ) 

SpreadCAT-selftest:
	@-(cd $(SELFTEST_DIR)/SpreadCAT; make -s selftest)

SpreadCAT-selftest-clean:
	@-(cd $(SELFTEST_DIR)/SpreadCAT; make -s clean) 

SpreadCAT-selftest-full:
	@-(cd $(SELFTEST_DIR)/SpreadCAT; make -s selftest-full)

else
SpreadCAT-selftest-install::
	@echo "### There is no copy of SpreadCAT wihtin your SHARK directory."

SpreadCAT-selftest:    SpreadCAT-selftest-install
SpreadCAT-selftest-clean:    SpreadCAT-selftest-install
SpreadCAT-selftest-full:    SpreadCAT-selftest-install
endif



ifdef FLAG_PAGODE
Pagode-selftest-install:
	@-(cd $(SELFTEST_DIR)/Pagode; make -s ) 

Pagode-selftest:
	@-(cd $(SELFTEST_DIR)/Pagode; make -s selftest)

Pagode-selftest-full:
	@-(cd $(SELFTEST_DIR)/Pagode; make -s selftest-full)

Pagode-selftest-clean:
	@-(cd $(SELFTEST_DIR)/Pagode; make -s clean) 
else
Pagode-selftest-install::
	@echo "### There is no copy of Pagode wihtin your SHARK directory."

Pagode-selftest::    Pagode-selftest-install
Pagode-selftest-clean::    Pagode-selftest-install
Pagode-selftest-full::    Pagode-selftest-install
endif


ifdef FLAG_SIMANN
SimAnn-selftest-install:
	@-(cd $(SELFTEST_DIR)/SimAnn; make -s ) 

SimAnn-selftest:
	@-(cd $(SELFTEST_DIR)/SimAnn; make -s selftest)

SimAnn-selftest-full:
	@-(cd $(SELFTEST_DIR)/SimAnn; make -s selftest-full)

SimAnn-selftest-clean:
	@-(cd $(SELFTEST_DIR)/SimAnn; make -s clean) 
else
SimAnn-selftest-install::
	@echo "### There is no copy of SimAnn wihtin your SHARK directory."

SimAnn-selftest::    SimAnn-selftest-install
SimAnn-selftest-clean::    SimAnn-selftest-install
SimAnn-selftest-full::    SimAnn-selftest-install
endif


ifdef FLAG_MOEALIB
MOO-EALib-selftest-install:
	@-(cd $(SELFTEST_DIR)/MOO-EALib; make -s )

MOO-EALib-selftest:
	@-(cd $(SELFTEST_DIR)/MOO-EALib; make -s selftest)

MOO-EALib-selftest-full:
	@-(cd $(SELFTEST_DIR)/MOO-EALib; make -s selftest-full)

MOO-EALib-selftest-clean:
	@-(cd $(SELFTEST_DIR)/MOO-EALib; make -s clean)
else
MOO-EALib-selftest-install::
	@echo "### There is no copy of MOO-EALib wihtin your SHARK directory."

MOO-EALib-selftest::    MOO-EALib-selftest-install
MOO-EALib-selftest-clean::    MOO-EALib-selftest-install
MOO-EALib-selftest-full::    MOO-EALib-selftest-install
endif

endif