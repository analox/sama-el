############### determine the host type
#
HOST			= Linux
ARCH_MACH		= i686

SHELL 			= /bin/sh

############### name of package
#
PACKAGE			= TestData
VERSION         	= 1.4.4

############### variables set by configure script
#
CPPFLAGS 		= 
X_CFLAGS 		=  -I/usr/X11R6/include
X_LIBS 			=  -L/usr/X11R6/lib
X_EXTRA_LIBS		= 
X_PRE_LIBS 		=  -lSM -lICE
DEFS 			= -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" 
LDFLAGS 		= 
CC				= gcc
CPP 			= gcc -E
CXX 			= g++
MAKE			= /usr/bin/make
FIND			= /usr/bin/find
PERL			= /usr/bin/perl
LIBEXT			= .so
CFLAGS_ND		= -Wall -O3 -DNDEBUG  -fomit-frame-pointer  -finline-functions -ffast-math
CFLAGS_D		= -Wall -g -DDEBUG
CFLAGS			= -D__LINUX__ -fpic
LD_SHARED		= -shared

############### directories
#
sharkdir		= $(shell cd ..; pwd)
top_inc_dir		= $(sharkdir)/include

basedir			= .
srcdir			= $(basedir)/src
mandir			= $(basedir)/man
libdir			= $(basedir)/lib
includedir		= $(basedir)/include

############### make commands
#
.PHONY:		all debug library library-debug examples examples-debug \
		man everything clean dist distclean depend independ \
		win2unix unix2dos


all: library

everything: library library-debug examples examples-debug man

library:
	@( cd RandomDistr; $(MAKE) -s ) 
	@( cd TimeSeries; $(MAKE) -s ) 

debug: library-debug

library-debug:	
	( cd RandomDistr; $(MAKE) -s library-debug )
	( cd TimeSeries; $(MAKE) -s library-debug )	

examples: library
	@echo " ## Building examples for TestData libraries... "
	@( cd RandomDistr; $(MAKE) -s examples )
	@( cd TimeSeries; $(MAKE) -s examples )
	@echo " ## Done."

examples-debug: library-debug
	@echo " ## Building examples for TestData libraries (debug)... "
	@( cd RandomDistr; $(MAKE) -s examples-debug )
	@( cd TimeSeries; $(MAKE) -s examples-debug )
	@echo " ## Done."

man:
	@echo " ## No manual available yet."
#	@echo " ## ... Omitted."

###############	clean up
#
clean:
	( cd RandomDistr; $(MAKE) -s clean )
	( cd TimeSeries; $(MAKE) -s clean )

distclean: clean
	( cd RandomDistr; $(MAKE) -s distclean )
	( cd TimeSeries; $(MAKE) -s distclean )

win2unix:
	@printf "### Windows => Unix/Linux - Removing carriage returns from all files..."
	@for i in `find . -type f -name "*" -print`; do mv $$i tmp ; tr -d "\r" < tmp > $$i ; rm -f tmp ; done
	@echo "done."

unix2dos:
	@printf "### Transforming all Visual C++ project and workspace files to Windows format... "
	@for i in `$(FIND) $(BASE_DIR) -type f \(\
	-name "*.dsp" -o            \
	-name "*.dsw"               \) -print`; \
	do $(PERL) -pi -e 's/\n/\r\n/;' $$i; done
	@echo "Done."

newconf:
	rm -f config.*

confclean:
	newconf
	rm -f configure
	rm -f Makefile

ARCNAME := $(PACKAGE)-$(VERSION)

dist: clean
	@echo 
	@echo "### Creating tgz distribution archive"
	@( cd RandomDistr; $(MAKE) -s dist )
	@( cd RandomDistr; mv RandomDistr*.tgz ../ )
	@echo 
	@echo 
	@( cd TimeSeries; $(MAKE) -s dist )
	@( cd TimeSeries; mv  TimeSeries*.tgz ../ )
	@echo "### Done."
	@echo

dist-bzip2: clean
	@echo 
	@echo "### Creating bzip2 distribution archive"
	@( cd RandomDistr; $(MAKE) -s dist-bzip2 )
	@( cd RandomDistr; mv RandomDistr*.bz2 ../ )
	@echo 
	@echo 
	@( cd TimeSeries; $(MAKE) -s dist-bzip2 )
	@( cd TimeSeries; mv  TimeSeries*.bz2 ../ )
	@echo "### Done."
	@echo

###############	dependencies
#
depend: 
	( cd RandomDistr; $(MAKE) -s depend )
	( cd TimeSeries; $(MAKE) -s depend )

independ:
	( cd RandomDistr; $(MAKE) -s independ )
	( cd TimeSeries; $(MAKE) -s independ )


# DO NOT DELETE



