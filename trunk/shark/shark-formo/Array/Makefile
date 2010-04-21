############### determine the host type
#
HOST			= Linux
ARCH_MACH		= i686

SHELL 			= /bin/sh

############### name of package
#
PACKAGE			= Array
PACKAGE_MD		= $(PACKAGE)_$(ARCH_MACH)
VERSION         	= 1.4.4

############### variables set by configure script
#
LDLIBS 			= 
CPPFLAGS 		= 
X_CFLAGS 		=  -I/usr/X11R6/include
X_LIBS 			=  -L/usr/X11R6/lib
FLAG_PAGODE		= 
LIBS_PVM 		= 
X_EXTRA_LIBS		= 
X_PRE_LIBS 		=  -lSM -lICE
DEFS 			= -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" 
LDFLAGS 		= 
CC			= gcc
CPP 			= gcc -E
CXX 			= g++
FIND			= /usr/bin/find
MAKE			= /usr/bin/make
PERL			= /usr/bin/perl
LIBEXT			= .so
CFLAGS_ND		= -Wall -O3 -DNDEBUG  -fomit-frame-pointer  -finline-functions -ffast-math
CFLAGS_D		= -Wall -g -DDEBUG
CFLAGS			= -D__LINUX__ -fpic
LD_SHARED		= -shared

############### directories
#
sharkdir		= $(shell cd ..; pwd)
shark_includedir	= $(sharkdir)/include

basedir			= .
srcdir			= $(basedir)/src
mandir			= $(basedir)/man
libdir			= $(basedir)/lib
includedir		= $(basedir)/include

####### objects needed to generate package library
#
OBJS	      		= $(srcdir)/ReadArray_$(ARCH_MACH).o
OBJS_DEBUG		= $(OBJS:%.o=%.debug.o)

############### additional include directories
#
CFLAGS		   += 	-I$(shark_includedir) $(DEFS)
CXXFLAGS	    = 	$(CFLAGS)


.PHONY:		all debug library library-debug examples examples-debug \
		man doc everything clean dist distclean depend independ \
		dos2unix distfilter

all: \
	library


everything: 		\
	library 	\
	library-debug 	\
	examples 	\
	examples-debug	\
	man		\
	doc


library: 	\
	lib 	\
	$(libdir)/lib$(PACKAGE_MD)$(LIBEXT)


lib:
	@-mkdir -p $(libdir)


$(libdir)/lib$(PACKAGE_MD)$(LIBEXT): $(OBJS)
	$(CXX) $(LD_SHARED) $(X_LIBS) -o $@ $(OBJS)


debug: \
	library-debug

library-debug: \
	$(libdir)/debug $(libdir)/debug/lib$(PACKAGE_MD)$(LIBEXT)


$(libdir)/debug: 
	@-mkdir -p $(libdir)/debug


$(libdir)/debug/lib$(PACKAGE_MD)$(LIBEXT): $(OBJS_DEBUG)
	$(CXX) $(LD_SHARED) $(X_LIBS) -o $@ $(OBJS_DEBUG)


examples: \
	library
	@printf " ## Building examples for library Array... "
	@( cd examples; $(MAKE) -s )
	@echo "Done."


examples-debug: \
	library-debug
	@printf " ## Building examples for library Array (debug)... "
	@( cd examples; $(MAKE) -s debug )
	@echo "Done."


man:
	@printf " ## Creating documentation for Array library... "
	@( cd man; $(MAKE) -s )
	@echo "Done."

doc:
	@printf " ## Creating reference for Array library... "
	@( cd doc; $(MAKE) -s )
	@echo "Done."

dos2unix:
	@printf "### Removing carriage returns from all source files..."
	@for i in `find . -type f \( -name "*.cpp" -o -name "*.h" -o -name "Makefile.in" -o -name "configure.in" \) -print`; do mv $$i tmp ; tr -d "\r" < tmp > $$i ; rm -f tmp ; done
	@echo "done."
	@echo "### Removing carriage returns in directory \"examples\"..."
	@( cd examples; $(MAKE) -s dos2unix )
	@echo "### Removing carriage returns in directory \"man\"..."
	@( cd man; $(MAKE) -s dos2unix )

####### transform Visual C++ project- and workspace-files to dos-format by
#       adding carriage returns if necessary
#
unix2dos:
	@printf "### Transforming all Visual C++ project and workspace files to Windows format... "
	@for i in `$(FIND) $(BASE_DIR) -type f \(\
	-name "*.dsp" -o            \
	-name "*.dsw"               \) -print`; \
	do $(PERL) -pi -e 's/\n/\r\n/;' $$i; done
	@echo "Done."


###############	clean up
#

clean:
	@-$(RM) $(OBJS) $(OBJS_DEBUG)
	@-$(RM) Makefile.bak
	@( cd examples; $(MAKE) -s clean )
	@( cd man; $(MAKE) -s clean )
	@( cd doc; $(MAKE) -s clean )


distclean: \
	clean
	@-$(RM) $(libdir)/lib$(PACKAGE_MD)$(LIBEXT)
	@-$(RM) $(libdir)/debug/lib$(PACKAGE_MD)$(LIBEXT)
	@( cd examples; $(MAKE) -s distclean )
	@( cd man; $(MAKE) -s distclean )
	@( cd doc; $(MAKE) -s distclean )


distfilter: \
	distclean 
	@-$(FIND) . \( 	-name "*.o"   -o 	\
			-name "*$(LIBEXT)"  -o 	\
			-name "*.a"   -o 	\
			-name "*tab*" -o 	\
			-name "*~"    -o 	\
			-name "*%"    -o 	\
			-name "Makefile" -o 	\
			-name "configure" -o 	\
			-name "config.*" -o 	\
			-name "*.dvi" -o 	\
			-name "*.ps" -o 	\
			-name "*.log" -o 	\
			-name "*.aux" -o 	\
			-name "*.toc" -o 	\
			-name "exclude" -o 	\
			-name "*.tgz"    -o \
			-name "*.tar.bz2"    -o \
			-name "*.tar"    -o \
			-name "Array.shared" -o \
			-name "*.bak" \) | cut -c3- >> exclude
	@-$(FIND) . -type d \( -name "lib" -o  -name "debug" -o -name "CVS" \) | cut -c3- >> exclude


ARCNAME	:= $(PACKAGE)-$(VERSION)

dist: 	\
	distfilter
	@echo 
	@echo "### Creating tgz distribution archive"
	@-(cd ..; tar -chf- Array/* -X Array/exclude | gzip - > $(ARCNAME).tgz )
	@-mv ../$(ARCNAME).tgz .
	@-rm -f exclude
	@echo "### Done."
	@echo


dist-bzip2: \
	distfilter
	@echo 
	@echo "### Creating bzip2 distribution archive"
	@-(cd ..; tar -chf- Array/* -X Array/exclude | bzip2 > $(ARCNAME).tar.bz2 )
	@-mv ../$(ARCNAME).tar.bz2 .
	@-rm -f exclude
	@echo "### Done."
	@echo


###############	dependencies
#
depend:
	@printf " ## Array library dependencies... "
	@makedepend -o_$(ARCH_MACH).o -Y -I$(includedir) -I$(top_inc_dir) \
	$(srcdir)/*.cpp 2> /dev/null
	@makedepend -a -o_$(ARCH_MACH).debug.o -Y -I$(includedir) \
	-I$(top_inc_dir) $(srcdir)/*.cpp 2> /dev/null
	@echo "Done."

independ:
	@printf " ## Array library dependencies... "
	@makedepend
	@echo "Done."


.SUFFIXES:	.o $(LIBEXT) .a .c .cpp .h

###############	implicit rules
#
%_$(ARCH_MACH).o: %.c
	$(COMPILE.cc) $(CFLAGS_ND) -o $@ $<

%_$(ARCH_MACH).o: %.cpp
	$(COMPILE.cc) $(CFLAGS_ND) -o $@ $<

%_$(ARCH_MACH).debug.o:	%.c
	$(COMPILE.c) $(CFLAGS_D) -o $@ $<

%_$(ARCH_MACH).debug.o:	%.cpp
	$(COMPILE.cc) $(CFLAGS_D) -o $@ $<


# DO NOT DELETE
