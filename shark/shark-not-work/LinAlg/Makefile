############### determine the host type
#
HOST			= Linux
ARCH_MACH		= i686

SHELL 			= /bin/sh

############### name of package
#
PACKAGE			= LinAlg
VERSION         	= 1.4.4

############### variables set by configure script
#
LDLIBS 			= 
CPPFLAGS 		= 
X_CFLAGS 		= 
X_LIBS 			= 
X_EXTRA_LIBS		= 
X_PRE_LIBS 		= 
DEFS 			= -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DYYTEXT_POINTER=1 -DX_DISPLAY_MISSING=1 
LDFLAGS 		= 
CC				= gcc
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
basedir			= .
srcdir			= $(basedir)/src
mandir			= $(basedir)/man
libdir			= $(basedir)/lib
includedir		= $(basedir)/include

sharkdir		= $(shell cd ..; pwd)
top_inc_dir		= $(sharkdir)/include

####### objects needed to generate package library
#
OBJS	      	=	$(srcdir)/bfgs_$(ARCH_MACH).o \
			$(srcdir)/bfgs2_$(ARCH_MACH).o \
			$(srcdir)/cblnsrch_$(ARCH_MACH).o \
			$(srcdir)/detsymm_$(ARCH_MACH).o \
			$(srcdir)/discrimAnalysis_$(ARCH_MACH).o \
			$(srcdir)/dlinmin_$(ARCH_MACH).o \
			$(srcdir)/eigenerr_$(ARCH_MACH).o \
			$(srcdir)/eigensort_$(ARCH_MACH).o \
			$(srcdir)/eigensymm_$(ARCH_MACH).o \
			$(srcdir)/eigensymmJacobi_$(ARCH_MACH).o \
			$(srcdir)/eigensymmJacobi2_$(ARCH_MACH).o \
			$(srcdir)/fft_$(ARCH_MACH).o \
			$(srcdir)/g_inverse_$(ARCH_MACH).o \
			$(srcdir)/invert_$(ARCH_MACH).o \
			$(srcdir)/linalg_$(ARCH_MACH).o \
			$(srcdir)/LinearClassifier_$(ARCH_MACH).o \
			$(srcdir)/linearRegress_$(ARCH_MACH).o \
			$(srcdir)/LinearRegression_$(ARCH_MACH).o \
			$(srcdir)/linmin_$(ARCH_MACH).o \
			$(srcdir)/lnsrch_$(ARCH_MACH).o \
			$(srcdir)/PCA_$(ARCH_MACH).o \
			$(srcdir)/rank_$(ARCH_MACH).o \
			$(srcdir)/rankDecomp_$(ARCH_MACH).o \
			$(srcdir)/svd_$(ARCH_MACH).o \
			$(srcdir)/svdrank_$(ARCH_MACH).o \
			$(srcdir)/svdsort_$(ARCH_MACH).o 
OBJS_DEBUG		= $(OBJS:%.o=%.debug.o)

############### additional include directories
#
#X_CFLAGS	   += -I$(top_inc_dir)

############## os specific settings
#
CFLAGS		   += $(X_CFLAGS) -I$(top_inc_dir) $(DEFS)
CXXFLAGS		= $(CFLAGS)

.PHONY:		all debug library library-debug examples examples-debug \
		man doc everything clean dist distclean depend independ \
		unix2dos

all: library

everything: library library-debug examples examples-debug man doc

library: lib $(libdir)/lib$(PACKAGE)_$(ARCH_MACH)${LIBEXT}

lib:
	@-mkdir -p $(libdir)

$(libdir)/lib$(PACKAGE)_$(ARCH_MACH)${LIBEXT}: $(OBJS)
		$(CXX) $(LD_SHARED) $(X_LIBS) -o $@ $(OBJS)

debug: library-debug

library-debug:	$(libdir)/debug $(libdir)/debug/lib$(PACKAGE)_$(ARCH_MACH)${LIBEXT}

$(libdir)/debug:
	@-mkdir -p $(libdir)/debug

$(libdir)/debug/lib$(PACKAGE)_$(ARCH_MACH)${LIBEXT}: $(OBJS_DEBUG)
		$(CXX) $(LD_SHARED) $(X_LIBS) -o $@ $(OBJS_DEBUG)

examples: library
	@printf " ## Building examples for library LinAlg... "
	@( cd examples; $(MAKE) -s )
	@echo "Done."

examples-debug: library-debug
	@printf " ## Building examples for library LinAlg (debug)... "
	@( cd examples; $(MAKE) -s debug )
	@echo "Done."

man:
	@echo " ## No manual available yet."
#	( cd man; $(MAKE) -s )

doc:
	@printf " ## Creating reference for LinAlg library... "
	@( cd doc; $(MAKE) -s )
	@echo "Done."

###############	clean up
#
clean:
		-$(RM) $(OBJS) $(OBJS_DEBUG)
		-$(RM) Makefile.bak
		( cd examples; $(MAKE) -s clean )
		( cd doc; $(MAKE) -s clean )

distclean: clean
		-$(RM) $(libdir)/lib$(PACKAGE)_$(ARCH_MACH)${LIBEXT}
		-$(RM) $(libdir)/debug/lib$(PACKAGE)_$(ARCH_MACH)${LIBEXT}
		( cd examples; $(MAKE) -s distclean )
		( cd doc; $(MAKE) -s distclean )

ARCNAME := $(PACKAGE)-$(VERSION)

dist: clean
	@echo 
	@echo "### Creating tgz distribution archive"
	@-$(FIND) . \( -name "*.o"   -o \
			    -name "*${LIBEXT}"  -o \
			    -name "*.a"   -o \
				-name "*tab*" -o \
				-name "*~"    -o \
				-name "config.*"    -o \
				-name "Makefile"    -o \
				-name "configure"    -o \
				-name "*.tgz"    -o \
				-name "*.tar.bz2"    -o \
				-name "*.tar"    -o \
	            -name "*%"    -o \
	            -name "*.bak" \) | cut -c3- >> exclude
	@-$(FIND) . -type d \( -name "lib" -o  -name "debug" \) | cut -c3- >> exclude
	@-(cd ..; tar -chf- $(PACKAGE)/* -X $(PACKAGE)/exclude | gzip - > $(ARCNAME).tgz )
	@-mv ../$(ARCNAME).tgz .
	@-rm -f exclude
	@echo "### Done."
	@echo

dist-bzip2: clean
	@echo 
	@echo "### Creating bzip2 distribution archive"
	@-$(FIND) . \( -name "*.o"   -o \
			    -name "*${LIBEXT}"  -o \
			    -name "*.a"   -o \
				-name "*tab*" -o \
				-name "*~"    -o \
				-name "config.*"    -o \
				-name "Makefile"    -o \
				-name "configure"    -o \
				-name "*.tgz"    -o \
				-name "*.tar.bz2"    -o \
				-name "*.tar"    -o \
	            -name "*%"    -o \
	            -name "*.bak" \) | cut -c3- >> exclude
	@-$(FIND) . -type d \( -name "lib" -o  -name "debug" \) | cut -c3- >> exclude
	@-(cd ..; tar -chf- $(PACKAGE)/* -X $(PACKAGE)/exclude | bzip2 > $(ARCNAME).tar.bz2 )
	@-mv ../$(ARCNAME).tar.bz2 .
	@-rm -f exclude
	@echo "### Done."
	@echo


unix2dos:
	@printf "### Transforming all Visual C++ project and workspace files to Windows format... "
	@for i in `$(FIND) $(BASE_DIR) -type f \(\
	-name "*.dsp" -o            \
	-name "*.dsw"               \) -print`; \
	do $(PERL) -pi -e 's/\n/\r\n/;' $$i; done
	@echo "Done."

###############	dependencies
#
depend: dep_text
	@makedepend -o_$(ARCH_MACH).o -Y -I$(includedir) -I$(top_inc_dir) \
	$(srcdir)/*.cpp 2> /dev/null
	@makedepend -a -o_$(ARCH_MACH).debug.o -Y -I$(includedir) \
	-I$(top_inc_dir) $(srcdir)/*.cpp 2> /dev/null
	@echo "Done."


independ: dep_text
	@makedepend
	@echo "Done."

dep_text:
	@printf " ## LinAlg library dependencies... "

.SUFFIXES:	.o ${LIBEXT} .a .c .cpp .h

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


