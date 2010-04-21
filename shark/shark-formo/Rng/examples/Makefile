# Rng examples makefile

############### determine the host type
#
HHOST			= Linux
ARCH_MACH		= i686

############### variables set by configure script
#
LDLIBS 			= 
CPPFLAGS 		= 
X_CFLAGS 		=  -I/usr/X11R6/include
ifeq ($(HOST),SunOS)
	X_CFLAGS   := $(X_CFLAGS:-I%=-isystem %)
endif
X_LIBS 			=  -L/usr/X11R6/lib
X_EXTRA_LIBS		= 
X_PRE_LIBS 		=  -lSM -lICE
DEFS 			= -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" 
LDFLAGS 		= 
CC			= gcc
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
sharkdir		= $(shell cd ../..; pwd)
top_includedir		= $(sharkdir)/include
basedir			= ..
libdir			= $(basedir)/lib


############### additional include directories
#
#X_CFLAGS	+=	-I$(top_includedir)	
X_LIBS_ND	 =	-L$(libdir)

X_LIBS_D	= $(X_LIBS_ND:%/lib=%/lib/debug)
X_LIBS_D       := $(X_LIBS) $(X_LIBS_D)
X_LIBS_ND      := $(X_LIBS) $(X_LIBS_ND)

############### additional settings
#
LDLIBS         := -lRng_$(ARCH_MACH)

############## os specific settings
#

CFLAGS		   += $(X_CFLAGS) -I$(top_includedir) $(DEFS)

############### example progs
#
PROGS			= 	RngTest_$(ARCH_MACH) \
				RngTest2_$(ARCH_MACH) \
				RngTestSimple_$(ARCH_MACH)

PROGS_DEBUG		= $(PROGS:%=%_debug)

###############	implicit rules
#
%.o: %.cpp
	$(CXX) $(CFLAGS) $(CFLAGS_ND) -o $@ -c $<

%.debug.o: %.cpp
	$(CXX) $(CFLAGS) $(CFLAGS_D) -o $@ -c $<

all: $(PROGS)

debug: $(PROGS_DEBUG)

############### suffixes list
#
.SUFFIXES:
.SUFFIXES:	.o ${LIBEXT} .a .c .h

############### make program
#
%_$(ARCH_MACH): %.o
	$(CXX) $(X_LIBS_ND) -o $@ $< $(LDLIBS)

%_$(ARCH_MACH)_debug: %.debug.o
	$(CXX) $(X_LIBS_D) -o $@ $< $(LDLIBS)

clean:
	-rm -f *.o *~ raw red white pca

distclean: clean
	-rm -f $(PROGS) $(PROGS_DEBUG)

unix2dos:
	@printf "### Transforming all Visual C++ project and workspace files to Windows format... "
	@for i in `$(FIND) $(BASE_DIR) -type f \(\
	-name "*.dsp" -o            \
	-name "*.dsw"               \) -print`; \
	do $(PERL) -pi -e 's/\n/\r\n/;' $$i; done
	@echo "Done."









