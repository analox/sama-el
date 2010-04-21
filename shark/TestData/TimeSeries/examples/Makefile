# TimeSeries examples makefile

############### determine the host type
#
HOST			= Linux
ARCH_MACH		= i686

############### variables set by configure script
#
LDLIBS 			= 
CPPFLAGS 		= 
X_CFLAGS 		= 
ifeq ($(HOST),SunOS)
	X_CFLAGS   := $(X_CFLAGS:-I%=-isystem %)
endif
X_LIBS 			= 
X_EXTRA_LIBS		= 
X_PRE_LIBS 		= 
DEFS 			= -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DYYTEXT_POINTER=1 -DX_DISPLAY_MISSING=1 
LDFLAGS 		= 
CC				= gcc
CPP 			= gcc -E
CXX 			= g++
MAKE 			= /usr/bin/make
FIND 			= /usr/bin/find
PERL 			= /usr/bin/perl
LIBEXT			= .so
CFLAGS_ND		= -Wall -O3 -DNDEBUG  -fomit-frame-pointer  -finline-functions -ffast-math
CFLAGS_D		= -Wall -g -DDEBUG
CFLAGS			= -D__LINUX__ -fpic
LD_SHARED		= -shared

############### directories
#
sharkdir		= $(shell cd ../../..; pwd)
basedir			= ..
top_inc_dir		= $(sharkdir)/include
libdir			= $(basedir)/lib


############### additional include directories
#
#X_CFLAGS	   += 		-I$(top_inc_dir)	

X_LIBS_ND       = -L$(libdir) -L../../../Rng/lib -L../../../LinAlg/lib/ -L../../../Array/lib/

X_LIBS_D	    = $(X_LIBS_ND:%/lib=%/lib/debug)
X_LIBS_D       := $(X_LIBS) $(X_LIBS_D)
X_LIBS_ND      := $(X_LIBS) $(X_LIBS_ND)

############### additional settings
#
LDLIBS         := -lLinAlg_$(ARCH_MACH) -lArray_$(ARCH_MACH) -lRng_$(ARCH_MACH) -lTimeSeries_$(ARCH_MACH) $(LDLIBS)

############## os specific settings
#
CFLAGS		   += $(X_CFLAGS) $(DEFS) -I$(top_inc_dir)

############### example progs
#
PROGS			= 	embeddingTest_$(ARCH_MACH) \
				lorenz84_$(ARCH_MACH)

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









