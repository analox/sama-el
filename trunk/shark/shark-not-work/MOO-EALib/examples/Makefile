# MOO_EALib examples makefile

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
FLAG_SPREADCAT		= 
DEFS 			= -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DYYTEXT_POINTER=1 -DX_DISPLAY_MISSING=1 
LDFLAGS 		= 
CC			= gcc
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
sharkdir		= $(shell cd ../..; pwd)
basedir			= ..
top_inc_dir		= $(sharkdir)/include
ealib_libdir		= $(sharkdir)/EALib/lib
mooealib_libdir		= ../lib
rng_libdir		= $(sharkdir)/Rng/lib
array_libdir		= $(sharkdir)/Array/lib
linalg_libdir		= $(sharkdir)/LinAlg/lib
reclam_libdir		= $(sharkdir)/ReClaM/lib

ifdef FLAG_SPREADCAT
	spreadcat_libdir	= $(sharkdir)/SpreadCAT/lib
endif

############### additional include directories
#
X_CFLAGS	   	+= -I$(top_inc_dir)	

X_LIBS_ND       	= -L$(ealib_libdir) -L$(mooealib_libdir) -L$(rng_libdir) -L$(reclam_libdir) -L$(array_libdir) -L$(linalg_libdir) 
ifdef FLAG_SPREADCAT
	X_LIBS_ND	+=  -L$(spreadcat_libdir)
endif
X_LIBS_D	= $(X_LIBS_ND:%/lib=%/lib/debug)
X_LIBS_D       := $(X_LIBS) $(X_LIBS_D)
X_LIBS_ND      := $(X_LIBS) $(X_LIBS_ND)

############### additional settings
#
LDLIBS         := -lEALib_$(ARCH_MACH) -lMOO-EALib_$(ARCH_MACH) -lRng_$(ARCH_MACH) -lReClaM_$(ARCH_MACH) -lArray_$(ARCH_MACH) -lLinAlg_$(ARCH_MACH) $(LDLIBS)
ifdef FLAG_SPREADCAT
	LDLIBS += -lSpreadCAT_$(ARCH_MACH) 
endif  
############## os specific settings
#
CFLAGS	       += 	-Wall
CFLAGS_ND   	= -O3 -DNDEBUG \
			-fomit-frame-pointer \
			-finline-functions \
			-ffast-math
CFLAGS_D		= -g -DDEBUG

ifeq ($(HOST),Linux)
    DEFS       += -D__LINUX__
endif

ifeq ($(HOST),SunOS)
    DEFS	   += -D__SOLARIS__
    X_LIBS_ND  := $(X_LIBS_ND) $(X_LIBS_ND:-L%=-R%)
	X_LIBS_D   := $(X_LIBS_D) $(X_LIBS_D:-L%=-R%)
endif

CFLAGS		   += $(X_CFLAGS) $(DEFS)

############### example progs
#
PROGS			= 	DWA_$(ARCH_MACH)		\
				MOGA_$(ARCH_MACH)		\
				NSGA2_$(ARCH_MACH)		\
				PAES_$(ARCH_MACH)		\
				VEGA_$(ARCH_MACH)		\
				NSGA2example_$(ARCH_MACH)	\
				realcoded-NSGA2_$(ARCH_MACH)

PROGS_SCAT		= 	DWA-SCAT_$(ARCH_MACH)		\
				MOGA-SCAT_$(ARCH_MACH)		\
				NSGA2-SCAT_$(ARCH_MACH)		\
				PAES-SCAT_$(ARCH_MACH)		\
				VEGA-SCAT_$(ARCH_MACH)

PROGS_DEBUG		= $(PROGS:%=%_debug)
PROGS_SCAT_DEBUG	= $(PROGS_SCAT:%=%_debug)
###############	implicit rules
#
%.o: %.cpp
	$(CXX) $(CFLAGS) $(CFLAGS_ND) -o $@ -c $<

%.debug.o: %.cpp
	$(CXX) $(CFLAGS) $(CFLAGS_D) -o $@ -c $<

ifdef FLAG_SPREADCAT
all: $(PROGS_SCAT) $(PROGS)
debug: $(PROGS_SCAT_DEBUG) $(PROGS_DEBUG)
else
all:  $(PROGS)
	@echo " ## SpreadCAT library not found, not all examples are installed"
debug: $(PROGS_DEBUG)
	@echo " ## SpreadCAT library not found, not all examples are installed"
endif
############### suffixes list
#
.SUFFIXES:
.SUFFIXES:	.o $(LIBEXT) .a .c .h

############### make program
#
%_$(ARCH_MACH): %.o
	$(CXX) $(X_LIBS_ND) -o $@ $< $(LDLIBS)

%_$(ARCH_MACH)_debug: %.debug.o
	$(CXX) $(X_LIBS_D) -o $@ $< $(LDLIBS)

clean:
	-rm -f *.o *~ raw red white pca

distclean: clean
	-rm -f $(PROGS) $(PROGS_DEBUG) $(PROGS_SCAT) $(PROGS_SCAT_DEBUG) Archive1.txt

dos2unix:
	@printf "### Removing carriage returns from all source files..."
	@for i in `find . -type f \( -name "*.cpp" -o -name "Makefile.in" \) -print`; do mv $$i tmp ; tr -d "\r" < tmp > $$i ; rm -f tmp ; done
	@echo "done."

unix2dos:
	@printf "### Transforming all Visual C++ project and workspace files to Windows format... "
	@for i in `$(FIND) $(BASE_DIR) -type f \(\
	-name "*.dsp" -o            \
	-name "*.dsw"               \) -print`; \
	do $(PERL) -pi -e 's/\n/\r\n/;' $$i; done
	@echo "Done."









