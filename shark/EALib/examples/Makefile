# EALib examples makefile

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
 X_CFLAGS		:= $(X_CFLAGS:-I%=-isystem %)
endif
X_LIBS 			= 
X_EXTRA_LIBS		= 
X_PRE_LIBS 		= 
FLAG_PAGODE		= 
FLAG_SPREADCAT		= 
PVM_ROOT 		= 
LIBS_PVM 		= 
PVM_LDLIBS		= 
DEFS 			= -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DYYTEXT_POINTER=1 -DX_DISPLAY_MISSING=1 
LDFLAGS 		= 
CC			= gcc
CPP 			= gcc -E
CXX 			= g++
MAKE			= /usr/bin/make
FIND 			= /usr/bin/find
PERL 			= /usr/bin/perl
LIBEXT			= .so
CFLAGS_ND		= -Wall -O3 -DNDEBUG  -fomit-frame-pointer  -finline-functions -ffast-math
CFLAGS_D		= -Wall -g -DDEBUG
CFLAGS			= -D__LINUX__ -fpic
LD_SHARED		= -shared

ifndef LDLIBS
	USE_MOTIF	=	no
else
	USE_MOTIF	=	yes
endif

############### directories
#
sharkdir		= $(shell cd ../..; pwd)
basedir			= ..
top_inc_dir		= $(sharkdir)/include
ealib_libdir		= $(basedir)/lib
array_libdir		= $(sharkdir)/Array/lib
rng_libdir		= $(sharkdir)/Rng/lib
linalg_libdir		= $(sharkdir)/LinAlg/lib
reclam_libdir		= $(sharkdir)/ReClaM/lib

ifdef FLAG_PAGODE
 pagode_libdir		= $(sharkdir)/Pagode/lib
 widgetdir		= $(sharkdir)/Pagode/plot_widgets/plot_widgets
endif

ifdef FLAG_SPREADCAT
 spreadcat_libdir	= $(sharkdir)/SpreadCAT/lib
endif

############### additional include directories
#

ifeq ($(USE_MOTIF), yes) 
 X_LIBS_ND	 = -L$(ealib_libdir) -L$(array_libdir) -L$(rng_libdir) -L$(linalg_libdir) -L$(reclam_libdir)

 ifdef FLAG_PAGODE
  X_LIBS_ND	+=  -L$(pagode_libdir)
 endif

 ifdef FLAG_SPREADCAT
  X_LIBS_ND	+=  -L$(spreadcat_libdir)
 endif
else
  X_LIBS_ND	= -L$(ealib_libdir) -L$(array_libdir) -L$(rng_libdir) -L$(linalg_libdir) -L$(reclam_libdir)
endif

ifdef LIBS_PVM 
 X_CFLAGS       += -I${PVM_ROOT}/include
 X_LIBS_ND      += $(LIBS_PVM)
endif

X_LIBS_D	= $(X_LIBS_ND:%/lib=%/lib/debug)
X_LIBS_D       := $(X_LIBS) $(X_LIBS_D)
X_LIBS_ND      := $(X_LIBS) $(X_LIBS_ND)

############### additional settings
#

ifdef LIBS_PVM 
 LDLIBS		:= -lArray_$(ARCH_MACH) -lRng_$(ARCH_MACH) -lLinAlg_$(ARCH_MACH) -lEALibPVM_$(ARCH_MACH) -lReClaM_$(ARCH_MACH) $(LDLIBS)

 LDLIBS		+= $(PVM_LDLIBS)
else
 LDLIBS   	:= -lArray_$(ARCH_MACH) -lRng_$(ARCH_MACH) -lLinAlg_$(ARCH_MACH) -lEALib_$(ARCH_MACH) -lReClaM_$(ARCH_MACH) $(LDLIBS)
endif

ifeq ($(USE_MOTIF), yes)
 ifdef FLAG_PAGODE
  LDLIBS 	+= -lPagode_$(ARCH_MACH)
 endif

 ifdef FLAG_SPREADCAT
  LDLIBS 	+= -lSpreadCAT_$(ARCH_MACH)
 endif
endif

##############  specific settings
#

ifdef LIBS_PVM 
 DEFS		+= -DPVM_EXISTS=1
endif

CFLAGS		+= -I$(top_inc_dir) $(DEFS) $(X_CFLAGS)

############### example progs
#

ifeq ($(USE_MOTIF), yes) 
 PROGS		= 	ackleyES_$(ARCH_MACH) \
			countingOnes_$(ARCH_MACH) \
			integerES_$(ARCH_MACH) \
			sphereGA_$(ARCH_MACH) \
			sphereOstermeier_$(ARCH_MACH) \
			steadyState_$(ARCH_MACH) \
			cmaExample_$(ARCH_MACH) \
			cmaplus_$(ARCH_MACH)

 ifdef FLAG_SPREADCAT
  ifdef FLAG_PAGODE
   PROGS	+=	esAdap_$(ARCH_MACH)
  endif
 endif
else
 PROGS		= 	ackleyES_$(ARCH_MACH) \
			countingOnes_$(ARCH_MACH) \
			integerES_$(ARCH_MACH) \
			sphereGA_$(ARCH_MACH) \
			sphereOstermeier_$(ARCH_MACH) \
			steadyState_$(ARCH_MACH) \
	                cmaExample_$(ARCH_MACH) \
			cmaplus_$(ARCH_MACH)

endif

ifdef LIBS_PVM 	
 PROGS		+= 	pvmSphereGA_$(ARCH_MACH)
endif

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









