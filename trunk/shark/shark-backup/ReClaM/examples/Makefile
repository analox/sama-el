# ReClaM examples makefile

############### determine the host type
#
HOST			= Linux
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
FLAG_PAGODE		= 
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

ifndef LDLIBS
	USE_MOTIF	=	no
else
	USE_MOTIF	=	yes
endif

############### directories
#
sharkdir			= $(shell cd ../..; pwd)
top_inc_dir		= $(sharkdir)/include
basedir			= ..
reclam_libdir		= $(basedir)/lib
array_libdir		= $(sharkdir)/Array/lib
rng_libdir		= $(sharkdir)/Rng/lib
linalg_libdir		= $(sharkdir)/LinAlg/lib
mixture_libdir		= $(sharkdir)/Mixture/lib

ifdef FLAG_PAGODE
	pagode_libdir		= $(sharkdir)/Pagode/lib
	plotwidget_dir		= $(sharkdir)/Pagode/plot_widgets/plot_widgets
endif

############### additional include directories
#
ifeq ($(USE_MOTIF), yes)
	X_CFLAGS	   += -I$(top_inc_dir)
	X_LIBS_ND       = -L$(reclam_libdir) -L$(array_libdir) -L$(rng_libdir) -L$(linalg_libdir) -L$(mixture_libdir) 
ifdef FLAG_PAGODE
	X_CFLAGS	+= -I$(plotwidget_dir)
	X_LIBS_ND       += -L$(pagode_libdir)
endif
else
	X_CFLAGS	   += -I$(top_inc_dir)
	X_LIBS_ND       = -L$(reclam_libdir) -L$(array_libdir) -L$(rng_libdir) -L$(linalg_libdir) -L$(mixture_libdir)
endif

X_LIBS_D	    = $(X_LIBS_ND:%/lib=%/lib/debug)
X_LIBS_D       := $(X_LIBS) $(X_LIBS_D)
X_LIBS_ND      := $(X_LIBS) $(X_LIBS_ND)

############### additional settings
#
ifeq ($(USE_MOTIF), yes)
	LDLIBS         := -lArray_$(ARCH_MACH) -lRng_$(ARCH_MACH) -lReClaM_$(ARCH_MACH) -lMixture_$(ARCH_MACH) -lLinAlg_$(ARCH_MACH) $(LDLIBS)
ifdef FLAG_PAGODE
	LDLIBS += -lPagode_$(ARCH_MACH)     
endif	
else
	LDLIBS         := -lArray_$(ARCH_MACH) -lRng_$(ARCH_MACH) -lReClaM_$(ARCH_MACH) -lMixture_$(ARCH_MACH) -lLinAlg_$(ARCH_MACH) $(LDLIBS)
endif

############## os specific settings
#

CFLAGS		   += $(X_CFLAGS) $(DEFS) -I$(top_inc_dir)

############### example progs
#
ifeq ($(USE_MOTIF), yes)
	PROGS			=	simpleRBFNet_$(ARCH_MACH) \
					simpleFFNet_$(ARCH_MACH) \
					simpleRNNet_$(ARCH_MACH) \
					simpleFFNetSource_$(ARCH_MACH) \
					simpleMSERNNet_$(ARCH_MACH) \
					SVMtest_$(ARCH_MACH) \
					KernelOptimization_$(ARCH_MACH)
ifdef FLAG_PAGODE
	PROGS		+=  proben1BFGS_$(ARCH_MACH) proben1_$(ARCH_MACH)
endif
else
	PROGS			=	simpleRBFNet_$(ARCH_MACH) \
					simpleFFNet_$(ARCH_MACH) \
					simpleRNNet_$(ARCH_MACH) \
					simpleFFNetSource_$(ARCH_MACH) \
					simpleMSERNNet_$(ARCH_MACH) \
					SVMtest_$(ARCH_MACH) \
					KernelOptimization_$(ARCH_MACH)
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









