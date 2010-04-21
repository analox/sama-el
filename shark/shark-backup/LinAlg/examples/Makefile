# LinAlg examples makefile

HOST			= Linux
ARCH_MACH		= i686

SHELL 			= /bin/sh

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
top_inc_dir		= $(sharkdir)/include

basedir			= ..
includedir		= $(top_inc_dir)
libdir_linalg		= $(sharkdir)/LinAlg/lib
libdir_array		= $(sharkdir)/Array/lib
libdir_rng		= $(sharkdir)/Rng/lib


############### additional include directories
#
#X_CFLAGS	   += -I$(top_inc_dir)
X_LIBS_ND       = -L$(libdir_linalg) -L$(libdir_array) -L$(libdir_rng)
X_LIBS_D	    = -L$(libdir_linalg)/debug -L$(libdir_array)/debug -L$(libdir_rng)/debug

############### additional settings
#
LDLIBS         := -lLinAlg_$(ARCH_MACH) -lArray_$(ARCH_MACH) -lRng_$(ARCH_MACH)

############## os specific settings
#
DEFS		   += -DFUNCPROTO -DXTFUNCPROTO -DUSE_DIRENT

CFLAGS		   += $(X_CFLAGS) -I$(top_inc_dir) $(DEFS)

############### example progs
#
PROGS			= pca_test_$(ARCH_MACH) \
			  rank_test_$(ARCH_MACH) \
			  rank_decomp_test_$(ARCH_MACH) \
			  svd_test_$(ARCH_MACH) \
			  svdrank_test_$(ARCH_MACH) \
			  svdsort_test_$(ARCH_MACH) \
			  eigensymm_test_$(ARCH_MACH) \
			  eigensymmJacobi_test_$(ARCH_MACH) \
			  eigensymmJacobi2_test_$(ARCH_MACH) \
			  eigensort_test_$(ARCH_MACH) \
			  eigenerr_test_$(ARCH_MACH) \
			  bfgs_test_$(ARCH_MACH) \
			  bfgs2_test_$(ARCH_MACH) \
			  g_inverse_matrix_$(ARCH_MACH) \
			  detsymm_test_$(ARCH_MACH) \
			  lnsrch_test_$(ARCH_MACH) \
			  cblnsrch_test_$(ARCH_MACH) \
			  linmin_test_$(ARCH_MACH) \
			  dlinmin_test_$(ARCH_MACH) \
			  linalg_simple_test_$(ARCH_MACH) \
			  covar_corrcoef_test_$(ARCH_MACH) \
                          linearRegression_test_$(ARCH_MACH) \
			  linearClassifier_test_$(ARCH_MACH) \
			  fft_test_$(ARCH_MACH)


PROGS_DEBUG		= $(PROGS:%=%_debug)

all: $(PROGS)

debug: $(PROGS_DEBUG)

###############	implicit rules
#
%.o: %.cpp
	$(CXX) $(CFLAGS) $(CFLAGS_ND) -o $@ -c $<

%.debug.o: %.cpp
	$(CXX) $(CFLAGS) $(CFLAGS_D) -o $@ -c $<

############### suffixes list
#
.SUFFIXES:
.SUFFIXES:	.o ${LIBEXT} .a .c .h

.PHONY:		all clean depend

############### make program
#
%_$(ARCH_MACH): %.o
	$(CXX) $(X_LIBS_ND) -o $@ $< $(LDLIBS)

%_$(ARCH_MACH)_debug: %.debug.o
	$(CXX) $(X_LIBS_D) -o $@ $< $(LDLIBS)

clean:
	-rm -f *~ *.o  config.status config.cache config.log

distclean: clean
	-rm -f $(PROGS) $(PROGS_DEBUG)

unix2dos:
	@printf "### Transforming all Visual C++ project and workspace files to Windows format... "
	@for i in `$(FIND) $(BASE_DIR) -type f \(\
	-name "*.dsp" -o            \
	-name "*.dsw"               \) -print`; \
	do $(PERL) -pi -e 's/\n/\r\n/;' $$i; done
	@echo "Done."








