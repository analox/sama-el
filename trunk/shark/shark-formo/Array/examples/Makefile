# Array examples makefile

############### determine the host type
#
HOST			= Linux
ARCH_MACH		= i686

############### variables set by configure script
#
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
libdir			= $(basedir)/lib


############### additional include directories
#

############### additional settings
#
LDLIBS         = -lArray_$(ARCH_MACH)
LDFLAGS        := -L$(sharkdir)/Array/lib
LDFLAGS_DEBUG  := -L$(sharkdir)/Array/lib/debug

############## os specific settings
#
CFLAGS		   += 	-I$(top_inc_dir) $(DEFS)
CXXFLAGS	    = 	$(CFLAGS)

############### example progs
#
PROGS			= mean-var_$(ARCH_MACH)
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
	$(CXX) $(LDFLAGS) -o $@ $< $(LDLIBS)

%_$(ARCH_MACH)_debug: %.debug.o
	$(CXX) $(LDFLAGS_DEBUG) -o $@ $< $(LDLIBS)

clean:
	-rm -f *.o *~ raw red white pca

distclean: clean
	-rm -f $(PROGS) $(PROGS_DEBUG)

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









