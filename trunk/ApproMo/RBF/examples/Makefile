ARCH_MACH = $(shell uname -m)
OS_TYPE   = $(shell uname -s)

############### name of target
#
TARGET		= ./RBFTest


############### determine machine type
#

########################################################################
##########        Sun Solaris Path                            ##########
########################################################################
ifeq "${ARCH_MACH}" "sun4u"

CPPFLAGS := $(CPPFLAGS) -D__SOLARIS__

############### C++ compiler commands
#
CC		= purify /usr/bin/gcc
CCC		= purify /usr/bin/g++

############### Other Packages
#
#SHARK  = $(SHARKDIR) 
SHARK = /home/ysong/dlim/shark
#APPROMO = $(APPROMO)
APPROMO = /home/ysong/dlim/ApproMo/

EVONET = ../../../EvoNet

############### libraries needed to generate program
#
LDLIBS		:=  	-lEALib_${ARCH_MACH} 		\
			-lRng_${ARCH_MACH}		\
			-lArray_${ARCH_MACH} 		\
			-lLinAlg_sun4u			\
			-lm				\
			-lX11 -lXt -lXext

endif

########################################################################
##########        Linux Cluster Environment                   ##########
########################################################################

ifeq "${ARCH_MACH}" "i686"
CPPFLAGS := $(CPPFLAGS) -D__LINUX__ -DEBUG

############### C++ compiler commands
#
CCC		= /usr/bin/g++ 
CC		= /usr/bin/gcc

############### Other Packages
#
SHARK = $(SHARKDIR)
SISL  = $(SISLDIR)
EVONET = $(EVONETDIR)
EVONET20 = $(EVONET2)

############### libraries needed to generate program
#
LDLIBS		:=  	-lRng_${ARCH_MACH} 		\
			-lArray_${ARCH_MACH} 		\
			-lLinAlg_${ARCH_MACH} 		\
			-lReClaM_${ARCH_MACH}		\
			-lEALib_${ARCH_MACH} 		\
			-lm	

endif

############### objects needed to generate target
#
#Main_RTypeTrim
#Main_FrontWingTrim
OBJS = 	RBFTest.o				\
	../../RBF/src/rbf.o			\


############### PVM specifications ##############
#
# CPPVM_ROOT = /opt/cppvm
# PVMLIB 		= 	-lpvm3 -lgpvm3 -lnsl 
# PVMINCLUDE 	= -I$(PVM_ROOT)/include -I$(CPPVM_ROOT)/h
# PVMLDPATH  	= -L$(PVM_ROOT)/lib/$(PVM_ARCH) 

############### FORTRAN compiler commands
#
# CCF		= g77
# CCFLAGS		=  $(FFLAGS)
# COMPILE.f	= $(CCF) $(CCFLAGS) $(CPPFLAGS) -c -w
# LINK.f		= $(CCF) $(CCFLAGS) $(CPPFLAGS) $(LDFLAGS)

############### Purify related flags
# PFLAGS	= -messages=first -leaks-at-exit=yes

PURIFY	= purify $(PFLAGS)


CCFLAGS		=  $(CPPFLAGS)
COMPILE.cc	=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  -c -w
LINK.cc		=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  $(LDFLAGS)



# for Solaris

############### options for preprocessor (include directories)
#
CPPFLAGS	:=\
		$(CPPFLAGS)  			\
		-Wall                           \
		-I/home/ysong/dlim/shark/include/Array	\
		-I/home/ysong/dlim/shark/LinAlg/include	\
		-I/home/ysong/dlim/shark/include/EALib	\
		-I/home/ysong/dlim/shark/include/Rng		\
		-I/home/ysong/dlim/shark/include/		\
		-I../include			

#-I${PVM_ROOT}/include 		\



############### options for compiler
#

CCFLAGS:= -g -Wall -pedantic
# -g -pg 

CXXFLAGS:=$(CXXFLAGS) $(CFLAGS)

############### options for linker
#

############### options for linker
#

LDFLAGS		:= $(LDFLAGS) 				\
		   -L/home/dlim/shark/lib_i686/ 		\
	           -L../lib				
#		   -L$(TCLTK)/				\
#		    $(PVMLDPATH)                        \
#		   -L$(SISL)/lib 			\
#	  	   -L$(CPPVM_ROOT)/lib/LINUX/		\
#		   -L/usr/openwin/lib -L/usr/dt/lib	\
#		   -L./					\
#		   -L../Kriging/lib

###############	make targets
#
all:		$(TARGET) $(INCLUDE)


###############	do not modify the following lines ###############

###############	implicit rules
#


.c:
		$(LINK.c) -o $@ $< $(LDLIBS)
.c.o:
		$(COMPILE.c) -o $@ $<
.cpp:
		$(LINK.cc) -o $@ $< $(LDLIBS)
.cpp.o:
		$(COMPILE.cc) -o $@ $<
.f:
		$(LINK.f) -o $@ $< $(LDLIBS)
.f.o:
		$(COMPILE.f) -o $@ $<
############### suffixes list
#
.SUFFIXES:	.o .so .a .c .cpp .h

.PHONY:		all clean depend

############### make shared object (dynamic link library)
#
$(TARGET):	$(OBJS)
		$(LINK.cc) -o $@ $(OBJS) $(LDLIBS)

###############	clean up
#
clean:
		-$(RM) $(OBJS)
		-$(RM) Makefile.bak RBFTest *~

###############	update dependencies
#
depend:
		makedepend -Y -- $(CPPFLAGS) *.c *.cc *.cpp 2> /dev/null

cleandepend:
		makedepend

# DO NOT DELETE
