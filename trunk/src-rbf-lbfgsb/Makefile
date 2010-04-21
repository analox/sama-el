ARCH_MACH = $(shell uname -m)

############### name of target
#
TARGET		= ./sama

########################################################################
##########        Linux Cluster Environment                   ##########
########################################################################

CPPFLAGS := $(CPPFLAGS) -D__LINUX__

############### C++ compiler commands
#
CCC		= /opt/gcc32/bin/g++
CC		= /opt/gcc32/bin/gcc

############### Other Packages
#
SHARKDIR = /home/nghia/dlim/shark

############### libraries needed to generate program
#
LDLIBS		:=  	-lArray_${ARCH_MACH} 		\
			-lEALib_${ARCH_MACH} 		\
			-lRng_${ARCH_MACH}		\
			-lLinAlg_${ARCH_MACH}		\
			-lRegression			\
			-lDatabase			\
			-lfifo				\
			-llbfgsb			\
			-lInterpolateRBF		\
			-lKriging                       \
                        -lLAPACK                        \
                        -lBLAS                          \
                        -lg2c  

############### objects needed to generate target
#
OBJS = 	sama.o			

PURIFY	= purify $(PFLAGS)

CCFLAGS		=  $(CPPFLAGS)
COMPILE.cc	=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  -c -w
LINK.cc		=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  $(LDFLAGS)



############### options for preprocessor (include directories)
#
CPPFLAGS	:= \
		-Wall 				\
		-I${SHARKDIR}/include		\
		-I${SHARKDIR}/include/Array	\
		-I${SHARKDIR}/include/EALib	\
		-I./include		\
		$(CPPFLAGS)


############### options for compiler
#

CCFLAGS:= -g -Wall -pedantic

CXXFLAGS:=$(CXXFLAGS) $(CFLAGS)


############### options for linker
#

LDFLAGS		:=  -L${SHARKDIR}/lib_i686	\
		    -L./lib			\
		    $(LDFLAGS)

###############	make targets
#
all:		$(TARGET) $(INCLUDE)

.c:
		$(LINK.c) -o $@ $< $(LDLIBS)
.c.o:
		$(COMPILE.c) -o $@ $<
.cpp:
		$(LINK.cc) -o $@ $< $(LDLIBS)
.cpp.o:
		$(COMPILE.cc) -o $@ $<
		
############### suffixes list
#
.SUFFIXES:	.o .so .a .c .cpp .h

.PHONY:		all clean

############### make shared object (dynamic link library)
#
$(TARGET):	$(OBJS)
		$(LINK.cc) -o $@ $(OBJS) $(LDLIBS)
		-$(RM) $(OBJS)

###############	clean up
#
clean:
		-$(RM) Makefile.bak sama $(OBJS) *~


############### update dependencies
#
depend:
	makedepend -Y -- $(CPPFLAGS) *.c *.cc *.cpp 2> /dev/null

cleandepend:
	makedepend

