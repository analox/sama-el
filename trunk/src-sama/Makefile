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
CCC		= ../gcc32/bin/g++
CC		= ../gcc32/bin/gcc

############### Other Packages
#
SHARKDIR = ../shark

############### libraries needed to generate program
#
LDLIBS		:=  	-lArray_${ARCH_MACH} 		\
			-lEALib_${ARCH_MACH} 		\
			-lRng_${ARCH_MACH}		\
			-lLinAlg_${ARCH_MACH}		\
			-lfifo				\
			-llbfgsb			\
			-lInterpolateRBF		\
			-lRegression			\
			-lDatabase			\
			-lKriging                       \
                        -lLAPACK                        \
                        -lBLAS                          \
                        -lg2c  
############### objects needed to generate target
#
#OBJS = 	sama.o tools/distribution/WeightCal.o tools/distribution/UniformDistFunc.o tools/distribution/GaussianDistFunc.o \
	tools/database/DBEntry.o tools/database/ILDatabase.o tools/utilities.o OSS2/oss2.o

OBJS = sama-airfoil.o tools/distribution/WeightCal.o tools/distribution/UniformDistFunc.o tools/distribution/GaussianDistFunc.o \
	tools/database/DBEntry.o tools/database/ILDatabase.o tools/utilities.o OSS2/oss2.o

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
		-I../gcc32/lib \
		-I.	\
		$(CPPFLAGS)


############### options for compiler
#

CCFLAGS:= -g -Wall -pedantic -O3

CXXFLAGS:=$(CXXFLAGS) $(CFLAGS)


############### options for linker
#

LDFLAGS		:=  -L${SHARKDIR}/lib_${ARCH_MACH}	\
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

