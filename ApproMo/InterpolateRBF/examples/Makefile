ARCH_MACH = $(shell uname -m)

############### name of target
#
TARGET		= ./InterpolateRBFTest


CPPFLAGS := $(CPPFLAGS) -D__LINUX__ -DEBUG

############### C++ compiler commands
#
CCC		= /usr/bin/g++
CC		= /usr/bin/gcc

############### Other Packages
#
SHARKDIR = /home/dlim/dlim/shark

############### libraries needed to generate program
#
LDLIBS		:=  	-lRng_${ARCH_MACH} 		\
			-lArray_${ARCH_MACH} 		\
			-lLinAlg_${ARCH_MACH} 		\
			-lEALib_${ARCH_MACH} 		\
			-lm	


############### objects needed to generate target
#
OBJS = 	InterpolateRBFTest.o				\
	../../InterpolateRBF/src/InterpolateRBF.o	

############### Purify related flags
PURIFY	= purify $(PFLAGS)


CCFLAGS		=  $(CPPFLAGS)
COMPILE.cc	=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  -c -w
LINK.cc		=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  $(LDFLAGS)




############### options for preprocessor (include directories)
#
CPPFLAGS	:=\
		$(CPPFLAGS)  			\
		-pedantic                           \
		-I${SHARKDIR}/include/Array	\
		-I${SHARKDIR}/LinAlg/include	\
		-I${SHARKDIR}/include/EALib	\
		-I${SHARKDIR}/include/Rng	\
		-I${SHARKDIR}/include/		\
		-I../include			



############### options for compiler
#
CCFLAGS:= -g -pedantic

CXXFLAGS:=$(CXXFLAGS) $(CFLAGS)

############### options for linker
#

LDFLAGS		:= $(LDFLAGS) 				\
		   -L${SHARKDIR}/lib_i686/ 	

###############	make targets
#
all:		$(TARGET) $(INCLUDE)


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

############### suffixes list
#
.SUFFIXES:	.o .so .a .c .cpp .h

.PHONY:		all clean

############### make shared object (dynamic link library)
#
$(TARGET):	$(OBJS)
		$(LINK.cc) -o $@ $(OBJS) $(LDLIBS)

###############	clean up
#
clean:
		-$(RM) $(OBJS)
		-$(RM) Makefile.bak InterpolateRBFTest *~


############### update dependencies
#
depend:		
	makedepend -Y -- $(CPPFLAGS) *.c *.cc *.cpp 2> /dev/null

cleandepend:	
	makedepend

