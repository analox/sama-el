ARCH_MACH = $(shell uname -m)
OS_TYPE   = $(shell uname -s)

############### name of target
#
TARGET		= ./myDaceExam




########################################################################
##########        Linux Cluster Environment                   ##########
########################################################################

CPPFLAGS := $(CPPFLAGS) -D__LINUX__ -DEBUG

############### C++ compiler commands
#
CCC		= /usr/bin/g++32
CC		= /usr/bin/gcc32

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
			-lm	\
			-lLAPACK			\
			-lBLAS				\
			-lg2c	\
			-lDatabase

############### objects needed to generate target
#
#Main_RTypeTrim
#Main_FrontWingTrim
OBJS = 	myDaceExam.o					\
	../src/rngs.o					\
	../src/rvgs.o					\
	../src/krig.o					\
	../src/krigify.o 				\
	../src/gamma.o 					\
	../src/ParamEstimate.o 				\
	../src/DirectSearch.o				\
	../src/PatternSearch.o 				\
	../src/CompassSearch.o 				\
	../src/approx.o 				\
	../src/objective.o				\
	../src/KrigApprox.o



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
		-I/home/ysong/dlim/ApproMo/Database/include \
		-I/home/ysong/dlim/shark/include/Array	\
		-I/home/ysong/dlim/ApproMo/NeuralNet/include	\
		-I/home/ysong/dlim/ApproMo/Database/include 	\
		-I/home/ysong/dlim/shark/include/ReClaM 	\
		-I/home/ysong/dlim/shark/FileUtil/include	\
		-I/home/ysong/dlim/shark/Basics/include	\
		-I/home/ysong/dlim/shark/Math/include		\
		-I/home/ysong/dlim/shark/LinAlg/include	\
		-I/home/ysong/dlim/shark/include/EALib	\
		-I/home/ysong/dlim/shark/include/Rng		\
		-I/home/ysong/dlim/shark/include/		\
		-I../include			

#-I${PVM_ROOT}/include 		\

#		-I$(PAGODE_HOME)plot_widgets/plot_widgets \


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
		    -L/usr/lib		\
		    -L/home/ysong/dlim/shark/lib_i686/ 		\
		    -L/home/ysong/dlim/ApproMo/Database/lib	\
	            -L../lib				

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
		-$(RM) Makefile.bak myDaceExam *~

###############	update dependencies
#
depend:
		makedepend -Y -- $(CPPFLAGS) *.c *.cc *.cpp 2> /dev/null

cleandepend:
		makedepend

