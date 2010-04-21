ARCH_MACH = $(shell uname -m)
OS_TYPE   = $(shell uname -s)

############### name of target
#
TARGET		= ./Approximation


############### determine machine type
#


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
SHARK = /home//ysong/dlim/shark
APPROMO = /home//ysong/dlim/ApproMo

############### libraries needed to generate program
#
LDLIBS		:=  	-lEALib_${ARCH_MACH} 		\
			-lRng_${ARCH_MACH} 		\
			-lArray_${ARCH_MACH} 		\
			-lLinAlg_${ARCH_MACH} 		\
			-lReClaM_${ARCH_MACH}		

endif

############### objects needed to generate target
#
#Main_RTypeTrim
#Main_FrontWingTrim
OBJS = 	ApproximationTest.o				\
	../src/Database.o				




PURIFY	= purify $(PFLAGS)


CCFLAGS		=  $(CPPFLAGS)
COMPILE.cc	=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  -c 
LINK.cc		=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  $(LDFLAGS)



# for Solaris

############### options for preprocessor (include directories)
#
CPPFLAGS	:=\
		$(CPPFLAGS)  			\
		-Wall                           \
		-I../include			\
		-I/home//ysong/dlim/shark/include/Array	\
		-I/home//ysong/dlim/shark/NeuralNet/include\
		-I/home//ysong/dlim/shark/Database/include \
		-I/home//ysong/dlim/shark/include/ReClaM 	\
		-I/home//ysong/dlim/shark/include/EALib	




############### options for compiler
#

CCFLAGS:= -g -Wall -pedantic
# -g -pg 

CXXFLAGS:=$(CXXFLAGS) $(CFLAGS)

############### options for linker
#

############### options for linker
#

LDFLAGS		:= -L/home//ysong/dlim/shark/lib_i686/ 	\
		   $(LDFLAGS)	

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
		-$(RM) Makefile.bak

###############	update dependencies
#
depend:
		makedepend -Y -- $(CPPFLAGS) *.c *.cc *.cpp 2> /dev/null

cleandepend:
		makedepend

# DO NOT DELETE
