ARCH_MACH = $(shell uname -m)
OS_TYPE   = $(shell uname -s)

############### name of target
#
TARGET		= ./Neural3DScan


############### determine machine type
#

########################################################################
##########        Sun Solaris Path                            ##########
########################################################################
ifeq "${ARCH_MACH}" "sun4u"

CPPFLAGS := $(CPPFLAGS) -D__SOLARIS__

############### C++ compiler commands
#
CC		= purify /usr/bin/gcc32
CCC		= purify /usr/bin/g++32

############### Other Packages
#
SHARK  = $(SHARKDIR) 
SHARK = /home/dlim/shark32
SISL   = $(SISLDIR)
SISL = /opt2/SISL/
TCLTK = /usr/local/lib
APPROMO = /home/dlim/Software/EL-Tec/ApproMo
############### libraries needed to generate program
#
LDLIBS		:=  	-lEALib_${ARCH_MACH} 		\
			-lRng_${ARCH_MACH}		\
			-lArray_${ARCH_MACH} 		\
			-lLinAlg_sun4u			\
			-lsisl  			\
			-lReClaM_${ARCH_MACH} 		\
			-ltcl8.3				\
			-ltk8.3 				\
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
CCC		= /usr/bin/g++32 
CC		= /usr/bin/gcc32

############### Other Packages
#
SHARK = $(SHARKDIR)
SISL  = $(SISLDIR)

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
#Neural3DScan.o					\

OBJS =  Neural3DScan.o 	\
	/home/dlim/Software/EL-Tec/ApproMo/Database/src/Database.o			\
	/home/dlim/Software/EL-Tec/ApproMo/NeuralNet/src/EvoNet.o			\
	/home/dlim/Software/EL-Tec/ApproMo/NeuralNet/src/Ops_Basic.o			\
	/home/dlim/Software/EL-Tec/ApproMo/NeuralNet/src/FFOps_OL_neu.o 		\
	/home/dlim/Software/EL-Tec/ApproMo/NeuralNet/src/FFOps_Basic_neu.o		\
	/home/dlim/Software/EL-Tec/ApproMo/NeuralNet/src/FFManip.o


CCFLAGS		=  $(CPPFLAGS)
COMPILE.cc	=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  -c
LINK.cc		=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  $(LDFLAGS)



# for Solaris

############### options for preprocessor (include directories)
#
CPPFLAGS	:=\
		$(CPPFLAGS)  			\
		-Wall                           \
		-I/home/dlim/shark32/include/		\
		-I/home/dlim/shark32/include/Array	\
		-I/home/dlim/Software/EL-Tec/ApproMo/NeuralNet/include\
		-I/home/dlim/Software/EL-Tec/ApproMo/Database/include \
		-I/home/dlim/shark32/include/ReClaM 	\
		-I/home/dlim/shark32/include/EALib	


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
		   -L$(SHARK)/lib_i686/ 		
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
		-$(RM) Makefile.bak

###############	update dependencies
#
depend:
		makedepend -Y -- $(CPPFLAGS) *.c *.cc *.cpp 2> /dev/null

cleandepend:
		makedepend

# DO NOT DELETE


