ARCH_MACH = $(shell uname -m)
OS_TYPE   = $(shell uname -s)

############### name of target
#
TARGET		= ./FIFOTest


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
SHARKDIR = /home/dlim/dlim/shark
APPROMO = /home/dlim/dlim/ApproMo


############### libraries needed to generate program
#
LDLIBS		:=  	-lArray_${ARCH_MACH} 		

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
SHARKDIR = /home/dlim/dlim/shark
APPROMO = /home/dlim/dlim/ApproMo

############### libraries needed to generate program
#
LDLIBS		:=  	-lArray_${ARCH_MACH} 		
			
endif

############### objects needed to generate target
#
OBJS =  ../src/fifo.o	FIFOTest.o			

PURIFY	= purify $(PFLAGS)

CCFLAGS		=  $(CPPFLAGS)
COMPILE.cc	=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  -c -w
LINK.cc		=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  $(LDFLAGS)


# for Solaris

############### options for preprocessor (include directories)
#
CPPFLAGS	:=\
		$(CPPFLAGS)  			\
		-Wall 		\
		-I../include                          \
		-I${SHARKDIR}/include		


############### options for compiler
#

CCFLAGS:= -g -Wall -pedantic

CXXFLAGS:=$(CXXFLAGS) $(CFLAGS)


############### options for linker
#

LDFLAGS		:= $(LDFLAGS) 			\
   		    -L${SHARKDIR}/lib_i686	
		
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
		-$(RM) ../src/fifo.o Makefile.bak FIFOTest *~
#		-$(RM) $(OBJS)


###############	update dependencies
#
depend:
		makedepend -Y -- $(CPPFLAGS) *.c *.cc *.cpp 2> /dev/null

cleandepend:
		makedepend

# DO NOT DELETE
