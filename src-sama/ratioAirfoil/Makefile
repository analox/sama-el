# GNU make version 3.79.1
# Makefile for clcdgrad

OBJSEX= \
	ratiograd.o \
	example.o
SRCSEX= \
	ratiograd.f \
	example.c


OBJS= ratiograd.o ratioairfoilObj.o
SRCS= ratiograd.f ratioairfoilObj.c

OBJSIN= ratiogradin.o ratioairfoilObjIn.o
SRCSIN= ratiogradin.f ratioairfoilObjIn.c

OBJSNS= NStest.o
SRCSNS= NStest.c

OBJS_ShapePress= ratiograd.o ShapePress.o
SRCS_ShapePress= ratiograd.f ShapePress.c

CC=gcc

CFLAGS=
LDLIBS= -lg2c -lm
LDFLAGS= 
RM = /bin/rm -f

.c.o:
	$(CC) -c $(CFLAGS) $<

.f.o:
	$(CC) -c $(CFLAGS) $<

all: example

clean:
	$(RM) *.o *~ core

example: $(OBJSEX)
	$(CC) -o $@ $(LDFLAGS) $(OBJSEX) $(LDLIBS)

ratioairfoilObjIn: $(OBJSIN)
	$(CC) -o $@ $(LDFLAGS) $(OBJSIN) $(LDLIBS)

ratioairfoilObj: $(OBJS)
	$(CC) -o $@ $(LDFLAGS) $(OBJS) $(LDLIBS)

NSTest: $(OBJSNS)
	$(CC) -o $@ $(LDFLAGS) $(OBJSOBJSNS) /staff/asysong/Netsolve/NetSolve-1.4.1/lib/i686_pc_linux_gnu/libnetsolve.a $(LDLIBS)

ShapePress: $(OBJS_ShapePress)
	$(CC) -o $@ $(LDFLAGS) $(OBJS_ShapePress) $(LDLIBS)
	