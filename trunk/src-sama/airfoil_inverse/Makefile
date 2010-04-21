# GNU make version 3.79.1
# Makefile for clcdgrad

OBJS= clcdgrad.o airfoilObj.o
SRCS= clcdgrad.f airfoilObj.cpp

OBJS_Ex= pdgrad.o example.o
SRCS_Ex= pdgrad.f example.c

OBJS_PRESSUREIN= pdgradin.o airfoilPressureProIn.o
SRCS_PRESSUREIN= pdgradin.f airfoilPressureProIn.c

OBJS_PRESSUREIN2= airfoilPressureProIn2.o NSPressIN.a
SRCS_PRESSUREIN2= airfoilPressureProIn2.c

OBJS_lana= airfoilana.o lana.o pdgradin.o
SRCS_lana= airfoilana.c lana.c pdgradin.f

OBJS_ShapePressIN= pdgradin.o ShapePress.o
SRCS_ShapePressIN= pdgradin.f ShapePress.c

CC=gcc

CFLAGS=
LDLIBS= -lg2c -lm
LDFLAGS=
RM = /bin/rm -f

.c.o:
	$(CC) -c $(CFLAGS) $<

.f.o:
	$(CC) -c $(CFLAGS) $<

all: airfoilObj

clean:
	$(RM) *.o *~ core

example: $(OBJS_Ex)
	$(CC) -o $@ $(LDFLAGS) $(OBJS_Ex) $(LDLIBS)

airfoilPressureIN: $(OBJS_PRESSUREIN)
	$(CC) -o $@ $(LDFLAGS) $(OBJS_PRESSUREIN) $(LDLIBS)

airfoilPressureIN2: $(OBJS_PRESSUREIN2)
	$(CC) -o $@ $(LDFLAGS) $(OBJS_PRESSUREIN2) $(LDLIBS)

lana: $(OBJS_lana)
	$(CC) -o $@ $(LDFLAGS) $(OBJS_lana) $(LDLIBS)

ShapePressIN: $(OBJS_ShapePressIN)
	$(CC) -o $@ $(LDFLAGS) $(OBJS_ShapePressIN) $(LDLIBS)

