# Author: Michael Casey
# (C) MIT Media Labs, January 1994, All Rights Reserved

INCDIR = /usr/tmp
LIBDIR = /usr/local/lib

LIBS = -laudio -lmidi /usr/lib/libcomplib.sgimath.a -lftn -lm

OBJS = main.o hrtf.o procbuf.o astream.o mstream.o libstuff.o arraylib.o

.SUFFIXES: .cc


CC = CC
cc = cc

CFLAGS = -g
APPNAME = 3Daudio

all: ${APPNAME}

${APPNAME}: ${OBJS} 
	${CC} ${CFLAGS} -o ${APPNAME} ${OBJS} -L${LIBDIR} ${LIBS}

.c.o:
	${cc} ${CFLAGS} -I${INCDIR} -c $<

.cc.o:
	${CC} ${CFLAGS} -I${INCDIR} -c $<
	
depends:
	gcc -I${INCDIR} -MM *.cc *.c > depends.make

install: ${APPNAME}
	cp ${APPNAME} ${BINDIR}
	
clean:
	rm -f *.o ${APPNAME}

include depends.make


