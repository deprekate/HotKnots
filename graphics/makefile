SHELL	= /bin/sh

CC	= g++ #-mpentiumpro -Wall  # insert your favorite C compiler
BINDIR	= ../bin
MANDIR	= ../man/man1

COPT    = -O  -g
LIB = -lX11
LIB_PATH = -L/usr/lib -L/usr/X11R6/lib
X11_INCLUDE = -I/usr/X11R6/include

all:	Plot PLOTbpseq.o graphics.o PlotRna.o

graphics.o: graphics.c
	$(CC) -c -g  $(X11_INCLUDE) $(LIB_PATH) graphics.c

PlotRna.o: PlotRna.c
	$(CC) -c -g  PlotRna.c

Plot: PLOTbpseq.o PlotRna.o graphics.o
	$(CC) -o Plot PLOTbpseq.o PlotRna.o graphics.o $(LIB_PATH) $(LIB)

PLOTbpseq.o: PLOTbpseq.c graphics.h
	$(CC) -c -g  $(X11_INCLUDE)  PLOTbpseq.c

clean:
	rm -f *.o
