# Edit CFLAGS to suit your compiler;  
# make install will copy executables to BINDIR.

SHELL	= /bin/sh

CC	= gcc #-mpentiumpro -Wall  # insert your favorite C compiler
CFLAGS  = -g

LE_LIB =  -lsimfold -lLEModel 
#LE_INCLUDE = -I../LE  -I../LE/H   -I./
LE_INCLUDE = -I../LE  -I../LE/H   -I../simfold/include -I../simfold/src/simfold  -I../simfold/src/common  -I./
LE_PATH = -L../LE -L../simfold

X11_LIB = -lX11
X11_PATH = -L/usr/lib -L/usr/X11R6/lib 
X11_INCLUDE = -I/usr/X11R6/include
 

LDFLAGS = $(CFLAGS)  $(X11_PATH) $(LE_PATH) $(LE_LIB) $(X11_LIB)  -lm 
PROG = HotKnots

FILES   = HotKnot.c\
          HotKnotEnergy.c\
          hotspot.c\
          goodStem.c\
          utils.c\
          sc.cpp\
          score.cpp
          

OBJ     = HotKnot.o\
          HotKnotEnergy.o\
          hotspot.o\
          goodStem.o\
          sc.o\
          utils.o\
          score.o
          

default: $(PROG) # computeEnergy

$(PROG): $(OBJ)
	g++ -o $@ $(OBJ) $(LDFLAGS)
.SUFFIXES: .cpp .c
.cpp.o:
	g++ -c $(X11_INCLUDE) $(LE_INCLUDE) $*.cpp $(CFLAGS)
.c.o:
	g++ -c $(X11_INCLUDE) $(LE_INCLUDE) $*.c $(CFLAGS)

computeEnergy: computeEnergy.o
	g++ -o computeEnergy computeEnergy.o $(LDFLAGS)

clean:
	rm -f *.o $(PROG) computeEnergy






