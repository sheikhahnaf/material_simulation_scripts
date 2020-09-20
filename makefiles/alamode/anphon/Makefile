.SUFFIXES: .h .cpp
#-----------------------------------------------
# Makefile for the program 'ANPHON'.
# Please modify the variables properly.
# We recommend to use Intel c++ compiler.
#-----------------------------------------------

# Use gcc to activate openmp parallelization
MPICXX = mpiCC
CXXFLAGS = -O2 -std=c++11 -fopenmp
INCLUDE = -I../include -I$(HOME)/include -I/usr/lib/x86_64-linux-gnu/openmpi/include

LINKER = ${MPICXX}

## Mavericks
#LDFLAGS = -framework vecLib
## Yosemite
#LDFLAGS = -framework Accelerate

LAPACK =  -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi -lgomp -llapack -lblas
LIBS = ${LAPACK}  -lfftw3

#----------------------------------------------
# General rules
#----------------------------------------------

.cpp.o:
	${MPICXX} ${CXXFLAGS} ${INCLUDE} -c $<

PROG = anphon

CXXSRC = $(wildcard *.cpp)

OBJS= ${CXXSRC:.cpp=.o}

default: anphon

all: ${PROG}

anphon: ${OBJS}
	${LINKER} ${LDFLAGS} -o $@ ${OBJS} ${LIBS}

clean:
	rm -f ${OBJS}

.PHONY: clean

depend:
	gcc -MM *.cpp > .depend

-include .depend
