.SUFFIXES: .h .cpp
#-----------------------------------------------
# Makefile for the program 'ALM'.
# Please modify the variables properly.
# We recommend to use Intel c++ compiler.
#-----------------------------------------------

CXX = icpc 
CXXFLAGS = -O2 -xHOST -qopenmp -std=c++11
INCLUDE = -I../include -I$(HOME)/include -I/usr/local/include

CXXL = ${CXX}
LDFLAGS = -mkl -L/usr/local/lib -lsymspg

LAPACK = 
LIBS = ${LAPACK}

#-----------------------------------------------
# General rules 
#-----------------------------------------------

.cpp.o:
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $<

PROG = alm

CXXSRC= $(wildcard *.cpp)

OBJS= ${CXXSRC:.cpp=.o}

default: alm

all: ${PROG}

alm: ${OBJS}
	${CXXL} ${LDFLAGS} -o $@ ${OBJS} ${LIBS}

clean:
	rm -f ${OBJS}

.PHONY: clean

depend:
	gcc -MM *.cpp > .depend

-include .depend
