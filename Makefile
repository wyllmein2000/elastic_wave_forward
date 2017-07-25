# Date: 15 May 2014
#
MPIROOT = /usr/local/mpich-install
#MPIROOT = /home/wuyx0a/mpich3
INC = -I$(MPIROOT)/include
LIB = -L$(MPIROOT)/lib
BIN = ../../bin

FC = gfortran
#CC = g++
#CC = mpic++
CC = $(MPIROOT)/bin/mpicxx
OPT = -c -msse2 -O3

SUB = ../acous2d_new/sub_tools.f90

OBJSUB = shotdomain.o wavelet.o getpara.o numpy.o array.o inout.o
OBJACOU2D = main_acous.o  Afd2d.o $(OBJSUB)
OBJELAS2D = main_forward.o fd2d.o $(OBJSUB)

#ALL = elas3danaly beachball acou2d elas2d
ALL = beachball acou2d elas2d

all: $(ALL)

elas3danaly: main_analyt.f90
	$(FC) -o $(BIN)/elas3danaly main_analyt.f90 $(SUB)

beachball: beachball.f90
	$(FC) -o $(BIN)/beachball beachball.f90

elas2d: $(OBJELAS2D)
	$(CC) -fopenmp -o $(BIN)/elas2d $(OBJELAS2D) $(LIB)

acou2d: $(OBJACOU2D)
	$(CC) -fopenmp -o $(BIN)/acou2d $(OBJACOU2D) $(LIB)

.cpp.o:
	$(CC) $(OPT) $(INC) -fopenmp $<

clean:
	rm -f *.o
