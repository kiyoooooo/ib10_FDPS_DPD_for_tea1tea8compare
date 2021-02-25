#PS_PATH = -I../fdps
PS_PATH = -I../../src/
#PS_PATH = -I../../FDPS-master/src/

#CC = g++
#CC = icc
#CC = mpic++
CC = mpiicpc
#CC = /opt/local/bin/mpic++-mpich-mp
#CC = time mpicxx
#CC = time mpicc
CXXFLAGS = -std=c++11
#CXXFLAGS = -std=icc
CXXFLAGS += -O3
CXXFLAGS += -Wall
CXXFLAGS += -ffast-math
CXXFLAGS += -funroll-loops

#CXXFLAGS += -g -O0

CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL


#CXXFLAGS += -DNOSE_HOOVER
CXXFLAGS += -DDISSIPATIVE_RANDOM

#CXXFLAGS += -DSMALL_SOLVENT

#CXXFLAGS += -DNANOSLIT
#CXXFLAGS += -DNANOTUBE


SRC = janus.cpp
PROGRAM = $(SRC:%.cpp=%.out)

$(PROGRAM):$(SRC)
	$(CC) $(MULEXP) $(PS_PATH) $(CXXFLAGS) -o $@ $<

clean:
	rm $(PROGRAM) 


