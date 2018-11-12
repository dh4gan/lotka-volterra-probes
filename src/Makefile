#####################################################
###                                               ###
### Makefile for lotka-volterra-probes            ###
###                                               ###
### Duncan H. Forgan 12/11/2018          	  ###
###       				          ###
###                                               ###
#####################################################

# Compiler variables:
CC     = g++

VERSION = $(shell git describe --tags --abbrev=0)

# Git version
GIT_HASH = $(shell git describe --abbrev=9 --dirty --always)

# Define compiler flags
CFLAGS = -O3 -DVERSION=\"$(VERSION)\" -DGIT_HASH=\"${GIT_HASH}\" -std=c++11

# Create object files:
%.o: %.cpp
	$(CC) $(CFLAGS) -c $<
%.o: %.f90
	$(CC) $(CFLAGS) -c $<

# Source files (.F90)
SOURCESA = main.cpp Vector3D.cpp parFile.cpp Body.cpp Edge.cpp Vertex.cpp LKVertex.cpp Graph.cpp
OBJECTSA    = $(SOURCESA:.cpp=.o)

# Create executable files:
build: lotka_volterra_probes

lotka_volterra_probes:  $(OBJECTSA)
	$(CC) $(CFLAGS) -o $@ $(OBJECTSA)

# Clean statements:
clean: 
	\rm *.o lotka_volterra_probes

# End Makefile
