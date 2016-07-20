# Makefile

RM = rm
binfile = ssmst
objs = tree.o mstree.o

CXX = g++
AR = ar

CFLAGS = -O3 -fPIC -Wall -L./ -I./

.SUFFIXES: .cpp .o

$(binfile): $(objs)
	$(CXX) -o $(binfile) $(objs) -L. -lslidesort -fopenmp -lmslidesort

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) $(objs) *~

tree.o: olmst.h
mstree.o: olmst.h

