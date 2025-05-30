CC=g++
DEBUG=-g
FLAGS=-Wall -O2
LIBS=-lm
OBJS=mismatches.o utils.o string.o


all: run_kernel



utils.o: utils.h utils.cc config.h
	$(CC) $(DEBUG) $(FLAGS) -c utils.cc

string.o: string.h string.cc
	$(CC) $(DEBUG) $(FLAGS) -c string.cc

simgraph.o: simgraph.h simgraph.cc string.h
	$(CC) $(DEBUG) $(FLAGS) -c simgraph.cc

mismatches.o: mismatches.h mismatches.cc utils.h
	$(CC) $(DEBUG) $(FLAGS) -c mismatches.cc



gkernel.o: gkernel.h gkernel.cc config.h
	$(CC) $(DEBUG) $(FLAGS) -c gkernel.cc

run_kernel: run_kernel.cc gkernel.o simgraph.o $(OBJS) 
	$(CC) $(DEBUG) $(FLAGS) $(LIBS) run_kernel.cc gkernel.o simgraph.o $(OBJS) -o $@



clean:
	rm -rf run_kernel *.o *.*~ *~ core *.dSYM

