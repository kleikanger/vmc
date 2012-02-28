LFLAGS=-llapack -lcblas
CFLAGS=-O4
DEBUG=-Wall -g
PAT=/home/karleik/masterProgging/vmc
NPROC=-n 2
CC=mpic++
GPROFOUT=gprof_testres.txt

vmcsolver.out: all
	$(CC) $(CFLAGS) vmcsolver.o slaterMatrix.o orbital.o ipdist.o zignor.o zigrandom.o -o vmcsolver.out $(LFLAGS)

all: vmcsolver.o slaterMatrix.o orbital.o ipdist.o zignor.o zigrandom.o

run: vmcsolver.out
	mpiexec $(NPROC) vmcsolver.out

profile: all
	$(CC) $(CFLAGS) -pg -g2 vmcsolver.o slaterMatrix.o orbital.o ipdist.o zignor.o zigrandom.o -o vmcsolver_gprof.out $(LFLAGS)
	mpiexec $(NPROC) vmcsolver_gprof.out 
	gprof ./vmcsolver_gprof.out gmon.out > $(GPROFOUT)
	less $(GPROFOUT)

vmcsolver.o: vmcsolver.h vmcsolver.cpp
	$(CC) $(CFLAGS) -c vmcsolver.h vmcsolver.cpp $(LFLAGS)

slaterMatrix.o: $(PAT)/QDslater/slaterMatrix.h $(PAT)/QDslater/slaterMatrix.cpp orbital.o
	$(CC) $(CFLAGS) -c $(PAT)/QDslater/slaterMatrix.h $(PAT)/QDslater/slaterMatrix.cpp $(LFLAGS)

orbital.o: $(PAT)/orbital/orbital.h $(PAT)/orbital/orbital.cpp
	$(CC) $(CFLAGS) -c $(PAT)/orbital/orbital.h $(PAT)/orbital/orbital.cpp $(LFLAGS)

ipdist.o: $(PAT)/ipdist/ipdist.h $(PAT)/ipdist/ipdist.cpp
	$(CC) $(CFLAGS) -c $(PAT)/ipdist/ipdist.h $(PAT)/ipdist/ipdist.cpp $(LFLAGS)

zignor.o: $(PAT)/ziggurat/zignor.h $(PAT)/ziggurat/zignor.c
	$(CC) $(CFLAGS) -c $(PAT)/ziggurat/zignor.h $(PAT)/ziggurat/zignor.c

zigrandom.o: $(PAT)/ziggurat/zigrandom.h $(PAT)/ziggurat/zigrandom.c
	$(CC) $(CFLAGS) -c $(PAT)/ziggurat/zigrandom.h $(PAT)/ziggurat/zigrandom.c

clear:
	rm *.out *.o
