LFLAGS=-llapack -lcblas
CFLAGS=-O4
DEBUG=-Wall -g
PAT=/home/karleik/masterProgging/vmc
CC=mpic++

vmcsolver.out: vmcsolver.o slaterMatrix.o orbital.o ipdist.o zignor.o zigrandom.o
	$(CC) $(CFLAGS) vmcsolver.o slaterMatrix.o orbital.o ipdist.o zignor.o zigrandom.o -o vmcsolver.out $(LFLAGS)

all: vmcsolver.o slaterMatrix.o orbital.o ipdist.o zignor.o zigrandom.o

run: vmcsolver.out
	$(CC) $(CFLAGS) vmcsolver.o slaterMatrix.o orbital.o ipdist.o zignor.o zigrandom.o -o vmcsolver.out $(LFLAGS)
	mpiexec -n 2 vmcsolver.out

profile: all
	c++ -O4  -pg -g2 vmcsolver.o slaterMatrix.o orbital.o ipdist.o zignor.o zigrandom.o -o vmcsolver.out -lcblas -llapack
	./vmcsolver.out
	gprof ./vmcsolver.out gmon.out > testres.txt
	less testres.txt

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
