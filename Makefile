LFLAGS=-llapack -lcblas 
CFLAGS=-O3  -I/opt/intel/mkl/include/ -I/usr/include/i386-linux-gnu/ # -pg -g2 # -L/opt/intel/mkl/lib/32
DEBUG=-Wall -g
PAT=/home/karleik/masterProgging/vmc
NPROC=-n 2
CC=mpic++
EXEC=mpirun
GPROFOUT=gprof_testres.txt

vmcmain.out: all
	$(CC) $(CFLAGS) vmcmain.o sampler.o walker.o slaterMatrix.o orbital.o ipdist.o \
		zignor.o zigrandom.o newmatrix.o -o runVMC.out $(LFLAGS)

all: vmcmain.o sampler.o walker.o slaterMatrix.o orbital.o ipdist.o zignor.o zigrandom.o newmatrix.o

run: vmcmain.out
	$(EXEC) $(NPROC) runVMC.out

profile: all
	$(CC) $(CFLAGS) vmcmain.o sampler.o walker.o slaterMatrix.o orbital.o ipdist.o \
		zignor.o zigrandom.o newmatrix.o -o vmcmain_gprof.out $(LFLAGS)
	$(EXEC) $(NPROC) vmcmain_gprof.out 
	gprof ./vmcmain_gprof.out gmon.out > $(GPROFOUT)
	less $(GPROFOUT)

vmcmain.o: $(PAT)/vmcmain/vmcmain.cpp $(PAT)/definitions/vmcmain_Def.h
	$(CC) $(CFLAGS) -c $(PAT)/vmcmain/vmcmain.cpp 

sampler.o: $(PAT)/sampler/sampler.h $(PAT)/sampler/sampler.cpp walker.o $(PAT)/definitions/sampler_Def.h
	$(CC) $(CFLAGS) -c $(PAT)/sampler/sampler.h $(PAT)/sampler/sampler.cpp

walker.o: $(PAT)/walker/walker.h $(PAT)/walker/walker.cpp $(PAT)/definitions/randomNumberGenerators.h  slaterMatrix.o orbital.o ipdist.o zignor.o zigrandom.o newmatrix.o
	$(CC) $(CFLAGS) -c $(PAT)/walker/walker.h $(PAT)/walker/walker.cpp 

slaterMatrix.o: $(PAT)/QDslater/slaterMatrix.h $(PAT)/QDslater/slaterMatrix.cpp orbital.o newmatrix.o
	$(CC) $(CFLAGS) -c $(PAT)/QDslater/slaterMatrix.h $(PAT)/QDslater/slaterMatrix.cpp 

orbital.o: $(PAT)/orbital/orbital.h $(PAT)/orbital/orbital.cpp
	$(CC) $(CFLAGS) -c $(PAT)/orbital/orbital.h $(PAT)/orbital/orbital.cpp 

ipdist.o: $(PAT)/ipdist/ipdist.h $(PAT)/ipdist/ipdist.cpp newmatrix.o
	$(CC) $(CFLAGS) -c $(PAT)/ipdist/ipdist.h $(PAT)/ipdist/ipdist.cpp 

zignor.o: $(PAT)/ziggurat/zignor.h $(PAT)/ziggurat/zignor.c
	$(CC) $(CFLAGS) -c $(PAT)/ziggurat/zignor.h $(PAT)/ziggurat/zignor.c

zigrandom.o: $(PAT)/ziggurat/zigrandom.h $(PAT)/ziggurat/zigrandom.c
	$(CC) $(CFLAGS) -c $(PAT)/ziggurat/zigrandom.h $(PAT)/ziggurat/zigrandom.c

newmatrix.o: $(PAT)/newmatrix/newmatrix.h $(PAT)/newmatrix/newmatrix.cpp
	$(CC) $(CFLAGS) -c $(PAT)/newmatrix/newmatrix.h $(PAT)/newmatrix/newmatrix.cpp

clear:
	rm *.out *.o
