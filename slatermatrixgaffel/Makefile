test: slaterMatrix.cpp slaterMatrix.h ../orbital/orbital.cpp ../orbital/orbital.h ../lib/lib.cpp test.cpp ../lib/lib.h 
	c++ -o test.out -Wall test.cpp slaterMatrix.cpp slaterMatrix.h ../orbital/orbital.h ../orbital/orbital.cpp ../lib/lib.h ../lib/lib.cpp  -llapack -lcblas

run: test
	make test && ./test.out

profile: test
	g++ -pg -g2 -o test.out test.cpp slaterMatrix.cpp slaterMatrix.h ../orbital/orbital.h ../orbital/orbital.cpp ../lib/lib.h ../lib/lib.cpp  -llapack -lcblas && ./test.out  && gprof ./test.out gmon.out > testres.txt && less testres.txt 
