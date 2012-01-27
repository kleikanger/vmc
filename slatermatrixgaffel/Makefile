test: slaterMatrix.cpp slaterMatrix.h ../orbital/orbital.cpp ../orbital/orbital.h ../lib/lib.cpp test.cpp ../lib/lib.h 
	c++ -o test.out -Wall test.cpp slaterMatrix.cpp slaterMatrix.h ../orbital/orbital.h ../orbital/orbital.cpp ../lib/lib.h ../lib/lib.cpp  -llapack -lcblas

run: test
	make test && ./test.out
