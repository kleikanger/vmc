test: slaterMatrix.cpp slaterMatrix.h ../orbital/orbital.cpp ../orbital/orbital.h ../lib/lib.cpp test.cpp ../lib/lib.h 
	c++ -o test.out -Wall  test.cpp slaterMatrix.cpp slaterMatrix.h ../orbital/orbital.h ../orbital/orbital.cpp ../lib/lib.h ../lib/lib.cpp ../ipdist/ipdist.h ../ipdist/ipdist.cpp  -llapack -lcblas

run: test
	make test && ./test.out

profile: test
	g++ -pg -g2 -o test.out test.cpp slaterMatrix.cpp slaterMatrix.h ../orbital/orbital.h ../orbital/orbital.cpp ../lib/lib.h ../lib/lib.cpp ../ipdist/ipdist.h ../ipdist/ipdist.cpp  -llapack -lcblas && ./test.out  && gprof ./test.out gmon.out > testres.txt && less testres.txt 


test2: slaterMatrix.cpp slaterMatrix.h ../orbital/orbital.cpp ../orbital/orbital.h ../lib/lib.cpp test2.cpp ../lib/lib.h 
	c++ -O3 -o test2.out -Wall  test2.cpp slaterMatrix.cpp slaterMatrix.h ../orbital/orbital.h ../orbital/orbital.cpp ../lib/lib.h ../lib/lib.cpp ../ipdist/ipdist.h ../ipdist/ipdist.cpp  -llapack -lcblas && ./test2.out
