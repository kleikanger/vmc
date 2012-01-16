vmc: lib/lib.cpp lib/lib.h vmc_be.cpp slaterMatrix/slaterMatrix.h slaterMatrix/slaterMatrix.cpp orbital/orbital.h orbital/orbital.cpp 
	c++ -o vmc_be.out -Wall vmc_be.cpp slaterMatrix/slaterMatrix.cpp slaterMatrix/slaterMatrix.h orbital/orbital.h orbital/orbital.cpp lib/lib.cpp 

clang: vmc 
	clang -o vmc_be slaterMatrix/slaterMatrix.cpp slaterMatrix/slaterMatrix.h orbital/orbital.h orbital/orbital.cpp lib/lib.cpp 
	
run: vmc
	./vmc_be.out

