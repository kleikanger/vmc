all: *.cpp lib/lib.cpp
	c++ -o main.exe gaussianQuadrature.cpp eigenPack.cpp main.cpp lib/lib.cpp

