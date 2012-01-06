vmc: lib/lib.cpp lib/lib.h vmc.cpp
	c++ -o vmc.out  lib/lib.cpp vmc.cpp

run: vmc
	make vmc
	/.vmc.out
