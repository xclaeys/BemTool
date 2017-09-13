GCC  = g++ -O3 
INCLUDE = -I/usr/include/eigen3/ -I./bemtool/

all: test2D clean

######################################################

test2D: test2D.o
	$(GCC) test2D.o -o test2D

test2D.o: test2D.cpp
	$(GCC) $(INCLUDE) -c test2D.cpp -o test2D.o

######################################################

test3D: test3D.o
	$(GCC) test3D.o -o test3D

test3D.o: test3D.cpp
	$(GCC) $(INCLUDE) -c test3D.cpp -o test3D.o

######################################################

mtf: mtf.o
	$(GCC) mtf.o -o mtf

mtf.o: mtf.cpp
	$(GCC) $(INCLUDE) -c mtf.cpp -o mtf.o

######################################################



clean:
	rm *.o
