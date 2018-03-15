GCC  = g++ -O3
INCLUDE = -I/usr/include/eigen3/ -I./bemtool/

all: testMax clean

######################################################

test: test.o
	$(GCC) test.o -o test

test.o: test.cpp
	$(GCC) $(INCLUDE) -c test.cpp -o test.o

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

l_beti: l_beti.o
	$(GCC) l_beti.o -o l_beti

l_beti.o: l_beti.cpp
	$(GCC) $(INCLUDE) -c l_beti.cpp -o l_beti.o

######################################################

testMax: testMax.o
	$(GCC) testMax.o -o testMax

testMax.o: testMax.cpp
	$(GCC) $(INCLUDE) -c testMax.cpp -o testMax.o

######################################################



clean:
	rm *.o
