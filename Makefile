a.out: read.cpp
	g++ -g -O2  read.cpp  -larmadillo

check: check.cpp
	g++ -g -O3  -march=native $^ /home/radek/Downloads/millepede/target/mpnum.o  /home/radek/Downloads/millepede/target/Dbandmatrix.o   -o $@  -larmadillo -lgfortran -fopenmp

#/home/radek/Downloads/millepede/target/pede.o  /home/radek/Downloads/millepede/target/mphistab.o
