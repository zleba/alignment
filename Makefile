a.out: read.cpp
	g++ -g -O2  read.cpp  -larmadillo


armaInc=arma/install/include
armaLib=arma/install/lib64
mklLib=/opt/intel/2018/mkl/lib/intel64

check: check.cpp
#	g++ -std=c++11 -g -O3 -I$(armaInc)  -DARMA_DONT_USE_WRAPPER $^ pede/mpnum.o  pede/Dbandmatrix.o   -o $@ -Wl,-rpath,$(mklLib) -L$(mklLib) -lmkl_rt      -lgfortran -fopenmp
	g++ -std=c++11 -g -fopenmp  -O3 -I$(armaInc)  -DARMA_DONT_USE_WRAPPER $^ pede/mpnum.o  pede/Dbandmatrix.o   -o $@  -lopenblas  -llapack     -lgfortran 

#/home/radek/Downloads/millepede/target/pede.o  /home/radek/Downloads/millepede/target/mphistab.o

# /opt/intel/2019/mkl/lib/intel64/libmkl_rt.so  -L$(mklLib) -lmkl_rt 
