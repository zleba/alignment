
inc=inc
armaInc=arma/install/include
armaLib=arma/install/lib64
mklLib=/opt/intel/2018/mkl/lib/intel64
mklLib=/afs/desy.de/SL/7/x86_64/opt/intel/2020/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin

runFit: obj/read.o obj/fit.o
	g++ -std=c++11 -g -O3 -I$(armaInc) -I$(inc) -DARMA_DONT_USE_WRAPPER $^  -o $@ -Wl,-rpath,$(mklLib) -L$(mklLib) -lmkl_rt  # -lgfortran -fopenmp


obj/read.o: src/read.cc
	g++ -c  -std=c++11 -g -O3 -I$(armaInc) -I$(inc)  -DARMA_DONT_USE_WRAPPER $^  -o $@ 

obj/fit.o: src/fit.cc
	g++ -c  -std=c++11 -g -O3 -I$(armaInc) -I$(inc) -DARMA_DONT_USE_WRAPPER $^  -o $@ 





a.out: read.cpp
#	g++ -g -O2  read.cpp    -larmadillo
	g++ -std=c++11 -g -O3 -I$(armaInc)  -DARMA_DONT_USE_WRAPPER $^  -o $@ -Wl,-rpath,$(mklLib) -L$(mklLib) -lmkl_rt  # -lgfortran -fopenmp


check: check.cpp
#	g++ -std=c++11 -g -O3 -I$(armaInc)  -DARMA_DONT_USE_WRAPPER $^ pede/mpnum.o  pede/Dbandmatrix.o   -o $@ -Wl,-rpath,$(mklLib) -L$(mklLib) -lmkl_rt      -lgfortran -fopenmp
	g++ -std=c++11 -g -fopenmp  -O3 -I$(armaInc)  -DARMA_DONT_USE_WRAPPER $^ pede/mpnum.o  pede/Dbandmatrix.o   -o $@  -lopenblas  -llapack     -lgfortran 

#/home/radek/Downloads/millepede/target/pede.o  /home/radek/Downloads/millepede/target/mphistab.o

# /opt/intel/2019/mkl/lib/intel64/libmkl_rt.so  -L$(mklLib) -lmkl_rt 
