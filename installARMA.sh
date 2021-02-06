
version=10.2.1

pwd=$PWD
mkdir -p $pwd/arma/install
cd $pwd/arma
#wget  https://sourceforge.net/projects/arma/files/armadillo-8.500.1.tar.xz/download  && tar xf download
wget  http://sourceforge.net/projects/arma/files/armadillo-${version}.tar.xz  && tar xf *.tar.xz
cd $pwd/arma/armadillo-${version}

cmake . -DCMAKE_INSTALL_PREFIX:PATH=$pwd/arma/install
make
make install
