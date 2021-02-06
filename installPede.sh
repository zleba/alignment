#svn checkout http://svnsrv.desy.de/public/MillepedeII/tags/V04-06-00 pede
#cd pede
#make pede

rm mp2-04-08-03.tgz
wget https://www.desy.de/~kleinwrt/MP2/tar/mp2-04-08-03.tgz
tar -xzvf mp2-04-08-03.tgz
mv mp2-04-08-03 pede
cd pede
make 
