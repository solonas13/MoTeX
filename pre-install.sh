#! /bin/sh

cd SMILEv1.47/
make
mv ./SigStat/bin/e-smile_shuffling ../smile
cd ..
mkdir SMILE
mv smile ./SMILE
mv alphabet ./SMILE

cd libdatrie-0.2.8
./configure --prefix="$(pwd)"/libdatrie
aclocal
automake
make
make install
mv libdatrie ..
cd ..

cd mpfr-3.1.2
./configure --prefix="$(pwd)"/libmpfr
aclocal
automake
make
make install
mv libmpfr ..
cd ..
