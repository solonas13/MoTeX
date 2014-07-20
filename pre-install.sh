#! /bin/sh

tar -xvf SMILEv1.47.tgz
xz -d libdatrie_0.2.8.orig.tar.xz
tar -xvf libdatrie_0.2.8.orig.tar
tar -xvf mpfr-3.1.2.tar.gz

cd SMILEv1.47/
make
mv ./SigStat/bin/e-smile_shuffling ../smile
cd ..
mkdir SMILE
mv smile ./SMILE
mv alphabet ./SMILE

cd libdatrie-0.2.8
./configure --prefix="$(pwd)"/libdatrie
make
make install
mv libdatrie ..
cd ..

cd mpfr-3.1.2
./configure --prefix="$(pwd)"/libmpfr
make
make install
mv libmpfr ..
cd ..
