#! /bin/sh

wget http://www-igm.univ-mlv.fr/~marsan/SMILEv1.47.tgz
wget http://www.lx.it.pt/~asmc/pub/software/RISO/alphabet
gunzip SMILEv1.47.tgz
tar -xvf SMILEv1.47.tar
rm SMILEv1.47.tar
cd SMILEv1.47/
make
mv ./SigStat/bin/e-smile_shuffling ../smile
cd ..
rm -r SMILEv1.47
mkdir SMILE
mv smile ./SMILE
mv alphabet ./SMILE
