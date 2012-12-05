#!/bin/bash
# use bash as shell
#$ -S /bin/bash
# join stout and stderr
#$ -j y 
#$ -R y
# change to working directory
#$ -cwd
# select parallel enviroment
#$ -pe mvapich 1056
#$ -q maxmem.q,standard.q 
#$ -l h_rt=24:00:00
#$ -binding linear:48
# load needed modules
source /etc/profile.d/modules.sh
module load sge

#  For icc
module load hydra intel/compiler/64/12.1/2011_sp1.11.339 mvapich2/intel/64/1.6-qlc libraries libdatrie mpfr/3.1.1/icc/12.1
#  For gcc
#  module load hydra gcc mpich2/ge/gcc/64/1.4.1p1 libraries libdatrie

# launch job, hostfile and number of cores are set automagicaly
mpistart ./motexMPI -a DNA -i ../data/hsapiens.full_gene_set.clean -o hsapiens.full_gene_set.clean.gcc.motex -d 0 -q 10 -l 12 -e 2
mpistart ./motexMPI -a DNA -i ../data/hsapiens.full_gene_set.clean -o hsapiens.full_gene_set.clean.gcc.motex -d 0 -q 10 -l 15 -e 3
mpistart ./motexMPI -a DNA -i ../data/hsapiens.full_gene_set.clean -o hsapiens.full_gene_set.clean.gcc.motex -d 0 -q 10 -l 20 -e 4
