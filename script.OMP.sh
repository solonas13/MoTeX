#!/bin/bash
# use bash as shell
#$ -S /bin/bash
# join stout and stderr
#$ -j y
# change to working directory
#$ -cwd
# select parallel enviroment
#$ -pe mvapich 48
#$ -q maxmem.q,standard.q 
#$ -l h_rt=1:00:00
#$ -binding linear:48
# load needed modules
source /etc/profile.d/modules.sh
module load sge
#  For icc
 module load intel/compiler/64/12.1/2011_sp1.11.339 libraries libdatrie mpfr/3.1.1/icc/12.1
#  For gcc
# module load gcc libraries libdatrie
# launch job, hostfile and number of cores are set automagicaly
time ./motexOMP -a DNA -i ./dnc_subtilis_330-30.seq -o bg.motex -d 0 -q 20 -l 10 -e 2 -t 48
