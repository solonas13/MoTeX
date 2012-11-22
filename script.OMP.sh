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
#$ -l h_rt=2:00:00
#$ -binding linear:48
# load needed modules
source /etc/profile.d/modules.sh
module load sge
#  For icc
module load intel/compiler/64/12.1/2011_sp1.11.339 libraries libdatrie
#  For gcc
#  module load gcc libraries libdatrie
# launch job, hostfile and number of cores are set automagicaly
./motexOMP -a DNA -i ../data/dnc_subtilis_330-30.seq -o output.motex -d 1 -q 10 -l 13 -e 2 -n 100 -t 48
