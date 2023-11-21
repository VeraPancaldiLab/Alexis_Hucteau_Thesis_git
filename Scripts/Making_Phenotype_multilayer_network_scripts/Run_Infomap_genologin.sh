#!/bin/bash
#SBATCH -p workq
#SBATCH -t 01-00:00:00 #Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds".
#SBATCH --cpus-per-task 16
#SBATCH --mem=100G

#Load modules
#Need gcc-11.3.0
module load compiler/gcc-11.3.0

module load bioinfo/infomap-v2.7.1

mkdir /home/ahucteau/work/Multilayer/R
infomap /home/ahucteau/work/Multilayer/R_multilayer_infomap.edges \
--seed 12345 -N 5 -M 16 -L 16 --inner-parallelization \
-d -v -f directed /home/ahucteau/work/Multilayer/R
