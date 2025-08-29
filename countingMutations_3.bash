#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=4000              # Memory per cpu
#SBATCH -o cMut.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e cMut.%N.%j.err.log      # File to which STDERR will be written

##########################
##
## This script should be run as: 
## sbatch countingMutations_3.bash "twoLoci_rec0.5kb_cre1kb_s-0.0425_epist" 10000
##
##
## by Isabel Alves - Feb 2025
#########################

wrkDir="/shared/home/ialves/CRE_evolution/SLiM"

cd $wrkDir
#model, eg twoLoci_rec_cre1kb_s0.3_epist_13_g10000
m=$1
#generation
g=$2

echo "Collecting model: $m"
echo "Collecting generation: $g"

cd sims_Feb2025/$m/
#going over the replicates
for k in {1..100}; do
    
    grep m1 ${m}_${k}_g${g}.segMutations | wc -l >> ${m}_g${g}_m1_sMut.out
    grep m1 ${m}_${k}_g${g}.fixedMutations | wc -l >> ${m}_g${g}_m1_fMut.out
    
    grep m2 ${m}_${k}_g${g}.segMutations | wc -l >> ${m}_g${g}_m2_sMut.out
    grep m2 ${m}_${k}_g${g}.fixedMutations | wc -l >> ${m}_g${g}_m2_fMut.out

    grep m3 ${m}_${k}_g${g}.segMutations | wc -l >> ${m}_g${g}_m3_sMut.out
    grep m3 ${m}_${k}_g${g}.fixedMutations | wc -l >> ${m}_g${g}_m3_fMut.out

done
