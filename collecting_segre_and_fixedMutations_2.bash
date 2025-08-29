#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=4000              # Memory per cpu
#SBATCH -o slim.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e slim.%N.%j.err.log      # File to which STDERR will be written

##################
##
## This script should be run like:
## batch ollecting_segre_and_fixedMutations_2.bash "twoLoci_rec0.5kb_cre1kb_s-0.0425_epist" 10000
##
## by Isabel Alves - Feb 2025
##
##################



wrkDir="/shared/home/ialves/CRE_evolution/SLiM"

cd $wrkDir
#model
m=$1
#generation
g=$2

echo "Collecting model: $m"
echo "Collecting generation: $g"

cd sims_Feb2025/$m/
#going over the replicates
for k in {1..100}; do
    
    csplit -z ${m}_${k}.out /'#OUT: 10000 10000'/ '{*}'
    mv xx01 ${m}_${k}_g${g}.segMut
    mv xx02 ${m}_${k}_g${g}.fixedMutations
    rm -f xx*
 
    mkdir ${m}_${k}_g${g}_tmp
    csplit -z --prefix ${m}_${k}_g${g}_tmp/${m}_${k}_g${g}_ ${m}_${k}_g${g}.segMut /Genomes/ '{*}'

    if [ -f "${m}_${k}_g${g}_tmp/${m}_${k}_g${g}_00" ] && [ -f "${m}_${k}_g${g}_tmp/${m}_${k}_g${g}_01" ]; 
    then
        mv ${m}_${k}_g${g}_tmp/${m}_${k}_g${g}_00 ${m}_${k}_g$g.segMutations
        mv ${m}_${k}_g${g}_tmp/${m}_${k}_g${g}_01 ${m}_${k}_g$g.genomes

    else
        echo "ERROR: no file found." 1>&2
    fi
    rm -f ${m}_${k}_g${g}.segMut
done
rm -rf *_tmp/
cd 
