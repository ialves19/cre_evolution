#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=4000              # Memory per cpu
#SBATCH -o slim.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e slim.%N.%j.err.log      # File to which STDERR will be written

wrkDir="/shared/home/ialves/CRE_evolution/SLiM"
outputDir="/shared/home/ialves/CRE_evolution/SLiM/sims_Feb2025"

cd $wrkDir
echo "$1 replicate nb: $2"

if [ ! -d $outputDir/$1 ]; 
then
mkdir $outputDir/$1
#mv $outputDir/${1}_${2}.out $outputDir/$1
#else 
#mv $outputDir/${1}_${2}.out $outputDir/$1
fi 

##provide the name without ext
./slim $outputDir/Models/$1.txt > $outputDir/$1/${1}_${2}.out 


