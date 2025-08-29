#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=4000              # Memory per cpu
#SBATCH -o slim.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e slim.%N.%j.err.log      # File to which STDERR will be written

wrkDir="/home/ialves/SLiM"

cd $wrkDir
echo "$1"

cd sims/$1/
for k in {1..25}; do
    csplit -z ${1}_${k}.out /OUT/ '{*}'
    for i in {1..5000}; do 
        if [ $i -lt 10 ]; then
            grep m1 xx0$i | wc -l >> m1_mut.out
            grep m2 xx0$i | wc -l >> m2_mut.out
            grep m3 xx0$i | wc -l >> m3_mut.out
        else 
            grep m1 xx$i | wc -l >> m1_mut.out
            grep m2 xx$i | wc -l >> m2_mut.out
            grep m3 xx$i | wc -l >> m3_mut.out
        fi
    done
    rm -f xx*
    mv m1_mut.out ${1}_${k}_m1_mut.out
    mv m2_mut.out ${1}_${k}_m2_mut.out
    mv m3_mut.out ${1}_${k}_m3_mut.out
done