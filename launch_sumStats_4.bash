#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 1                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=4000              # Memory per cpu
#SBATCH -o sumStats.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e sumStats.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

#################
##
## This script should be run as:
## sbatch launch_sumStats_4.bash "/shared/home/ialves/CRE_evolution/SLiM/sims_Feb2025" "twoLoci_rec0.5kb_cre1kb_s-0.3"
##
##
## by Isabel Alves - Feb2025
##
################

wrkDir=$1
model=$2

# Before running this script please activate the conda environment where you have tidyverse and ggplot2
# conda env list
# conda activate XXXXX
#
#for the computingSumStats_SLIM.R you need to provide the folder where the model-specific
# folders are ex: /home/ialves/SLiM/sims 
# and the model name 

Rscript --vanilla computingSumStats_SLIM.R $wrkDir $model


end=`date +%s`
runtime=$((end-start))
echo $runtime
