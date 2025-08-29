#!/bin/bash

#SBATCH -p fast                 # Partition to submit to (fast / slow)
#SBATCH -n 1                   # Number of task to start in parallel
#SBATCH -c 8                    # Number of cpu per task (per default 1 cpu with 2Go Ram per task)
#SBATCH --mem-per-cpu=2000              # Memory per cpu
#SBATCH -o pixy.%N.%j.out.log      # File to which STDOUT will be written
#SBATCH -e pixy.%N.%j.err.log      # File to which STDERR will be written

start=`date +%s`

conda activate pixy
module load bcftools

inputDir="/shared/home/ialves/demoHist_yeast3039/03-data"
vcfMatrix="full3039Matrix.AllPositions.DP10.GQ20.Mind20.99pNonMiss.ExcHet99.SNPs.vcf.gz"
outDir="/shared/home/ialves/demoHist_yeast3039/04-analysis/sumStats_pixy"
sampleFile="/shared/home/ialves/demoHist_yeast3039/04-analysis/ABBA-BABA/1279strainClades.txt"

cat $sampleFile | awk 'BEGIN{FS="\t"}{if ($2 != "NA") print $1}' | head -n -1 > $outDir/1107samplesABBA-BABA.txt
cat $sampleFile | awk 'BEGIN{FS="\t"}{if ($2 != "NA") print $0}' | head -n -1 > $outDir/1107samplesABBA-BABA_POPS.txt

# count initial nb of variable sites 
bcftools +counts $inputDir/$vcfMatrix

# remove everything that is not biallelic snps 
bcftools view -M2 -Ou $inputDir/$vcfMatrix | bcftools view -S $outDir/1107samplesABBA-BABA.txt -Oz -o $outDir/1107samples.AllPositionsMax2alleles.vcf.gz --threads 8

# count initial nb of variable sites 
bcftools +counts $outDir/1107samples.AllPositionsMax2alleles.vcf.gz

# indexing vcf 
bcftools index -t $outDir/1107samples.AllPositionsMax2alleles.vcf.gz --threads 8

# compute Pi xy and Dxy
pixy --stats pi fst dxy \
--vcf $outDir/1107samples.AllPositionsMax2alleles.vcf.gz \
--populations $outDir/1107samplesABBA-BABA_POPS.txt \
--output_folder $outDir \ 
--window_size 10000 \
--n_cores 8


end=`date +%s`
runtime=$((end-start))
#days
D=$((runtime / 60 / 60 / 24))
# all hours
H=$((runtime / 60 / 60))
H=$((H % 24))
# all minutes
M=$((runtime / 60))
M=$((M % 60))
# all seconds
S=$((runtime % 60))
echo "This task took: $D days, $H hours, $M minutes and $S seconds."