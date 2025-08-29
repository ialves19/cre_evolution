
## in the cluster
mkdir toCopy

find . -name '*_m2_fMut.out' -exec cp {} toCopy/  \;
find . -name '*_m3_fMut.out' -exec cp {} toCopy/  \;
find . -name 'sumStats_*' -exec cp {} toCopy/ \;
find . -name 'pi_*' -exec cp {} toCopy/ \;

## copiez le dossier toCopy dans votre machine ("en local")

## locally
# dans le dossier ou vous gardez les résultats
model=("twoLoci_rec500kb_cre1kb_s-0.001" "twoLoci_rec500kb_cre1kb_s-0.001_epist" "twoLoci_rec500kb_cre1kb_s-0.01"
 "twoLoci_rec500kb_cre1kb_s-0.01_epist" "twoLoci_rec500kb_cre1kb_s-0.0425" "twoLoci_rec500kb_cre1kb_s-0.0425_epist"
 "twoLoci_rec500kb_cre1kb_s-0.1" "twoLoci_rec500kb_cre1kb_s-0.1_epist" "twoLoci_rec500kb_cre1kb_s-0.3"
 "twoLoci_rec500kb_cre1kb_s-0.3_epist" "twoLoci_rec500kb_cre5kb_s-0.001" "twoLoci_rec500kb_cre5kb_s-0.001_epist"
 "twoLoci_rec500kb_cre5kb_s-0.01" "twoLoci_rec500kb_cre5kb_s-0.01_epist" "twoLoci_rec500kb_cre5kb_s-0.0425"
 "twoLoci_rec500kb_cre5kb_s-0.0425_epist" "twoLoci_rec500kb_cre5kb_s-0.1" "twoLoci_rec500kb_cre5kb_s-0.1_epist"
 "twoLoci_rec500kb_cre5kb_s-0.3" "twoLoci_rec500kb_cre5kb_s-0.3_epist" "twoLoci_rec500kb_cre50kb_s-0.001"
 "twoLoci_rec500kb_cre50kb_s-0.001_epist" "twoLoci_rec500kb_cre50kb_s-0.01" "twoLoci_rec500kb_cre50kb_s-0.01_epist"
 "twoLoci_rec500kb_cre50kb_s-0.0425" "twoLoci_rec500kb_cre50kb_s-0.0425_epist" "twoLoci_rec500kb_cre50kb_s-0.1"
 "twoLoci_rec500kb_cre50kb_s-0.1_epist" "twoLoci_rec500kb_cre50kb_s-0.3" "twoLoci_rec500kb_cre50kb_s-0.3_epist")
for m in ${model[@]}; 
do 
    mkdir $m
    cp toCopy/sumStats_${m}_haps.txt $m/
    cp toCopy/sumStats_${m}_perSite.txt $m/
    cp toCopy/${m}_g10000_m2_fMut.out $m/
    cp toCopy/${m}_g10000_m3_fMut.out $m/
done

