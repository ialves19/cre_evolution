# to check which runs did not work for each bootstrap 
# within BOOTSTRAP folder
b=98; 
for ite in `seq 1 100`; 
do 
    fileExists=`ls TwoDomest_EurOld_Taiw_stru_${b}_dSFS/TwoDomest_EurOld_Taiw_stru_${b}_${ite}/TwoDomest_EurOld_Taiw_stru_${b}_${ite}/TwoDomest_EurOld_Taiw_stru_${b}_${ite}.bestlhoods`; 
    if [ -z "$fileExists" ]; 
    then 
        echo "$ite not done"; 
    fi;
done

# cp est and tpl from the first bootstrap to all the other bootstraps 
# within BOOTSTRAP folder
for b in `seq 59 100`; 
    do 
        echo "Copying est and tpl file to directory TwoDomest_EurOld_Taiw_stru_${b}_dSFS"; 
        cp TwoDomest_EurOld_Taiw_stru_1_dSFS/TwoDomest_EurOld_Taiw_stru_1.est TwoDomest_EurOld_Taiw_stru_${b}_dSFS/TwoDomest_EurOld_Taiw_stru_${b}.est; 
        cp TwoDomest_EurOld_Taiw_stru_1_dSFS/TwoDomest_EurOld_Taiw_stru_1.tpl TwoDomest_EurOld_Taiw_stru_${b}_dSFS/TwoDomest_EurOld_Taiw_stru_${b}.tpl; 
done 

# check how many successful runs we have for a given bootstrap 
# used to fill in the table : Taiwanese_models_summary.xlsx
# within BOOTSTRAP folder
for boot in `seq 80 92`; 
do 
    echo "Analysing bootstrap #: ${boot}"; 
    find . -name  TwoDomest_EurOld_Taiw_stru_${boot}_*.bestlhoods | wc -l; 
done