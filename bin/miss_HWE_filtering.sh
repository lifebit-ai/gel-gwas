#Missingness & HWE filters
#module load plink/1.9

thres_m=$1
thres_HWE=$2
plinkfile=$3
phenofile=$4
phenocol=$5
output_tag=$6

#Filter missingness
plink --bfile $plinkfile --pheno $phenofile --pheno-name $phenocol --allow-no-sex --test-missing midp --out $output_tag --1 --keep-allele-order
 
awk '$5 < $thres_m {print}' $output_tag.missing > $output_tag.missing_FAIL

#Filter HWE 
plink --bfile $plinkfile --pheno $phenofile --pheno-name $phenocol --allow-no-sex --hwe $thres_HWE midp --out $output_tag.misHWEfiltered --make-just-bim --exclude $output_tag.missing_FAIL --1 --keep-allele-order --keep /re_gecip/BRS/thanos/aggV2_GRM_data/aggV2_bedmerge_30KSNPs_labkeyV9_unrelatedKING.fam
