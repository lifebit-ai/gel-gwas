#Autosomes GWAS filtering
for anc in ALL EUR SAS AFR; do
    for ANA in phen_ANA_C1_V2 phen_ANA_C2_V2;
    do
        >/re_gecip/BRS/gwas_workflow/GWASQC_filter_lists/aggV2_q0_005_wholecohortR9_${ANA}_${anc}.misHWEfiltered
        
        echo -e '#!/bin/bash' >/re_gecip/BRS/gwas_workflow/scripts/GWASQC_filter_lists_submission/filter_HWEmiss_pegasus_submission_${ANA}_${anc}.sh
        echo "#BSUB -q short
#BSUB -P bio
#BSUB -o /re_gecip/BRS/gwas_workflow/logfiles/GWASQC_filter_list_submission/filter_HWEmiss_pegasus_submission_${ANA}_${anc}.out
#BSUB -e /re_gecip/BRS/gwas_workflow/logfiles/GWASQC_filter_list_submission/filter_HWEmiss_pegasus_submission_${ANA}_${anc}.err
module load plink/1.9
chunklist=\$(cat /re_gecip/BRS/thanos/covid/chunklist)

for chunk in \$chunklist
do
sh /re_gecip/BRS/gwas_workflow/scripts/GWAS_QC_filtering/miss_HWE_filtering.sh 0.00001 0.00001 /re_gecip/BRS/gwas_workflow/plink_chunk_files/aggV2_q0_005_\${chunk}_wholecohortR9.filteredPASS /re_gecip/BRS/gwas_workflow/phenofiles/covid_2020_06_26_2020-06-30_14-58-48_aggV2_cohortR9_${anc}.phen ${ANA} /re_gecip/BRS/gwas_workflow/plink_chunk_files/aggV2_q0_005_\${chunk}_wholecohortR9_${ANA}_${anc}
cat /re_gecip/BRS/gwas_workflow/plink_chunk_files/aggV2_q0_005_\${chunk}_wholecohortR9_${ANA}_${anc}.misHWEfiltered.bim >> /re_gecip/BRS/gwas_workflow/GWASQC_filter_lists/aggV2_q0_005_wholecohortR9_${ANA}_${anc}.misHWEfiltered  
done
" >>/re_gecip/BRS/gwas_workflow/scripts/GWASQC_filter_lists_submission/filter_HWEmiss_pegasus_submission_${ANA}_${anc}.sh
        bsub </re_gecip/BRS/gwas_workflow/scripts/GWASQC_filter_lists_submission/filter_HWEmiss_pegasus_submission_${ANA}_${anc}.sh

    done
done
