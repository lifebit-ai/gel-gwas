#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/gel-gwas
========================================================================================
 lifebit-ai/gel-gwas GWAS pipeline built for Genomics England using SAIGE
 #### Homepage / Documentation
 https://github.com/lifebit-ai/gel-gwas
----------------------------------------------------------------------------------------
*/

/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/
Channel
  .fromPath(params.phenoFile)
  .ifEmpty { exit 1, "Pheno file not found: ${params.phenoFile}" }
  .set { phenoCh }
Channel
  .fromFilePairs("${params.plinkFile}",size:3, flat : true)
  .ifEmpty { exit 1, "PLINK files not found: ${params.plinkFile}" }
  .set { plinkCh }
Channel
  .fromPath(params.vcfsList)
  .ifEmpty { exit 1, "Cannot find CSV VCFs file : ${params.vcfsList}" }
  .splitCsv(skip:1)
  .map { chr, vcf, index -> [file(vcf).simpleName, chr, file(vcf), file(index)] }
  .set { vcfsCh }

/*--------------------------------------------------
  Pre-GWAS filtering - download, filter and convert VCFs
---------------------------------------------------*/

process pre_gwas_filtering {
  tag "$name"
  publishDir "${params.outdir}/pre_gwas_filtering", mode: 'copy'

  input:
  set val(name), val(chr), file(vcf), file(index) from vcfsCh

  output:
  set val(name), val(chr), file("${name}_filtered.vcf.gz"), file("${name}_filtered.vcf.gz.csi") into filteredVcfsCh
  file("${name}_filtered.{bed,bim,fam}") into plinkTestCh

  script:
  // TODO: (High priority) Only extract needed individuals from VCF files with `bcftools -S samples.txt` - get from samples file?
  // TODO: (Not required) `bcftools -T sites_to_extract.txt`
  """
  # Download, filter and convert (bcf or vcf.gz) -> vcf.gz
  bcftools view -q ${params.qFilter} $vcf -Oz -o ${name}_filtered.vcf.gz
  bcftools index ${name}_filtered.vcf.gz
 
  # Create PLINK binary from vcf.gz
  plink2 \
    --make-bed \
    --set-missing-var-ids @:#,\\\$r,\\\$a \
    --vcf ${name}_filtered.vcf.gz \
    --out ${name}_filtered \
    --vcf-half-call m \
    --double-id \
    --set-hh-missing \
    --new-id-max-allele-len 60 missing
  """
 # # # # # # # # # # # #  
 #add arguments:
thres_m=$1
thres_HWE=$2

add argument (could be fixed) also for plate_keys/IDs for unrelated individuals:
/re_gecip/BRS/thanos/aggV2_GRM_data/aggV2_bedmerge_30KSNPs_labkeyV9_unrelatedKING.fam

for mock test simply exlude "--keep /re_gecip/BRS/thanos/aggV2_GRM_data/aggV2_bedmerge_30KSNPs_labkeyV9_unrelatedKING.fam" OR

replace with mock list of IDs (say 10 first).
FID IID
ID1 ID1
ID2 ID2


default for both = 1e-5
###################
#Commands here:
#Filter missingness
plink --bfile ${name}_filtered --pheno $phenofile --pheno-name $phenocol --allow-no-sex --test-missing midp --out $output_tag --1 --keep-allele-order
 
awk '$5 < $thres_m {print}' $output_tag.missing > $output_tag.missing_FAIL

#Filter HWE 
plink --bfile $plinkfile --pheno $phenofile --pheno-name $phenocol --allow-no-sex --hwe $thres_HWE midp --out $output_tag.misHWEfiltered --make-just-bim --exclude $output_tag.missing_FAIL --1 --keep-allele-order --keep /re_gecip/BRS/thanos/aggV2_GRM_data/aggV2_bedmerge_30KSNPs_labkeyV9_unrelatedKING.fam

bcftools view  $vcf -Oz -o ${name}_filtered.vcf.gz
#awk add multi-allelic filter

${name}_filtered2.vcf.gz
 
 # # # # # # # # # # # # #
}

/*--------------------------------------------------
  GWAS Analysis 1 with SAIGE - Fit the null mixed-model
---------------------------------------------------*/

process gwas_1_fit_null_glmm {
  tag "$plink_GRM_snps"
  publishDir "${params.outdir}/gwas_1_fit_null_glmm", mode: 'copy'

  input:
  set val(plink_GRM_snps), file(bed), file(bim), file(fam) from plinkCh
  each file(phenoFile) from phenoCh

  output:
  file "*" into fit_null_glmm_results
  file ("step1_${params.phenoCol}_out.rda") into rdaCh
  file ("step1_${params.phenoCol}.varianceRatio.txt") into varianceRatioCh

  script:
  """
  step1_fitNULLGLMM.R \
    --plinkFile=${plink_GRM_snps} \
    --phenoFile=${phenoFile} \
    --phenoCol=${params.phenoCol} \
    --sampleIDColinphenoFile=IID \
    --traitType=${params.traitType} \
    --outputPrefix=step1_${params.phenoCol}_out \
    --outputPrefix_varRatio=step1_${params.phenoCol} \
    --nThreads=${task.cpus} ${params.saigeStep1ExtraFlags}
  """
}

/*--------------------------------------------------
  GWAS Analysis 2 with SAIGE - Perform mixed-model association testing
---------------------------------------------------*/

process gwas_2_spa_tests {
  tag "$name"
  publishDir "${params.outdir}/gwas_2_spa_tests", mode: 'copy'

  input:
  set val(name), val(chr), file(vcf), file(index) from filteredVcfsCh
  each file(rda) from rdaCh
  each file(varianceRatio) from varianceRatioCh

  output:
  file "*" into results

  script:
  """
  step2_SPAtests.R \
    --vcfFile=${vcf} \
    --vcfFileIndex=${index} \
    --vcfField=GT \
    --chrom=${chr} \
    --minMAC=20 \
    --sampleFile=day0_covid.samples \
    --GMMATmodelFile=${rda} \
    --varianceRatioFile=${varianceRatio} \
    --SAIGEOutputFile=${params.phenoCol}.${name}.SAIGE.gwas.txt \
    --numLinesOutput=2 \
    --IsOutputAFinCaseCtrl=TRUE \
    --IsDropMissingDosages=FALSE \
    --IsOutputNinCaseCtrl=TRUE \
    --IsOutputHetHomCountsinCaseCtrl=TRUE
  """
}

###############
add scripts for plotting
manhattan plot + qqplot + CI

process plotting

DONT add an argument for user. File is ready:

${params.phenoCol}.${name}.SAIGE.gwas.txt
