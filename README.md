# GEL GWAS

![](bin/covid_1_manhattan.png)

## Example usage
```bash
nextflow run main.nf \
  --plinkFile "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
  --phenoFile s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sample.phe \
  --phenoCol PHE \
  --vcfsList s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/vcfs.csv
```
