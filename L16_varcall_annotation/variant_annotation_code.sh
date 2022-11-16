
#### ANNOTATE THE SAMPLE

## download hg38 (UCSC) version of database
snpEff download -v hg38

### to execute snpeff we need to contain the memory
snpEff -Xmx4g ann -v hg38 results.vcf.gz >results_ann.vcf
