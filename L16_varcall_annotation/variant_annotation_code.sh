
#### ANNOTATE THE SAMPLE

mkdir -p /workspaces/class-variantcalling/analysis/variants/cache
cd /workspaces/class-variantcalling/analysis/variants

## download hg38 (UCSC) version of database
snpEff download -v hg38 -dataDir /workspaces/class-variantcalling/analysis/variants/cache

### to execute snpeff we need to contain the memory
snpEff -Xmx4g ann -dataDir /workspaces/class-variantcalling/analysis/variants/cache -v hg38 results.vcf.gz >results_ann.vcf


### filter variants

cat results_ann.vcf | grep "#CHROM" | cut -f 10-

grep "#" results_ann.vcf >filtered_variants.vcf
cat results_ann.vcf | grep HIGH | perl -nae 'if($F[10]=~/0\/0/ && $F[9]=~/1\/1/){print $_;}' >>filtered_variants.vcf
cat results_ann.vcf | grep HIGH | perl -nae 'if($F[10]=~/0\/0/ && $F[9]=~/0\/1/){print $_;}' >>filtered_variants.vcf

sudo conda install bioconda::snpsift

SnpSift extractFields \
-s "," -e "." \
filtered_variants.vcf \
"CHROM" "POS" "ID" "GEN[disease].GT" "GEN[normal].GT" ANN[*].GENE ANN[*].EFFECT

SnpSift extractFields \
-s "," -e "." \
filtered_variants.vcf \
"CHROM" "POS" "ID" "GEN[*].GT" ANN[0].GENE ANN[0].EFFECT