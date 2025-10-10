
#### ANNOTATE THE SAMPLE

mkdir -p /workspaces/class-variantcalling/analysis/variants/cache
cd /workspaces/class-variantcalling/analysis/variants

### to execute snpeff we need to contain the memory
snpEff -Xmx4g ann -v hg38 results.vcf.gz >results_ann.vcf


### filter variants

cat results_ann.vcf | grep "#CHROM" | cut -f 10-

## use SNPsift to filter those with impact HIGH where the control is homozygous for reference and case is either 
## heterozygous or homozygous for the alternative allele

SnpSift filter "ANN[0].IMPACT = 'HIGH' && isRef(GEN[1]) && (isVariant(GEN[0]) | isHom(GEN[0]))" results_ann.vcf >filtered_variants.vcf

SnpSift extractFields \
-s "," -e "." \
filtered_variants.vcf \
"CHROM" "POS" "ID" "GEN[disease].GT" "GEN[normal].GT" ANN[*].GENE ANN[*].EFFECT \
| column -t

SnpSift extractFields \
-s "," -e "." \
filtered_variants.vcf \
"CHROM" "POS" "ID" "REF" "ALT" "GEN[*].GT" ANN[0].GENE ANN[0].EFFECT \
| column -t