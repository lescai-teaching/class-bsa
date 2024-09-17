
### variant calling

cd /workspaces/class-variantcalling
mkdir -p analysis/variants
cd analysis/variants


## first single sample discovery

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /workspaces/class-variantcalling/analysis/alignment/normal_recal.bam \
   -O normal.g.vcf.gz \
   -ERC GVCF

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /workspaces/class-variantcalling/analysis/alignment/disease_recal.bam \
   -O disease.g.vcf.gz \
   -ERC GVCF

## then consolidate the 2 files

mkdir -p tmp

### on AMD64 this code ######
## combine the files into one
gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
      -V normal.g.vcf.gz \
      -V disease.g.vcf.gz \
      --genomicsdb-workspace-path compared_db \
      --tmp-dir /workspaces/class-variantcalling/analysis/variants/tmp \
      -L chr21

### on ARM64 (Mac M1 chip) this code
## combine the files into one
 gatk CombineGVCFs \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V normal.g.vcf.gz \
   -V disease.g.vcf.gz \
   -O cohort.g.vcf.gz

### on AMD64 this code ######
### finally we can call the genotypes jointly
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V gendb://compared_db \
   --dbsnp /workspaces/class-variantcalling/datasets_reference_only/gatkbundle/dbsnp_146.hg38_chr21.vcf.gz \
   -O results.vcf.gz



### on ARM64 (Mac M1 chip) this code
### finally we can call the genotypes jointly
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V cohort.g.vcf.gz \
   --dbsnp /workspaces/class-variantcalling/datasets_reference_only/gatkbundle/dbsnp_146.hg38_chr21.vcf.gz \
   -O results.vcf.gz