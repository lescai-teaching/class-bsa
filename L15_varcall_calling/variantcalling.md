# Resequencing Workflow - Calling the variants

The second step of our resequencing workflow deals with the identification of variant sites in the aligned data.

## Prepare your folders

As usual, we prepare our folders to keep the data in a tidy structure:

```{bash}
cd /config/workspace/variant_calling
mkdir -p variants
cd variants
```


## Identify variant sites

Like we have seen in class, the first step is identifying the variant sites and calling a first hypothesis of genotypes on a per-sample basis.

We therefore do this for each sample, first on the *normal* one

```{bash}
gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /config/workspace/variant_calling/alignment/normal_recal.bam \
   -O normal.g.vcf.gz \
   -ERC GVCF
```

Then on the sample we have indicated as *disease* case.


```{bash}
gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /config/workspace/variant_calling/alignment/disease_recal.bam \
   -O disease.g.vcf.gz \
   -ERC GVCF
```


## Combine data


Now the important part: we need to combine the previously generated datasets, in order to perform the joint-calling of the genotypes.

First we create a temporary folder, which is necessary for the computation

```{bash}
mkdir -p tmp
```
Then we put the data together, in a step that's necessary to create a local database of the samples we have called.


```{bash}
gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
    -V normal.g.vcf.gz \
    -V disease.g.vcf.gz \
    --genomicsdb-workspace-path compared_db \
    --tmp-dir /config/workspace/variant_calling/variants/tmp \
    -L chr21
```


## Calculate Genotypes


We then use the generated database in order to call the samples together:

```{bash}
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V gendb://compared_db \
   --dbsnp /config/workspace/datasets_class/reference/gatkbundle/dbsnp_146.hg38_chr21.vcf.gz \
   -O results.vcf.gz
```

The resulting VCF file will contain data on both individuals on separate columns.