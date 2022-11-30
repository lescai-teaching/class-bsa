# Resequencing Workflow - Calling the variants

The second step of our resequencing workflow deals with the identification of variant sites in the aligned data.

## Prepare your folders

As usual, we prepare our folders to keep the data in a tidy structure:

```{bash}
cd /config/workspace/dati_vscode/variant_calling
mkdir -p variants
cd variants
```


## Identify variant sites

Like we have seen in class, the first step is identifying the variant sites and calling a first hypothesis of genotypes on a per-sample basis.

We therefore do this for each sample, first on the *normal* one

```{bash}
gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /config/workspace/dati_vscode/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /config/workspace/dati_vscode/variant_calling/alignment/normal_recal.bam \
   -O normal.g.vcf.gz \
   -ERC GVCF
```

Then on the sample we have indicated as *disease* case.


```{bash}
gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /config/workspace/dati_vscode/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /config/workspace/dati_vscode/variant_calling/alignment/disease_recal.bam \
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

**The following code only works on AMD64 computers - NOT on Apple Silicon**

```{bash}
gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
      -V normal.g.vcf.gz \
      -V disease.g.vcf.gz \
      --genomicsdb-workspace-path compared_db \
      --tmp-dir /config/workspace/dati_vscode/variant_calling/variants/tmp \
      -L chr21
```


## Calculate Genotypes


We then use the generated database in order to call the samples together:

```{bash}
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /config/workspace/dati_vscode/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V gendb://compared_db \
   --dbsnp /config/workspace/dati_vscode/datasets_reference_only/gatkbundle/dbsnp_146.hg38_chr21.vcf.gz \
   -O results.vcf.gz
```

The resulting VCF file will contain data on both individuals on separate columns.


**The following code is specific for Apple Silicon**

First 

```{bash}
 gatk CombineGVCFs \
   -R /config/workspace/dati_vscode/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V normal.g.vcf.gz \
   -V disease.g.vcf.gz \
   -O cohort.g.vcf.gz
```


Then:

```{bash}
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /config/workspace/dati_vscode/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V cohort.g.vcf.gz \
   --dbsnp /config/workspace/dati_vscode/datasets_reference_only/gatkbundle/dbsnp_146.hg38_chr21.vcf.gz \
   -O results.vcf.gz
```