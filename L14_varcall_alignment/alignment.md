# Resequencing workflow - Alignments

This is the first of 3 parts dedicated to resequencing, i.e. determining the sequence variation of individuals on a reference-based workflow.

## Before you start

If you haven't done that already, you should first download the dataset required by cloning the appropriate git repository, as below:

```{bash}
cd /config/workspace/dati_vscode
git clone https://github.com/lescai-teaching/datasets_bsa-2022.git
git clone https://github.com/lescai-teaching/datasets_reference_only.git
```

Should this encounter problems (the dataset is minimal but several files need to be downloaded), you can download the packaged repository as follows:

```{bash}
cd /config/workspace/dati_vscode
wget https://github.com/lescai-teaching/datasets_bsa-2022/archive/refs/tags/1.0.0.tar.gz
tar -xvzf 1.0.0.tar.gz
```


## Prepare your folders

We are going to create a structured folder where to perform the variant calling:

```{bash}
mkdir -p variant_calling
cd variant_calling
mkdir -p raw_data
cd raw_data
```

Then, we are *not* copying the data folder: instead, we are creating symbolic links:


```{bash}
ln -s /config/workspace/dati_vscode/datasets_bsa-2022/germline_calling/reads/*.gz .
```

We can then prepare another subfolder where to save the results of this exercise:

```{bash}
cd ..
mkdir -p alignment
cd alignment
```


## Align the reads to the reference

We use **bwa** in order to align the reads to a subset of the Human genome, limited to chromosome 21.

We do that separately for each of the samples: first the normal sample

```{bash}
bwa mem \
-t 2 \
-R "@RG\tID:sim\tSM:normal\tPL:illumina\tLB:sim" \
/config/workspace/dati_vscode/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
/config/workspace/dati_vscode/variant_calling/raw_data/normal_1.000+disease_0.000_1.fq.gz \
/config/workspace/dati_vscode/variant_calling/raw_data/normal_1.000+disease_0.000_2.fq.gz \
| samtools view -@ 8 -bhS -o normal.bam -
```

Then our simulated disease case sample:


```{bash}
bwa mem \
-t 2 \
-R "@RG\tID:sim\tSM:disease\tPL:illumina\tLB:sim" \
/config/workspace/dati_vscode/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
/config/workspace/dati_vscode/variant_calling/raw_data/normal_0.000+disease_1.000_1.fq.gz \
/config/workspace/dati_vscode/variant_calling/raw_data/normal_0.000+disease_1.000_2.fq.gz \
| samtools view -@ 8 -bhS -o disease.bam -
```

At the end of both computation, we will have a **bam** file.

## Sort and index BAM files

The alignments need first to be sorted (this also saves space on disk)

```{bash}
samtools sort -o normal_sorted.bam normal.bam
samtools sort -o disease_sorted.bam disease.bam
```

and *indexed*: the index is very important, as it allows accessing to slices of the data by genomics coordinates.

```{bash}
samtools index normal_sorted.bam
samtools index disease_sorted.bam
```


## Mark duplicates

You might refer to the class slides, to better understand what does *marking duplicates* mean: we perform this with GATK separately on each sample. First the one we called *normal*, then the one we called *disease*:

```{bash}
gatk MarkDuplicates \
-I normal_sorted.bam \
-M normal_metrics.txt \
-O normal_md.bam

gatk MarkDuplicates \
-I disease_sorted.bam \
-M disease_metrics.txt \
-O disease_md.bam
```


## Base Quality Score Recalibration

Refer to the class slides for the meaning of *base recalibration*. Here we perform it in 2 steps.

### Calculate recalibration

First we calculate the shift in quality scores, and create a recalibration table:

```{bash}
gatk BaseRecalibrator \
   -I normal_md.bam \
   -R /config/workspace/dati_vscode/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   --known-sites /config/workspace/dati_vscode/datasets_reference_only/gatkbundle/dbsnp_144.hg38_chr21.vcf.gz \
   --known-sites /config/workspace/dati_vscode/datasets_reference_only/gatkbundle/Mills_and_1000G_gold_standard.indels.hg38_chr21.vcf.gz \
   -O normal_recal_data.table

gatk BaseRecalibrator \
   -I disease_md.bam \
   -R /config/workspace/dati_vscode/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   --known-sites /config/workspace/dati_vscode/datasets_reference_only/gatkbundle/dbsnp_144.hg38_chr21.vcf.gz \
   --known-sites /config/workspace/dati_vscode/datasets_reference_only/gatkbundle/Mills_and_1000G_gold_standard.indels.hg38_chr21.vcf.gz \
   -O disease_recal_data.table
```

### Apply recalibration to alignments

Then, we use the recalibration table to modify the quality scores in the BAM files, and write new recalibrated files:


```{bash}
gatk ApplyBQSR \
   -R /config/workspace/dati_vscode/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I normal_md.bam \
   --bqsr-recal-file normal_recal_data.table \
   -O normal_recal.bam

gatk ApplyBQSR \
   -R /config/workspace/dati_vscode/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I disease_md.bam \
   --bqsr-recal-file disease_recal_data.table \
   -O disease_recal.bam
```