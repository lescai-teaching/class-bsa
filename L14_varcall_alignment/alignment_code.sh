## clone repository
cd dati_vscode
git clone https://github.com/lescai-teaching/datasets_bsa-2022.git
git clone https://github.com/lescai-teaching/datasets_reference_only.git

## sym link so we do not change the repository itself

mkdir -p variant_calling
cd variant_calling
mkdir -p raw_data
cd raw_data
ln -s /config/workspace/dati_vscode/datasets_bsa-2022/germline_calling/reads/*.gz .
cd ..
mkdir -p alignment
cd alignment

## now we can perform the alignment with BWA

bwa mem \
-t 2 \
-R "@RG\tID:sim\tSM:normal\tPL:illumina\tLB:sim" \
/config/workspace/dati_vscode/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
/config/workspace/dati_vscode/variant_calling/raw_data/normal_1.000+disease_0.000_1.fq.gz \
/config/workspace/dati_vscode/variant_calling/raw_data/normal_1.000+disease_0.000_2.fq.gz \
| samtools view -@ 8 -bhS -o normal.bam -

## Real time: 152.349 sec; CPU: 257.270 sec

bwa mem \
-t 2 \
-R "@RG\tID:sim\tSM:disease\tPL:illumina\tLB:sim" \
/config/workspace/dati_vscode/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
/config/workspace/dati_vscode/variant_calling/raw_data/normal_0.000+disease_1.000_1.fq.gz \
/config/workspace/dati_vscode/variant_calling/raw_data/normal_0.000+disease_1.000_2.fq.gz \
| samtools view -@ 8 -bhS -o disease.bam -

## Real time: 157.120 sec; CPU: 273.662 sec


# sort the bam file
samtools sort -o normal_sorted.bam normal.bam
samtools sort -o disease_sorted.bam disease.bam

# index the bam file
samtools index normal_sorted.bam
samtools index disease_sorted.bam


# Marking duplicates

gatk MarkDuplicates \
-I normal_sorted.bam \
-M normal_metrics.txt \
-O normal_md.bam

gatk MarkDuplicates \
-I disease_sorted.bam \
-M disease_metrics.txt \
-O disease_md.bam


### recalibrating

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


#### Apply recalibration

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
