### locate data you need:

## datasets are in
## {where-you-put-it}/datasets_class/rna_sequencing/raw_data

cd /home/rstudio/data
mkdir -p rnaseq_exercise
cd rnaseq_exercise

## instead of copying - we create symbolic links to the original reads
mkdir reads
cd reads
ln -s /home/rstudio/data/datasets_class/rna_sequencing/raw_data/* .

### check if you can execute salmon by typing "salmon" on the terminal
## in my case
export PATH=${PATH}:/opt/software/salmon-1.5.2/bin/

## the index for the transcriptome is located in
## /home/rstudio/data/datasets_class/reference/trascriptome/chr21_transcripts_index

## now we can quantify all samples, by running a loop with salmon and the following


for sample in `ls *_1.fasta.gz`
do
index="/home/rstudio/data/datasets_class/reference/trascriptome/chr21_transcripts_index"
name=${sample%_1.fasta.gz}
echo "quantifying $name"
salmon quant \
 -p 2 \
 -i $index \
 -l IU \
 -1 "${name}_1.fasta.gz" -2 "${name}_2.fasta.gz" \
 --validateMappings \
 -o "${name}.quant"
echo -e "$name done now\n"
done

### let's inspect a quantification file

cd sample_01.quant
head quant.sf

## more information on the format of the output
## https://salmon.readthedocs.io/en/latest/file_formats.html

