
cd /workspaces/class-rnaseq
mkdir -p analysis
cd analysis

## instead of copying - we create symbolic links to the original reads
mkdir reads
cd reads
ln -s /workspaces/class-rnaseq/datasets-class-rnaseq/raw_data/* .

### check if you can execute salmon by typing "salmon" on the terminal
## sometimes it fails on RStudio terminal, on CodeSpaces while it works on GitPod
## export PATH=${PATH}:/usr/local/bin

## the index for the transcriptome is located in
## /workspaces/class-rnaseq/datasets_reference_only/trascriptome/chr21_transcripts_index

## now we can quantify all samples, by running a loop with salmon and the following


for sample in `ls *_1.fasta.gz`
do
index="/workspaces/class-rnaseq/datasets_reference_only/trascriptome/chr21_transcripts_index"
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
