mkdir -p markdown_exercise/{reads,variants}
cd markdown_exercise/reads
cp -v /home/rstudio/data/datasets_class/germline_calling/reads/* .

## to assess read quality we only get the first 5000 reads 
## otherwise it takes too long to read in python
bgzip -d -c normal_1.000+disease_0.000_1.fq.gz | head -n 20000 > normal_sample_1.fq
bgzip -d -c normal_1.000+disease_0.000_2.fq.gz | head -n 20000 > normal_sample_2.fq

cd ../variants
cp -v /home/rstudio/data/datasets_class/germline_calling/variants/*ann* .
cp -v /home/rstudio/data/variant_calling/variants/results_ann.vcf .

