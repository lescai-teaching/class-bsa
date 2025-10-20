
###########################################
## for the environment of variant calling #
###########################################


cd /workspaces/class-variantcalling

git clone https://github.com/lescai-teaching/class-bsa.git
git clone https://github.com/lescai-teaching/datasets-class-variantcalling.git
git clone https://github.com/lescai-teaching/reference_chr21.git
git clone https://github.com/lescai-teaching/reference_chr22.git

cd /workspaces/class-variantcalling/reference_chr22/
git lfs fetch --all
git lfs checkout

cd /workspaces/class-variantcalling

rm -rf */.git


###########################################
## for the environment of RNA sequencing  #
###########################################

cd /workspaces/class-rnaseq

git clone https://github.com/lescai-teaching/class-bsa.git
git clone https://github.com/lescai-teaching/datasets-class-rnaseq.git
git clone https://github.com/lescai-teaching/reference_chr21.git
git clone https://github.com/lescai-teaching/reference_chr22.git

cd /workspaces/class-rnaseq/reference_chr22/
git lfs fetch --all
git lfs checkout

cd /workspaces/class-rnaseq
rm -rf */.git