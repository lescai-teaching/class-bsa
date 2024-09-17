library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

setwd("/workspaces/class-rnaseq/analysis")
load("deseq2_analysis.RData")


###################################
## EXTRACT SIGNIFICANT GENES #####
###################################

universe <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = keys(org.Hs.eg.db),
                                  columns = c('ENTREZID','SYMBOL','ENSEMBL','ENSEMBLTRANS'),
                                  keytype = 'ENTREZID')

sig_genes <- resdata$gene[which(resdata$padj<0.05)]
entrez_genes_sig <- unique(universe[which(universe$ENSEMBL %in% sig_genes),]$ENTREZID)

pvalue_ens_genes <- resdata$padj[which(resdata$padj<0.05)]
names(pvalue_ens_genes)<-sig_genes

pvalue_entrez_genes <- resdata$padj[which(resdata$padj<0.05)]
names(pvalue_entrez_genes) <- entrez_genes_sig


###################################
## ENRICH GO ANALYSIS #####
###################################

ego <- enrichGO( gene = sig_genes,
                 universe = unique(tx2gene$GENEID),
                 OrgDb = org.Hs.eg.db,
                 keyType = 'ENSEMBL',
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05)

## writeLines(summary(ego), "enrich_GO_results.txt")

pdf("plots_ego.pdf")
dotplot(ego, showCategory=30)
dev.off()

pdf("plots_network-ego.pdf")
cnetplot(ego, foldChange=resdata$log2FoldChange[which(resdata$padj<0.5)])
dev.off()

###################################
## DISGNET ANALYSIS #####
###################################

## we need to unpack the file we have to read first
## use the terminal
## cd /home/rstudio/data/datasets_class/reference/trascriptome/
## gunzip all_gene_disease_associations.tsv.gz

gda <- read_tsv(gzfile("/workspaces/class-rnaseq/datasets_reference_only/trascriptome/all_gene_disease_associations.tsv.gz"))

disease2gene=gda[, c("diseaseId", "geneId")]
disease2name=gda[, c("diseaseId", "diseaseName")]

disgnet = enricher(entrez_genes_sig, TERM2GENE=disease2gene, TERM2NAME=disease2name)

## writeLines(summary(disgnet), "summary_disgnet.txt")

cnetplot(disgnet, foldChange=resdata$log2FoldChange[which(resdata$padj<0.5)])


save.image("deseq2_analysis.RData")