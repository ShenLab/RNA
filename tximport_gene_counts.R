#!/bin/env Rscript
# file name: tximport_gene_counts.R
# created by EMB
# creates gene counts table from kallisto files

library(tximport)
WD=getwd()
dir = file.path(WD,"kallisto_raw_output")
dir
tsv_files <- list.files(path=dir, pattern = "*.tsv", full.names=TRUE)
names(tsv_files) <- sub('\\..*$', '', basename(tsv_files)) 
genome <- sub(".*/(.*)_kallisto/to_release$","\\1", WD) 
genome

#check for valid genome name
if (genome == "human" | genome == "humanALL"){
  ensembl_dataset = "hsapiens_gene_ensembl"
} else if (genome == "mouse" | genome == "mouseALL"){
  ensembl_dataset = "mmusculus_gene_ensembl"
} else {
  stop("Genome has been entered incorrectly. [mouse|human]", call.=FALSE)
}


tx2g <- function(){
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = ensembl_dataset)
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                         "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                         ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  return(t2g)
}
t2g <- tx2g()

txi.kallisto <- tximport(tsv_files, type = "kallisto", tx2gene = t2g[,c(1,3)], ignoreTxVersion = TRUE)

gene_counts = round(txi.kallisto$counts)
gene_tpm = txi.kallisto$abundance
write.table( gene_counts, file="est_counts_genes_kallisto.txt", sep="\t", row.names=T, quote=F)
write.table( gene_tpm, file="tpm_values_genes_kallisto.txt", sep="\t", row.names=T, quote=F)
