#!/bin/env Rscript
# created by EMB
# Sleuth is memory intensive.  Use qrsh -l mem=40G
# must load Anaconda environment: source activate Sleuth
# execute above [organism]_kallisto directory
# include genome. Ex: "Rscript sleuth_DE_Analyses.R mouse"
# outputs DE_gene_list csv, Sleuth plots pdf

args <- commandArgs()
library("sleuth")

if (length(args)==0 | length(args)==1) {
  stop("Missing genome. Use Rscript sleuth_DE_Analyses.R [mouse|human]", call.=FALSE)
}

genome <- args[6]
print(args)
print(genome)

base_dir <- paste(genome, "kallisto", sep = "_") 
print(base_dir)


# manually insert samples to be compared in A vs. B groups: (ctrl is first)
A <- c("PA001","PA002", "PA003", "PA004")
B <- c("PA025", "PA026", "PA027", "PA028")
n <- length(A) # number of samples per group
samples <- c(A,B)


kallisto_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kallisto_dirs, sample=samples, condition = rep(c("ctrl", "exp"), each=n), stringsAsFactors=FALSE)

print(s2c)


if (genome == "human"){
  ensembl_dataset = "hsapiens_gene_ensembl"
}

if (genome == "mouse"){
  ensembl_dataset = "mmusculus_gene_ensembl"
}

if (genome != "human" & genome != "mouse") {
  stop("Genome has been entered incorrectly. Use Rscript sleuth_DE_Analysis.R [mouse|human]", call.=FALSE)
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

# create sleuth object, aggregating on gene-level, and do not filter out any genes
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g, aggregation_column = 'ext_gene', min_read = 0, min_prop = 0.5)
so <- sleuth_fit(so)
so <- sleuth_wt(so, 'conditionexp')
results_table <- sleuth_results(so, 'conditionexp', show_all = FALSE)

# now filter genes for plotting
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g, aggregation_column = 'ext_gene')
so <- sleuth_fit(so)
so <- sleuth_wt(so, 'conditionexp')

# write results to to_release
release_dir = file.path(base_dir, "to_release")
setwd(release_dir)

fnameStem <- paste( "Sleuth",
                    paste( sapply(samples, paste, collapse="-"),
                           collapse = "_" ),
                    sep="_" )

fname <- paste( fnameStem,
                "complete_DEG_list.csv",
                sep="-" )

results_ordered <- results_table[order(results_table$qval),]
write.csv(results_ordered , file=fname,row.names=F, quote=F)

fname <- paste( fnameStem,
                "qval_0.01_DEG_list.csv",
                sep="-" )


fname <- paste( fnameStem,
                "DEG_Plots.pdf",
                sep="-" )

pdf( file=fname )


top_30 <- results_ordered[1:30,1]

plot_ma(so, 'conditionexp', test_type = "wt", point_alpha = 0.6) + ggtitle("MA-Plot")
plot_volcano(so, 'conditionexp', test_type = "wt", point_alpha = 0.6) + ggtitle("Volcano Plot")
#plot_qq(so, 'conditionexp', test_type = "wt") + ggtitle("Q-Q Plot")
#plot_pca(so, color_by = 'condition',text_labels=TRUE) + ggtitle("PCA Plot")
#plot_pc_variance(so) + ggtitle("PCA Variances")

mat = sleuth:::spread_abundance_by(
  abund = so$obs_norm_filt,
  var = "scaled_reads_per_base",
  which_order = so$sample_to_covariates$sample)

# Calculate PCA, log2 transfer and normalization between samples
pc_x=1L
pc_y=2L
mat.pca <- prcomp(log2(t(mat)+1),center = TRUE,scale. = TRUE) 

#computation of variances
eigenvalues <- (mat.pca$sdev) ^ 2
var_explained <- eigenvalues * 100 / sum(eigenvalues)

#set label names
x_lab <- paste0('PC1 (', round(var_explained[1],digits=1))
x_lab <- paste0(x_lab, '%)')
y_lab <- paste0('PC2 (', round(var_explained[2],digits=1))
y_lab <- paste0(y_lab, '%)')

#Extract PC1 and PC2 to pcs, you can also change this PC1 and PC2 to others
pcs <- sleuth:::as_df(mat.pca$x[, c(pc_x, pc_y)])
pcs$sample <- rownames(pcs)
rownames(pcs) <- NULL

#add 'Group' information from experimental design
pcs <- dplyr::left_join(pcs, so$sample_to_covariates,by = 'sample')

#ggplot
pc_x <- paste0('PC', pc_x)
pc_y <- paste0('PC', pc_y)
ggplot(pcs, aes_string(pc_x, pc_y, colour = 'condition'))+
  geom_text(aes(label=sample),hjust=0, vjust=0)+
  ggtitle("PCA after Normalization (Log Transform)")+
  xlab(x_lab)+
  ylab(y_lab)+
  theme(text = element_text(size=10),axis.text = element_text(size=10),plot.title = element_text(size = 13)) 

plot_sample_heatmap(so) + ggtitle("Sample Heatmap")
plot_transcript_heatmap(so, transcripts = top_30) + ggtitle("Top 30 DE Genes Heatmap")
dev.off()

print("Sleuth is done.")

