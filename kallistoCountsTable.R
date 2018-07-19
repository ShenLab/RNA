#!/bin/env Rscript
# created by EMB
# file name: kallistoCountsTable.R
# merges estimated counts of abundance.tsv files from all samples in a project directory
# successor to  multimergRNA.R for FeatureCounts data
# executed in kallistoRelease.sh 

WD=getwd()
dir = file.path(WD,"kallisto_raw_output")
dir

# merging transcript estimated counts
multimerge = function(mypath){
filenames=list.files(path=mypath, pattern = "*.tsv", full.names=TRUE)
datalist = lapply(filenames, function(x){read.table(file=x, header=T, col.names = c("target_id", "NULL", "NULL", sub('\\..*$', '', basename(x)), "NULL"), colClasses = c("character", "NULL", "NULL", "numeric", "NULL"), skipNul = TRUE)})
# round output counts to integer
for (x in 1:length(datalist)) {
  datalist[[x]][,2] = round(datalist[[x]][,2])
}
Reduce(function(x,y) {merge(x,y,all=TRUE,by="target_id",sort=FALSE)}, datalist)
}

# merging transcript tpms
multimerge_TPM = function(mypath){
  filenames=list.files(path=mypath, pattern = "*.tsv", full.names=TRUE)
  datalist = lapply(filenames, function(x){read.table(file=x, header=T, col.names = c("target_id", "NULL", "NULL", "NULL", sub('\\..*$', '', basename(x))), colClasses = c("character", "NULL", "NULL", "NULL", "numeric"), skipNul = TRUE)})
  Reduce(function(x,y) {merge(x,y,all=TRUE,by="target_id",sort=FALSE)}, datalist)
}


merged <- multimerge(dir)
merged_TPM <- multimerge_TPM(dir)
write.table(merged, "countstable.txt", sep="\t", quote=FALSE, row.names=FALSE, na="0")
write.table(merged_TPM, "tpmtable.txt", sep="\t", quote=FALSE, row.names=FALSE, na="0")

