#!/bin/bash -l
#$ -cwd
# modified from generateSummaryTable.sh in star_code
# created by EMB
# executed in kallistoRelease.sh

outDir=$1
cd $outDir
echo "SampleID,Number_Reads,Pseudoaligned_Reads,Pseudoaligned_Ratio,Unique_Reads" > summary.csv
chmod 770 summary.csv

sample_dirs=`ls --hide="summary.csv"`

for i in $sample_dirs
do
    sampleID=`echo $i | cut -f1,3 -d '_'`
    numberReads=`grep "n_processed" $i/run_info.json | cut -f2 -d ":" | cut -c 2- | rev | cut -c 2- | rev`
    pseudoalignedReads=`grep "n_pseudoaligned" $i/run_info.json | cut -f2 -d ":" | cut -c 2- | rev | cut -c 2- | rev`
    pseudoalignedRatio=`grep "p_pseudoaligned" $i/run_info.json | cut -f2 -d ":" | cut -c 2- | rev | cut -c 2- | rev`
    uniqueReads=`grep "n_unique" $i/run_info.json | cut -f2 -d ":" | cut -c 2- | rev | cut -c 2- | rev`
    echo $sampleID", "$numberReads", "$pseudoalignedReads", "$pseudoalignedRatio", "$uniqueReads>>summary.csv

done

cd ..
