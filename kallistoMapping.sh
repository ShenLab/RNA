#!/bin/bash
#$ -cwd
# created by EMB
# File name: kallistoMapping.sh
 
USAGE="Usage: kallistoMapping.sh -s sample -g genome -l libtype -o output -f fastq1 -p fastq2"

KALLISTOBASE='/ifs/data/c2b2/ngs_lab/ngs/code/NGS/RNA_seq/KALLISTO'

while getopts "s:i:l:o:f:p:h" opt;
    do
    case "$opt" in
	s) sample="$OPTARG";;
        i) index="$OPTARG";;
        l) libtype="$OPTARG";;
        o) output="$OPTARG";;
        f) fastq1=("$OPTARG");;
        p) fastq2=("$OPTARG");;
	h) echo $USAGE; exit 1
    esac
done 
shift "$((OPTIND-1))"

# make output directory if it doesnt exist
if [[ ! -d $output ]]
then
    mkdir -p $output
fi


# check if SE (deprecated)
if [[ $fastq2 == "NA" ]]
   then
   if [[ $libtype == "N" ]]
      then 
      cmd="kallisto quant -i $index -o $output -b 100 -l 200 -s 20 --single $fastq1"
   else
      cmd="kallisto quant -i $index -o $output -b 100 -l 200 -s 20 --single $fastq1 --rf-stranded"
   fi
else   
   if [[ $libtype == "N" ]]
      then 
      cmd="kallisto quant -i $index -o $output -b 100 $fastq1 $fastq2"
   else
      cmd="kallisto quant -i $index -o $output -b 100 $fastq1 $fastq2 --rf-stranded"
   fi 
fi 

echo $fastq1
echo $fastq2
echo $cmd

$cmd


echo "Kallisto is done"
echo `date +"%H:%M:%S"`


