#!/bin/bash


fqList=$1
ref=$2
thread=$3
if [[ $ref == "" ]]; then
	$ref="~/resources/references/Mus_musculus.GRCm38.cdna.all.index"
fi

if [[ $thread == "" ]]; then
	$thread = 4
fi

for f in `cat $fqList`
do g=`echo $f | cut -f1 -d '_'`

## run kallisto with 4 CPU threads
kallisto quant -i $ref -o kallisto.sample.$g --single -t $thread -l 180 -s 20 -b 100 $f
echo "$f is done"
done
