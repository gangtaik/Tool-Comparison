#!/bin/bash

for i in $(cat ../hg19/ll.txt)
do
	more ${i}_fa.txt |sed -e '1d' > mid.txt
	more mid.txt |sed '1d' > end.txt
	more end.txt |awk -F "\t" '{OFS="\t"}{print $3":"$4"-"$5}' > ${i}_fa_hg38.txt

done

