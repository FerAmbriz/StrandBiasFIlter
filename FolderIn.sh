#!/bin/bash

d=($(ls -d */))
country=('COLOMBIA' 'GUATEMALA' 'PERU')

a=($(ls $d/*.txt))
head -1 $a | awk '{OFS="\t";print "'clinica'", $0}' | awk '{OFS="\t";print "'country'", $0}' | awk '{OFS="\t";print "'sample'", $0}'>> raw_data.txt

for i in ${d[@]};do
	prefix='HAPLOTYPECALLER_'; I=${i/#$prefix}; I=${I%/*}
	clinica=$I
	[[ ${country[@]} =~ (^|[[:space:]])$I($|[[:space:]]) ]] && country=$I || country='Mexico'
	#[[ ${country[@]} =~ (^|[[:space:]])$I($|[[:space:]]) ]] && clinica='NaN' || clinica=$I

	for j in $i/*.txt;do
		awk 'NR>1{OFS="\t";print "'$clinica'", $0}' $j | awk '{OFS="\t";print "'$country'", $0}' | awk '{OFS="\t";print "'$j'", $0}'>> raw_data.txt
		echo 'DONE'  $j
done
done
