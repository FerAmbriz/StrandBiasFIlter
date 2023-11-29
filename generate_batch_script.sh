#!/bin/bash

mkdir PNG
mkdir scripts

IGV="/home/ferambriz/Programs/IGV_Linux_2.16.0/igv.sh"

region=('chr15:89,786,965-89,787,374' 'chr16:23,652,741-23,652,948' 'chr2:47,596,797-47,596,968' 'chr2:47,630,194-47,630,325')
genes=('FANCI' 'PALB2' 'EPCAM' 'MSH2')
site=('chr15:89786999-89786999' 'chr16:23652916-23652916' 'chr2:47596828-47596828' 'chr2:47630224-47630224')

len_s=${#region[@]}

for j in /run/media/ferambriz/MIGUEL\ RUIZ/Bismark/norm/preprocessing/*.bam; do
  ID=${j%.*}; ID=${ID##*/};
  echo $ID
  for (( i=0; i<${len_s}; i=i+1 )); do
    echo ${region[i]} ${genes[i]} ${site[i]}
cat << EOF > scripts/igv_${ID}_${site[i]}_${genes[i]}.txt
new
genome hg19
load '$j'
snapshotDirectory /home/ferambriz/Documents/IGVscreen/PNG
goto ${region[i]}
region ${site[i]}
sort base
snapshot ${ID}_${site[i]}_${genes[i]}.png
exit
EOF
  done
done

for j in scripts/*.txt;
do
  $IGV -b $j
done 

