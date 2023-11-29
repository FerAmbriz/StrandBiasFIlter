#!/bin/bash
folder=$1

for i in small*.txt;do
  ID=${i%*.txt}; ID=${ID##*/};prefix='small_file'; ID=${ID/#$prefix}
  python StrandBias_strelka_haplotype_varscan.py $i haplotype . $ID
done

echo '============ Done strandbias ============'

python merge_custom.py $folder

echo '============ Done merge ================='
