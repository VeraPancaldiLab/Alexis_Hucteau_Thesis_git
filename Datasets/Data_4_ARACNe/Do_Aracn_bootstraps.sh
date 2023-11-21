#! /bin/bash

read -p 'File to analyse: ' File
read -p 'TF list: ' TFs
#TFs=$(echo ~/GitHub/Multiplex_DNAmet_PPI_Chrom_Coexp/DATA/TF_names_v_1.01.txt)
read -p 'Output_file: ' output
read -p 'pvalue threshold (ex 1E-8): ' pvalue

java -Xmx5G -jar ~/ARACNe-AP/dist/aracne.jar -e $File \
-o $output \
--tfs $TFs \
--pvalue $pvalue --seed 1 --calculateThreshold

for i in {1..100}
do
  date
  echo $i
  java -Xmx5G -jar ~/ARACNe-AP/dist/aracne.jar \
  -e $File \
  -o $output \
  --tfs $TFs \
  --pvalue $pvalue --seed $i
done

java -Xmx5G -jar ~/ARACNe-AP/dist/aracne.jar -o $output --consolidate

# ~/shutdown_o_clock.sh
