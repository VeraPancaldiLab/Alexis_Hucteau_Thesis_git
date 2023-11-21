#!/bin/bash

cd /mnt/SERVER-CRCT-STORAGE/CRCT18/Data\ Etudes\ Transcriptomic\ A\ NE\ PAS\ PARTAGER\ SANS\ AVIS\ JES/2023_U937_ATF4rep_AraC_D0_D2_D8_D12_JE508/Novaseq_250523/fastq/

Files=$(ls)

for fastq in $Files
do
  if [[ -f $fastq ]] && [[ "$fastq" != *"K562_111"* ]] && ([[ "$fastq" == *"K562"* ]] || [[ "$fastq" == *"MOLM13"* ]] || [[ "$fastq" == *"TUH"* ]]);
  then
    echo $fastq
    sshpass -p "135264/Qgmatp" scp $fastq ahucteau@genologin.toulouse.inra.fr:work/RNAseq/
  fi
done
