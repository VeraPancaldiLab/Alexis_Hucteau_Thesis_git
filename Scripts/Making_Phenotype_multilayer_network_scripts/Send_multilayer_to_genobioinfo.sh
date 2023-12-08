#!/bin/bash

cd ~/work/Final_mutlilayer/

Pheno=$(ls)

Comp="Comparison"

for pheno in $Pheno
do
  if [[ -d $pheno ]] && [[ "$pheno" != *"$Comp"* ]] && [[ "$pheno" != *"Scaled"* ]];
  then
    cd $pheno
    file=$pheno\_multilayer_infomap.edges
    echo $file
    sshpass -p PASSWORD scp $file ahucteau@genobioinfo.toulouse.inrae.fr:work/Final_Multilayers/$pheno/
    cd ..
  fi
done
