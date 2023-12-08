#!/bin/bash

cd ~/GitHub/Thesis_paper/Results/Multi_layer_pheno_like/Final_mutlilayer/

Pheno=$(ls)

Comp="Comparison"

for pheno in $Pheno
do
  if [[ -d $pheno ]] && [[ "$pheno" != *"$Comp"* ]] && [[ "$pheno" != *"Scaled"* ]] && [[ "$pheno" == *"_"* ]] ;
  then
    cd $pheno
    file=$pheno\_multilayer_infomap.edges
    echo $file
    sshpass -p PASSWORD scp $file ahucteau@genologin.toulouse.inra.fr:work/Final_Multilayers/$pheno/
    cd ..
  fi
done
