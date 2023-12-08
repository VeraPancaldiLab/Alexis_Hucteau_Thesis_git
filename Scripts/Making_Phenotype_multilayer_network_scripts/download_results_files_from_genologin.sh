#!/bin/bash

cd ~/GitHub/Thesis_paper/Results/Multi_layer_pheno_like/Final_mutlilayer/

Classes=$(ls)

Comp="Comparison"

for class in $Classes
do
  if [[ -d $class ]] && [[ "$class" != *"$Comp"* ]] && [[ "$class" != *"Scaled"* ]] && [[ "$class" == *"_"* ]] && [[ "$class" != *"esult"* ]] ;
  then
    echo $class
    echo ahucteau@genologin.toulouse.inra.fr:work/Multilayer/Final_Multilayers/$class/$class\_multilayer_infomap.tree
    sshpass -p PASSWORD scp ahucteau@genologin.toulouse.inra.fr:work/Final_Multilayers/$class/$class\_multilayer_infomap.tree ~/GitHub/Thesis_paper/Results/Multi_layer_pheno_like/Final_mutlilayer/$class/
  fi
done
