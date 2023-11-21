#!/bin/bash

cd ~/work/Final_Multilayers
Classes=$(ls -d */ | rev | cut -c2- | rev)

for class in $Classes
do
  if [[ -d $class ]] && [[ -z "$(ls -A ./$class)" ]];
  then
    script_file=$class\_run_infomap.sh
    output=/home/ahucteau/work/Final_Multilayers/$class
    echo $script_file
    echo $output
    touch $script_file
    echo '#!/bin/bash' > $script_file
    echo '#SBATCH -p workq' >> $script_file

    echo '#SBATCH -t 07:00:00' >> $script_file
    echo '#SBATCH --cpus-per-task 16' >> $script_file
    echo '#SBATCH --mem=100G' >> $script_file
    echo '#Load modules' >> $script_file
    echo '#Need gcc-11.3.0' >> $script_file
    echo 'module load compiler/gcc-11.3.0' >> $script_file
    echo 'module load bioinfo/infomap-v2.7.1' >> $script_file
    echo ' ' >> $script_file
    echo 'infomap /home/ahucteau/work/Final_Multilayers/'$class'_multilayer_infomap.edges \' >> $script_file
    echo '--seed 12345 -N 5 -M 16 -L 16 --inner-parallelization \' >> $script_file
    echo '-d -v -f directed '$output >> $script_file
    cat $script_file
    sbatch $script_file
    sleep 6h
  fi
done
