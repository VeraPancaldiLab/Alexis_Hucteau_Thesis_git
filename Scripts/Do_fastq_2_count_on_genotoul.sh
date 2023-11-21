#!/bin/bash

cd ~/work/RNAseq

nomenclature_1=$(echo _R1_001.fastq.gz)
nomenclature_2=$(echo _R2_001.fastq.gz)

nFiles=$(ls | wc -l)
nFiles=$(expr $nFiles / 2)
echo There are $nFiles samples to process

for fastq in *$nomenclature_1
do
  echo $fastq
  touch $fastq\_star_script.sh
  echo '#!/bin/bash' > $fastq\_star_script.sh
  echo '#SBATCH -p workq' >> $fastq\_star_script.sh

  echo '#SBATCH -t 4-00:00:00' >> $fastq\_star_script.sh
  echo '#SBATCH --cpus-per-task 8' >> $fastq\_star_script.sh
  echo '#SBATCH --mem-per-cpu=4000' >> $fastq\_star_script.sh

  echo '#Load modules' >> $fastq\_star_script.sh
  echo '#Need Miniconda3' >> $fastq\_star_script.sh
  echo 'module load system/Miniconda3' >> $fastq\_star_script.sh
  echo 'module load bioinfo/STAR-2.6.0c' >> $fastq\_star_script.sh
  echo 'STAR --quantMode GeneCounts --genomeDir ~/work/star-genome --runThreadN 8 \' >> $fastq\_star_script.sh
  echo '--readFilesIn '$fastq ${fastq%_R1_001.fastq.gz}_R2_001.fastq.gz '--readFilesCommand zcat --outFileNamePrefix '$fastq' \' >> $fastq\_star_script.sh
  echo '--outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate' >> $fastq\_star_script.sh
done
