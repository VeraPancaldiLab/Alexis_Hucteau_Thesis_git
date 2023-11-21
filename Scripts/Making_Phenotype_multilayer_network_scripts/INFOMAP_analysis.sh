#! /bin/bash

read -p 'File to analyse: ' File
read -p 'output file: ' Output


date

echo ~/Test_Git/mux/muxViz/src-exe/infomap-0.x/Infomap $File --input-format multilayer \
--clu --map --tree --expanded \
--seed 12345 --num-trials 25 \
-u --two-level --inner-parallelization -vvv \
--out-name $Output

~/Test_Git/mux/muxViz/src-exe/infomap-0.x/Infomap $File $Output --input-format multilayer \
--clu --map --tree --expanded \
--seed 12345 --num-trials 25 \
-u --two-level --inner-parallelization -vvv \
--out-name $File

date
