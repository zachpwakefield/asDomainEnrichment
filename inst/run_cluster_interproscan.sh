#!/bin/bash -l

#$ -P
#$ -pe omp
#$ -l h_rt=10:00:00
#$ -N
#$ -m ea
#$ -M
#$ -j y
output_location=
interproscan_location=

module load java/16.0.2 gcc/8.3.0 python3/3.8.10

for file in $output_location/bgoutFast.fa $output_location/fgoutFast.fa
do
  $interproscan_location -i ${file} -b ${file} -goterms -T $TMPDIR
done

module unload java/16.0.2 gcc/8.3.0 python3/3.8.10
