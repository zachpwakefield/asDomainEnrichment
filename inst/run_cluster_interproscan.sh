#!/bin/bash -l

#$ -P $1
#$ -pe omp $2
#$ -l h_rt=10:00:00
#$ -N $3
#$ -m ea
#$ -M $4
#$ -j y


module load java/16.0.2 gcc/8.3.0 python3/3.8.10

for file in $5/bgoutFast.fa $5/fgoutFast.fa
do
  $6 -i ${file} -b ${file} -goterms -T $TMPDIR
done

module unload java/16.0.2 gcc/8.3.0 python3/3.8.10
