#!/bin/bash -l

#$ -P evolution
#$ -pe omp 28
#$ -l h_rt=10:00:00
#$ -N qsubtry
#$ -m ea
#$ -M zachpwakefield@gmail.com
#$ -j y
output_location='/projectnb2/evolution/zwakefield/proteinImpacts/mpc_cpc/'
interproscan_location='/projectnb2/evolution/zwakefield/proteinImpacts/my_interproscan/interproscan-5.65-97.0/interproscan.sh'

module load java/16.0.2 gcc/8.3.0 python3/3.8.10

for file in $output_location/bgoutFast.fa $output_location/fgoutFast.fa
do
  $interproscan_location -i ${file} -b ${file} -goterms -T $TMPDIR
done

module unload java/16.0.2 gcc/8.3.0 python3/3.8.10
