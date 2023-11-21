#!/bin/bash -l

out_directory='/projectnb2/evolution/zwakefield/proteinImpacts/mpc_cpc/'
## SCC Settings
module load java/16.0.2 gcc/8.3.0 python3/3.8.10

for file in $out_directory/bgoutFast.fa $out_directory/fgoutFast.fa
do
  # ./my_interproscan/interproscan-5.65-97.0/interproscan.sh -i ${file} -f tsv -goterms -T $TMPDIR
  /projectnb2/evolution/zwakefield/proteinImpacts/my_interproscan/interproscan-5.65-97.0/interproscan.sh -i ${file} -b ${file} -goterms -T $TMPDIR
done

module unload java/16.0.2 gcc/8.3.0 python3/3.8.10
