#!/bin/bash -l

module load java/16.0.2 gcc/8.3.0 python3/3.8.10

for file in $1/bgoutFast.fa $1/fgoutFast.fa
do
  $2 -i ${file} -b ${file} -goterms -T $TMPDIR
done

module unload java/16.0.2 gcc/8.3.0 python3/3.8.10
