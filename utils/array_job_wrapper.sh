#!/bin/bash
#
#$ -cwd

my_script=$1
my_files=$2

first_j=1
last_j=$(ls $my_files | wc -l)

## Make one directory for every job ##
n=$first_j;
max=$last_j;
while [ "$n" -le "$max" ]; do
  mkdir "dir.$n"
  n=`expr "$n" + 1`;
done

n=$first_j;
for file in $my_files; do
  mv $file dir.${n}
  n=`expr "$n" + 1`;
done

## Go to directory and run job

/hpc/hub_oudenaarden/edann/bin/coverage_bias/utils/send_array_job.sh $my_script
