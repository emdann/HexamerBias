#!/bin/bash
#
#$ -cwd

my_script=$1
id=$(echo $my_script | sed 's,.*/,,')
shift
threads=$2
shift
my_files="$@"

first_j=1
last_j="$#"

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

## Run job array
echo "/hpc/hub_oudenaarden/edann/bin/coverage_bias/utils/send_array_job.sh $my_script" | qsub -cwd -t $first_j-$last_j -N array_${id} -pe threaded 10 -l h_rt=24:00:00 -l h_vmem=80G -l h_cpu=2:00:00
