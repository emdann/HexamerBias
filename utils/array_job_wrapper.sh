#!/bin/bash
#
#$ -cwd

my_script=$1
my_files=$2

## Make one directory for every job ##
n=$SGE_TASK_FIRST;
max=$SGE_TASK_LAST;
while [ "$n" -le "$max" ]; do
  mkdir "dir.$n"
  n=`expr "$n" + 1`;
done

n=$SGE_TASK_FIRST;
for file in $my_files; do
  mv $file dir.${n}
  n=`expr "$n" + 1`;
done

## Go to directory and run job
cd dir.$SGE_TASK_ID
input_file=$(ls)

$my_script $input_file
