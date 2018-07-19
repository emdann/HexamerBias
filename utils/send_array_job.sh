#!/bin/bash
#
#$ -cwd

my_script=$1

## Go to directory and run job
cd dir.$SGE_TASK_ID
input_file=$(ls)

$my_script $input_file
