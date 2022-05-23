#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --time=1:00:00
#SBATCH --mem=150G
#SBATCH -p cpu_bioscienze
#SBATCH --mail-user=marco.malatesta@unipr.it
#SBATCH --mail-type=BEGIN,END,FAIL


while getopts l:d: flag
do
    case "${flag}" in
        l) ligand=${OPTARG};;
        d) directory=${OPTARG};;
    esac
done
# cd $directory

for file in ${directory}/*.trg;
do
   # wait here if the number of jobs is 3 (or more)
   while (( (( $(jobs -p | wc -l) )) >= 20)) 
   do 
      sleep 5      # check again after 5 seconds
   done

total_ade.py $file $ligand &
done
wait

