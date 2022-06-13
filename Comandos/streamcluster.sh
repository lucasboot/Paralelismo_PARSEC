#!/bin/bash 
#SBATCH --job-name=OMP_hello 
#SBATCH --time=2-00:00
#SBATCH --cpus-per-task=8
#SBATCH --hint=compute_bound
export OMP_NUM_THREADS=32

pascalanalyzer -c 1:32 -i "10 20 128 1000000 200000 5000 none output.txt 32","10 20 128 16384 16384 1000 none output.txt 32","3 10 3 16 16 10 none output.txt 32" ./streamcluster -o main.json
