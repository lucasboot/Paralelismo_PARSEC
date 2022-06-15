#!/bin/bash 
#SBATCH --job-name=streamcluster
#SBATCH --time=2-00:00
#SBATCH --cpus-per-task=32
#SBATCH --hint=compute_bound
#SBATCH --partition=cluster,intel-512 --exclude=r1i0n[0-2,6-7,9-12,14-16],r1i1n[0-7,9-13,15-16],r1i2n[0-7,9-16],r1i3n[0-7,9-14]
#SBATCH --exclusive
export OMP_NUM_THREADS=32

pascalanalyzer -c 1:32 -i "10 20 128 1000000 200000 5000 none output_nartive.txt __nt__","10 20 128 16384 16384 1000 none output_large.txt __nt__"," 2 5 1 10 10 5 none output_test.txt __nt__" ./streamcluster -o main_noexclusive.json
