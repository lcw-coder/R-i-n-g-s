#!/bin/bash
#SBATCH -J f90
#SBATCH -o lmp.out
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 1000:00:00
#SBATCH -p parallel
#SBATCH --output=slurm-%x.%J.out

source /home/fwtop/.bashrc
module purge
module load Anaconda3/2022.05 
module load intel/2019b
module load Anaconda3/2022.05 
module load OpenMPI/4.0.3-GCC-9.3.0
unlimit -m ulimited
ulimit -s unlimited
ulimit -l unlimited

export OMP_NUM_THREADS=1
#export OMPI_MCA_btl_openib_allow_ib=1
#export OMPI_MCA_btl_openib_if_include="mlx5_0:1"

ifort ree.f90 -o cdft
./cdft

