#!/bin/bash 
#SBATCH -J WI-10
#SBATCH -A m1266
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -n 64
#SBATCH -C cpu
#SBATCH -q debug
##SBATCH -p <your partition>
#SBATCH --time 10
#SBATCH -o Report-%j.out

export CP2K_DATA_DIR=$CONDA_PREFIX/cp2k-data
export OMP_NUM_THREADS=2
export ASE_CP2K_COMMAND="mpirun -np 64 $CONDA_PREFIX/bin/cp2k_shell.psmp"

ulimit -s unlimited

python sim_pysages_cp2k.py
