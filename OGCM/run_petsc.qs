#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name NA
#SBATCH --nodes=16 --ntasks-per-node=40
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1

export OMP_NUM_THREADS=1

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_STACKSIZE=1000m
export LD_LIBRARY_PATH=/data4t/gusev/local/lib/:${LD_LIBRARY_PATH}
ulimit -s unlimited

mpirun ./inmsom > E/e.05 2>err.txt -ksp_type gmres -ksp_atol 1.0e-16 -ksp_rtol 1.0e-08 -pc_type asm -pc_asm_overlap 2 -sub_pc_type ilu -sub_pc_factor_levels 8

