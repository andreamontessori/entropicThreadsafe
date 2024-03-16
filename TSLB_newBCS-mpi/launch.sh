#!/bin/bash
#SBATCH --job-name rectTSLB
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=4     # Number of MPI ranks per node
#SBATCH --gres=gpu:4            # Number of requested gpus per node, can vary between 1 and 4
#SBATCH --time 01:00:00         # Walltime, format: HH:MM:SS
#SBATCH --mem=300000          # memory per node out of 494000MB
#SBATCH -A IscrC_recTSLB        # Name of the account
#SBATCH -p boost_usr_prod
#######SBATCH --qos=boost_qos_dbg
module load nvhpc/23.1
module load openmpi/4.1.4--nvhpc--23.1-cuda-11.8
#mpirun --mca orte_base_help_aggregate 0 --mca btl_base_warn_component_unused 0 -np 4 ./main.x 1 1 4 
mpirun -np 1 ./main.x 1 1 1
# 1 MPI task, 1 GPU
