#!/bin/bash
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH --job-name=mat_x_adm
#SBATCH --ntasks=2
#SBATCH --nodelist=slurm-node-1,

DIR=/home/admin/tutorial_examples/MatrixMultiplication_Flex-MPI

# exec with the selected number of processes
${DIR}/mat_x_adm.sh 2
