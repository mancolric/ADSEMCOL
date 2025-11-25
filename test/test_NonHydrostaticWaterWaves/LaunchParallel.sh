#!/bin/bash

# https://guiesbibtic.upf.edu/recerca/hpc/array-jobs

# NOTE: SBATCH parameters must be defined first. Do not define any variable before them.

# Editable options:
#SBATCH --job-name=Soliton1
#SBATCH --partition=cn1
#SBATCH --array=1-6%3
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G

# Do not modify:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 
#SBATCH --threads-per-core=1
# Note that these data are per job in the job array

# Outputs: %x is the job name, %a is the array id number:
# (https://stackoverflow.com/questions/50242293/using-sbatch-job-name-as-a-variable-in-file-output)
#SBATCH --output=../../temp/%x-%a.out
#SBATCH --error=../../temp/%x-%a.err

# Start:
echo "Starting on $(date)"
echo "Solving at $(hostname)"

# Load modules:
module load julia
module load bamg

# Read line $SLURM_TASK_ID of $SLURM_JOB_NAME.txt file to get input arguments:
ARGS=$( awk 'NR=='$SLURM_ARRAY_TASK_ID $SLURM_JOB_NAME.lp )

#Call Julia:
srun --exclusive julia -t$SLURM_CPUS_PER_TASK Launch.jl $ARGS

#Finish:
echo "Finished on $(date)"

