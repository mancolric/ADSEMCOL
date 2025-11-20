#!/bin/bash

# Before launching this file, change directory to ADSEMCOL folder.

# See https://login.scg.stanford.edu/faqs/cores/
# See https://stackoverflow.com/questions/65603381/slurm-nodes-tasks-cores-and-cpus

# Job and output names:
#SBATCH --job-name=Soliton
#SBATCH --output=temp/Soliton.out

# Select partition. For SLURM, a node is a computer part of a cluster. We use the node cn1 or the node cn2:
#SBATCH --nodes=1
#SBATCH --partition=cn1

# There is only one task, which is executing LaunchParallel.jl. This task will be executed using W workers (i.e., W parallel tasks for SLURM), and each worker (parallel task) requires T threads/cores.
# That is, in the node cn1 or cn2, at most "ntasks-per-node" problems will be solved in parallel. Each problem (task) requires "cpus-per-task" CPUs.
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=4

# Load default libraries for julia:
module load julia

# Run julia with desired parameters:
julia test/test_NonHydrostaticWaterWaves/LaunchParallel.jl


