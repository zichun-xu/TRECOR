#!/bin/bash

#SBATCH --array=1-500
#SBATCH --nodes=1
#SBATCH --output=/home/dxu/Rout/par-%J.out
#SBATCH --error=/home/dxu/Rout/par-%J.err
#SBATCH --cpus-per-task=1

#SBATCH --mail-type=ALL
#SBATCH --mail-user=zx87@uw.edu

ml fhR/4.4.0-foss-2023b
R CMD BATCH --no-save /home/dxu/Research_projects/tree_cov_reg/scripts/simu_run_methods.R Rout/simu_slurm${SLURM_ARRAY_TASK_ID}.Rout