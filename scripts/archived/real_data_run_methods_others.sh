#!/bin/bash

#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --output=/home/dxu/Rout/par-%J.out
#SBATCH --error=/home/dxu/Rout/par-%J.err
#SBATCH --cpus-per-task=1

#SBATCH --mail-type=ALL
#SBATCH --mail-user=zx87@uw.edu

ml fhR/4.4.0-foss-2023b
R CMD BATCH --no-save /home/dxu/Research_projects/tree_cov_reg/scripts/real_data_run_methods_others.R Rout/Yatsunenko_others_slurm${SLURM_ARRAY_TASK_ID}.Rout