#!/bin/bash
#SBATCH --mail-user=nathan.t.james@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:30:00
#SBATCH --mem=2G
#SBATCH --array=200-2500:100
#SBATCH --output=rstanarm_ord_time_logit_%a_levs.out

module restore R_basic
R --version

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

Rscript --no-save rstanarm_ord_time_logit_accre.r $SLURM_ARRAY_TASK_ID
