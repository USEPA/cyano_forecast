#!/bin/csh
#SBATCH -t 72:00:00
#SBATCH -J sbatch
#SBATCH --ntasks=1
#SBATCH --export=NONE
#SBATCH --partition=ord
#SBATCH --mem-per-cpu=8g

cd $SLURM_SUBMIT_DIR
module load intel/21.4
module load R/4.1.2
R --no-save <generate_week_assignments_tibble.R