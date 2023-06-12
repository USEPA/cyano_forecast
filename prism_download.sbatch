#!/bin/csh
#SBATCH -t 72:00:00
#SBATCH -J sbatch
#SBATCH --ntasks=1
#SBATCH --export=NONE
#SBATCH --partition=ord
#SBATCH --mem-per-cpu=128g

cd $SLURM_SUBMIT_DIR
module load intel/21.4
module load R/4.1.2
module load geos/3.10.2
module load gdal-3.4.2/intel-21
module load proj-8.2.1/intel-21
R --no-save <prism_download.R