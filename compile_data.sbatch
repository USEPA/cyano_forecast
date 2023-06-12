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
module load geos/3.10.1
module load gdal-3.3.3/intel-21
module load proj-7.2.1/intel-21
module load udunits-2.2.28/intel-21
R --no-save <compile_data.R