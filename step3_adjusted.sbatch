#!/bin/csh
#SBATCH -t 72:00:00
#SBATCH -J sbatch
#SBATCH --ntasks=1
#SBATCH --export=NONE
#SBATCH --partition=ord
#SBATCH --mem-per-cpu=128g
#SBATCH --gid=hab4cast
#SBATCH --account=hab4cast

cd $SLURM_SUBMIT_DIR
module load intel/19.0.3
module load gdal-2.4.1/intel-19.0
module load proj-5.2.0/intel-19.0
module load geos/3.7.1
export TMPDIR=/work/HAB4CAST/data_processing/data/OLCI_preprocessing/
module load R/3.5.2
R --no-save <cyanoCONUS_ice_step3_adjusted.R