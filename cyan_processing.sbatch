#!/bin/csh
#SBATCH -t 72:00:00
#SBATCH -J sbatch
#SBATCH --export=NONE
#SBATCH --partition=largemem
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=24GB
#SBATCH --gid=hab4cast
#SBATCH --account=hab4cast


cd $SLURM_SUBMIT_DIR
module load intel/21.4
module load R/4.1.2
module load gcc/6.1.0
module load geos/3.10.1
module load gdal-3.3.3/intel-21
module load proj-7.2.1/intel-21
module load udunits-2.2.28/intel-21
R --no-save <cyan_processing_conus.R