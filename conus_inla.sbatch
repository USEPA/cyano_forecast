#!/bin/csh
#SBATCH -t 72:00:00
#SBATCH -J sbatch
#SBATCH --ntasks=8
#SBATCH --export=NONE
#SBATCH --partition=largemem
#SBATCH --mem-per-cpu=64g
#SBATCH --gid=hab4cast
#SBATCH --account=hab4cast

cd $SLURM_SUBMIT_DIR
module load intel/21.4
module load R/4.2.0
module load gcc/9.1.0
module load geos/3.10.2
module load gdal-3.4.2/intel-21
module load proj-8.2.1/intel-21
module load udunits-2.2.28/intel-21
module load sqlite/3.38.3
R --no-save <conus_inla.R