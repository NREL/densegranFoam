#!/bin/bash
#SBATCH --ntasks=18
#SBATCH --time=4:00:00
#SBATCH --qos=high
#SBATCH --account=biorfeedflow  # Where to charge NREL Hours
#SBATCH --mail-user=syed.ahsan@nrel.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=logCons.%j.out

module purge
module load gcc
module load openmpi/3.1.3/intel-18.0.3
#source /home/sahsan/OpenFOAM/OpenFOAM-5.x/etc/bashrc

TMPDIR=/tmp/scratch/
#srun -n30 /home/sahsan/DEM_MFIX/MSDEM_NREL-master/AML_MSDEM_MPI_V1.0/src/mfixsolver -f mfix.dat
#application=`getApplication`
#blockMesh
#topoSet
#createBaffles
#setFields
#decomposePar
#$application > log &
./meshCreate
srun -n 18 fcicTwoPhaseEulerFoam -parallel 
