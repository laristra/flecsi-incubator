#!/bin/bash -l
#SBATCH -N 16
#SBATCH --ntasks=512
#SBATCH --ntasks-per-node=32
#SBATCH -t 2:00:00
#SBATCH -J hdf32
#SBATCH -p knl
export MPICH_MPIIO_HINTS="*:romio_cb_read=disable"
srun -N 16 --ntasks-per-node 32 ./hdf5proxy -size 412316860416 -nb_files 1
