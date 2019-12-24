#!/bin/bash -l
#SBATCH -N 1
#SBATCH -J hdf5test

dsize=$(( 10**6 ))
for ((inc=0; i<7; i++))
do
   runstring="srun -n 64 ./hdf5proxy -size $dsize -nb_files 4"
   echo $runstring
   echo -n "Data size is $dsize "
   eval "$runstring  |grep 'Elapsed time'"
   dsize=$(( $dsize*2 ))
done

