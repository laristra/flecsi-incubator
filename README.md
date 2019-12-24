#!/bin/sh

#set -v

# We need hdf5 built with mpi and C++
spack install hdf5+mpi+cxx

# Load the hdf5 module
spack load -r hdf5

# Compile the code
make

# Run the proxy with arguments for data size and number of files
dsize=$(( 10**6 ))
for ((inc=0; i<18; i++))
do
   runstring="mpirun -n 64 ./hdf5proxy -size $dsize -nb_files 4"
   echo $runstring
   echo -n "Data size is $dsize "
   eval "$runstring  |grep 'Elapsed time'"
   dsize=$(( $dsize*2 ))
done
