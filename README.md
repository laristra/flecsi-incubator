#!/bin/sh
data_sizes="100 400 800 1600"

#set -v

# We need hdf5 built with mpi and C++
spack install hdf5+mpi+cxx

# Load the hdf5 module
spack load -r hdf5

# Run the proxy with arguments for data size and number of files
for dsize in $data_sizes
do
   echo -n "Data size is $dsize "
   ./hdf5proxy -size $dsize -nb_files 1 |grep 'Elapsed time'
done
