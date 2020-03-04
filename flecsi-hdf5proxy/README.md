#!/bin/sh

#set -v



# We need hdf5 built with mpi and C++, see info later for setting up compilers.yaml and packages.yaml on Sierra first (and add %gcc7.3.1 to spack install)

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


On Sierra: 
cat ~/.spack/linux/compilers.yaml
compilers:
- compiler:
    spec: gcc@7.3.1
    paths:
      cc: gcc
      cxx: g++
      f77: gfortran
      fc: gfortran
    flags: {}
    operating_system: rhel7
    target: ppc64le
    modules:
    - gcc/7.3.1
    environment: {}
    extra_rpaths: []



cat ~/.spack/linux/packages.yaml
packages:
 spectrum-mpi:
  modules:
   spectrum-mpi@rolling-release:  spectrum-mpi/rolling-release
 cmake:
  modules:
   cmake@3.14.5 : cmake/3.14.5
 python:
  modules:
   python@3.7.2 : python/3.7.2
