#!/bin/bash
set -e

# Design parameters
# Nodes ==> 1, 16, 32, 64, and 96 nodes
# 32 proc per node and 64 proc per node
# 1 file for every 16 nodes
# Writing out 25% of memory per node
# Using striping factor of 8 -- see mpi_info setting in hdf5proxy 

# 1/4 of memory fails allocation -- changing to 1/32

module swap craype-haswell craype-mic-knl
module load cray-hdf5-parallel/1.10.5.2
module load cmake/3.16.2

module list

SCRATCH_DIR=/lustre/ttscratch1/gshipman
EXEC=/usr/projects/eap/users/gshipman/hdf5proxy/hdf5proxy
MEM_FRACTION=1/4

rm -rf ${SCRATCH_DIR}/hdf5proxy
mkdir ${SCRATCH_DIR}/hdf5proxy
cd ${SCRATCH_DIR}/hdf5proxy

echo "Memory per node is 96 GB. Writing out ${MEM_FRACTION} of node memory"
for ranks_per_node in 32 64
do
   echo ""
   for nodes in 1 16 32 64 96
   do
      dsizeGB=$(( nodes*96*${MEM_FRACTION} ))
      dsize=$(( nodes*96* 1024 * 1024 * 1024 *${MEM_FRACTION} ))
      #dsize=`echo '${nodes}*96*${MEM_FRACTION}*1024*1024*1024' | bc -l`
      nfiles=$(( nodes/16 ))
      ranks=$(( ${ranks_per_node}*${nodes} ))
      if [ ${nfiles} -eq 0 ]; then  nfiles=1; fi
      printf "%-38s" "Data size is $dsizeGB GiB Nodes is $nodes "
      runstring="srun -N $nodes --ntasks-per-node $ranks_per_node ${EXEC} -size $dsize -nb_files $nfiles"
      echo $runstring
      #echo -n "Data size is $dsizeGB GiB Nodes is $nodes	"
      #eval "$runstring  |grep 'Elapsed time'"
      BATCH_JOB="job${nodes}_${ranks_per_node}"
      echo "#!/bin/bash -l"                              >  $BATCH_JOB
      echo "#SBATCH -N $nodes"                           >> $BATCH_JOB
      echo "#SBATCH --ntasks=${ranks}"                   >> $BATCH_JOB
      echo "#SBATCH --ntasks-per-node=${ranks_per_node}" >> $BATCH_JOB
      echo "#SBATCH -t 1:00:00"                          >> $BATCH_JOB
      echo "#SBATCH -J hdf$nodes_$ranks_per_node"        >> $BATCH_JOB
      echo "#SBATCH -p knl"                              >> $BATCH_JOB
      echo "SCRATCH_DIR=${SCRATCH_DIR}"       >> $BATCH_JOB
      echo "EXEC=${EXEC}"      >> $BATCH_JOB
      echo "MEM_FRACTION=${MEM_FRACTION}"                           >> $BATCH_JOB
      echo "cd ${SCRATCH_DIR}/hdf5proxy"                 >> $BATCH_JOB
      echo "mkdir run${nodes}_${ranks_per_node}"         >> $BATCH_JOB
      echo "cd run${nodes}_${ranks_per_node}"            >> $BATCH_JOB
      #echo runstring="srun -N $nodes --ntasks-per-node $ranks_per_node ${EXEC} -size $dsize -nb_files $nfiles" >> $BATCH_JOB
      echo "echo \"$runstring \" "  >> $BATCH_JOB
      echo "$runstring  "            >> $BATCH_JOB
      echo $BATCH_JOB
      echo "==============="
      cat $BATCH_JOB
      echo ""
      sbatch < $BATCH_JOB
   done
done

echo ${SCRATCH_DIR}/hdf5proxy

#rm -rf ${SCRATCH_DIR}/hdf5proxy
