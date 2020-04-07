//========================================================================================
// (C) (or copyright) 2020. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001 for Los
// Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC
// for the U.S. Department of Energy/National Nuclear Security Administration. All rights
// in the program are reserved by Triad National Security, LLC, and the U.S. Department
// of Energy/National Nuclear Security Administration. The Government is granted for
// itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works, distribute copies to
// the public, perform publicly and display publicly, and to permit others to do so.
//========================================================================================

#include <mpi.h>
#include <hdf5.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <cstring>

#ifndef DEBUG
#define DEBUG 0
#endif

int RANKS_PER_FILE;
//long BUFFER_SIZE; // number of elements in double

bool create_hdf5_file(hid_t &hdf5_file_id, const std::string &file_name, MPI_Comm mpi_hdf5_comm)
{
  int rank;
  MPI_Comm_rank(mpi_hdf5_comm, &rank);

  hid_t file_creation_plist_id = H5P_DEFAULT;   // File creation property list
  hid_t file_access_plist_id   = H5P_DEFAULT;   // File access property list
  MPI_Info mpi_info  = MPI_INFO_NULL; // For MPI IO hints
  MPI_Info_create(&mpi_info);
  MPI_Info_set(mpi_info, "striping_factor", "8" );

  // Set up file access property list with parallel I/O access
  // H5Pcreate is a general property list create function
  // Here we are creating properties for file access
  file_access_plist_id = H5Pcreate(H5P_FILE_ACCESS);
  assert(file_access_plist_id >= 0);

  // Stores the MPI parameters -- comm, info -- in the property list
  int iret = H5Pset_fapl_mpio(file_access_plist_id, mpi_hdf5_comm, mpi_info);
  assert(iret >= 0);

  // Open the file collectively
  // H5F_ACC_TRUNC is for overwrite existing file if it exists. H5F_ACC_EXCL is no overwrite
  // 3rd argument is file creation property list. Using default here
  // 4th argument is the file access property list identifier
  hdf5_file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, file_creation_plist_id, file_access_plist_id);
  if (hdf5_file_id < 0) {
      if (rank == 0)
          std::cout << " H5Fcreate failed: " << hdf5_file_id << std::endl;
      return false;
  }
  if (DEBUG && rank == 0){
    std::cout << " create HDF5 file " << file_name
              //<< " file_id " << hdf5_file_id
              << std::endl;
  }

  // Terminates access to property list and frees all memory resources.
  iret = H5Pclose(file_access_plist_id);
  assert(iret >= 0);

  return true;
}

bool open_hdf5_file(hid_t &hdf5_file_id, const std::string &file_name, MPI_Comm mpi_hdf5_comm) {
  int rank;
  MPI_Comm_rank(mpi_hdf5_comm, &rank);

  hid_t file_access_plist_id   = H5P_DEFAULT;   // File access property list
  MPI_Info mpi_info  = MPI_INFO_NULL; // For MPI IO hints

  // Set up file access property list with parallel I/O access
  // H5Pcreate is a general property list create function
  // Here we are creating properties for file access
  file_access_plist_id = H5Pcreate(H5P_FILE_ACCESS);
  assert(file_access_plist_id >= 0);

  // Stores the MPI parameters -- comm, info -- in the property list
  int iret = H5Pset_fapl_mpio(file_access_plist_id, mpi_hdf5_comm, mpi_info);
  assert(iret >= 0);

  hdf5_file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDWR, file_access_plist_id);
  if (hdf5_file_id < 0) {
    if (rank == 0)
        std::cout << " H5Fopen failed: " << std::endl;
    return false;
  }
  if (DEBUG && rank == 0){
     std::cout << " open HDF5 file " << file_name
               //<< " file_id " << hdf5_file_id
               << std::endl;
  }
  return true;
}

bool close_hdf5_file(hid_t &hdf5_file_id, MPI_Comm mpi_hdf5_comm) {
  int rank;
  MPI_Comm_rank(mpi_hdf5_comm, &rank);

  assert(hdf5_file_id >= 0);
  herr_t status;
  status = H5Fflush(hdf5_file_id, H5F_SCOPE_LOCAL);
  assert(status >= 0);
  status = H5Fclose(hdf5_file_id);
  assert(status >= 0);
  if (DEBUG && rank == 0){
     std::cout << " close HDF5 file_id " << hdf5_file_id
               << std::endl;
  }
  hdf5_file_id = -1;
  return true;
}

bool create_hdf5_dataset(const hid_t hdf5_file_id, const std::string dataset_name, long ncount, MPI_Comm mpi_hdf5_comm)
{
  int rank;
  MPI_Comm_rank(mpi_hdf5_comm, &rank);

  const int ndims = 1;
  hsize_t dims[ndims];
  dims[0] = ncount;

  // 1st argument -- number of dimensions
  // 2nd argument -- array of current dimensions
  // 3rd argument -- maximum number of dimensions. NULL means that current is maximum
  // returns the dataspace id. Fortran interface has this inserted as arg 3 and arg 4 as err(max moves to arg 5).
  if (DEBUG  && rank == 0) std::cout << "creating dataset with dim " << ncount << std::endl;
  hid_t file_dataspace_id = H5Screate_simple(ndims, dims, NULL);

  // Creates a new dataset and links to file
  hid_t link_creation_plist_id    = H5P_DEFAULT; // Link creation property list
  hid_t dataset_creation_plist_id = H5P_DEFAULT; // Dataset creation property list
  hid_t dataset_access_plist_id   = H5P_DEFAULT; // Dataset access property list
  hid_t dataset_id = H5Dcreate2(hdf5_file_id,       // Arg 1: location identifier
    dataset_name.c_str(),                     // Arg 2: dataset name
    H5T_IEEE_F64LE,                           // Arg 3: datatype identifier
    file_dataspace_id,                        // Arg 4: dataspace identifier
    link_creation_plist_id,                   // Arg 5: link creation property list
    dataset_creation_plist_id,                // Arg 6: dataset creation property list
    dataset_access_plist_id);                 // Arg 7: dataset access property list

  if(dataset_id < 0) {
    if (rank == 0)
        std::cout << " H5Dcreate2 failed: " << dataset_id << std::endl;
    H5Sclose(file_dataspace_id);
    H5Fclose(hdf5_file_id);
    return false;
  }
  herr_t status;
  status = H5Dclose(dataset_id);
  assert(status == 0);
  status = H5Sclose(file_dataspace_id);
  assert(status == 0);
  status = H5Fflush(hdf5_file_id, H5F_SCOPE_LOCAL);
  assert(status == 0);
  return true;
}

bool write_data_to_hdf5(const hid_t &hdf5_file_id, const std::string dataset_name, const double *buffer, long ncount, long displs, MPI_Comm mpi_hdf5_comm) {
  herr_t status;
  int rank;
  MPI_Comm_rank(mpi_hdf5_comm, &rank);

  hid_t data_access_plist_id = H5P_DEFAULT;
  hid_t dataset_id = H5Dopen2(hdf5_file_id, dataset_name.c_str(), data_access_plist_id);
  if (dataset_id < 0) {
     if (rank == 0)
         std::cout << " H5Dopen2 failed: " << dataset_id << std::endl;
     H5Fclose(hdf5_file_id);
     return false;
  }

  const int ndims = 1;
  hsize_t count[1];
  hsize_t offset[1];
  count[0] = ncount;
  offset[0] = displs;
  hid_t mem_dataspace_id = H5Screate_simple(ndims, count, NULL);
  assert(mem_dataspace_id >= 0);

  /*
   * Select hyperslab in the file.
   */
  hid_t file_dataspace_id = H5Dget_space(dataset_id);
  assert(file_dataspace_id >= 0);
  if (DEBUG && rank == 0) std::cout << "Rank: "<< rank << " offset: " << displs << " count: " << ncount << std::endl;
  H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

  // Create property list for collective dataset write.
  hid_t xfer_plist_id = H5Pcreate(H5P_DATASET_XFER);
  assert(xfer_plist_id >= 0);

  status = H5Pset_dxpl_mpio(xfer_plist_id, H5FD_MPIO_COLLECTIVE);
  assert(status >= 0);

  // To write dataset independently use
  //    H5Pset_dxpl_mpio(xfer_plist_id, H5FD_MPIO_INDEPENDENT);
    
  status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, mem_dataspace_id, file_dataspace_id,
                    xfer_plist_id, buffer);
  assert(status >= 0);
  status = H5Sclose(mem_dataspace_id);
  assert(status >= 0);
  status = H5Sclose(file_dataspace_id);
  assert(status >= 0);
  status = H5Pclose(xfer_plist_id);
  assert(status >= 0);
  status = H5Dclose(dataset_id);
  assert(status >= 0);
  return true;
}

bool read_data_from_hdf5(const hid_t &hdf5_file_id, const std::string dataset_name, double *buffer, long ncount, long displs, MPI_Comm mpi_hdf5_comm) {
  herr_t status;
  int rank;
  MPI_Comm_rank(mpi_hdf5_comm, &rank);

  hid_t data_access_plist_id = H5P_DEFAULT;
  hid_t dataset_id = H5Dopen2(hdf5_file_id, dataset_name.c_str(), data_access_plist_id);
  if (dataset_id < 0) {
      if (rank == 0)
          std::cout << " H5Dopen2 failed: " << dataset_id << std::endl;
      H5Fclose(hdf5_file_id);
      return false;
  }

  const int ndims = 1;
  hsize_t count[1];
  hsize_t offset[1];
  count[0] = ncount;
  offset[0] = displs;
  hid_t mem_dataspace_id = H5Screate_simple(ndims, count, NULL);
  assert(mem_dataspace_id >= 0);

  /*
   * Select hyperslab in the file.
   */
  hid_t file_dataspace_id = H5Dget_space(dataset_id);
  assert(file_dataspace_id >= 0);

  status = H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
  assert(status >= 0);

  // Create property list for collective dataset write.
  hid_t xfer_plist_id = H5Pcreate(H5P_DATASET_XFER);
  assert(xfer_plist_id >= 0);

  status = H5Pset_dxpl_mpio(xfer_plist_id, H5FD_MPIO_COLLECTIVE);
  assert(status >= 0);
    
  status = H5Dread(dataset_id, H5T_IEEE_F64LE, mem_dataspace_id, file_dataspace_id,
                   xfer_plist_id, buffer);
  assert(status >= 0);
  status = H5Sclose(mem_dataspace_id);
  assert(status >= 0);
  status = H5Sclose(file_dataspace_id);
  assert(status >= 0);
  status = H5Pclose(xfer_plist_id);
  assert(status >= 0);
  H5Dclose(dataset_id);
  assert(status >= 0);
  return true;
}

int main(int argc, char** argv) {
  
  RANKS_PER_FILE = 4;
  long BUFFER_SIZE = 400;
  int nb_files = 1;
  
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-size")) {
      BUFFER_SIZE = atol(argv[++i]);
    }
    
    if (!strcmp(argv[i], "-nb_files")) {
      nb_files = atoi(argv[++i]);
    }
  }

  MPI_Init(&argc, &argv);

  int world_size, rank, new_world_size, new_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  BUFFER_SIZE = (BUFFER_SIZE / world_size) * world_size;
  nb_files = world_size / (world_size / nb_files);
  if (DEBUG && rank == 0) std::cout << "BUFFER_SIZE: " << BUFFER_SIZE << std::endl;
  if (DEBUG && rank == 0) std::cout << "#Files: " << nb_files << std::endl;
  assert(BUFFER_SIZE % world_size == 0);
  assert(world_size % nb_files == 0);
  RANKS_PER_FILE = world_size / nb_files;

  MPI_Comm new_comm;
  int new_color = rank / RANKS_PER_FILE;
  int nb_new_comms = world_size / RANKS_PER_FILE;
  MPI_Comm_split(MPI_COMM_WORLD, new_color, rank, &new_comm);

  MPI_Comm mpi_hdf5_comm  = new_comm;

  MPI_Comm_size(new_comm, &new_world_size);
  MPI_Comm_rank(new_comm, &new_rank);

  hid_t hdf5_file_id = -1;
  bool return_val = false;

  // initialize HDF5 library
  // this is only really required for the Fortran interface
  return_val = H5open();
  assert(return_val >= 0);

  // create hdf5 file
  if (DEBUG && rank == 0) std::cout << "Creating HDF5 file " << std::endl << world_size;

  std::string file_name = "checkpoint_" + std::to_string(new_color);
  return_val = create_hdf5_file(hdf5_file_id, file_name, mpi_hdf5_comm);
  assert(return_val == true);


#if 0
  long long begin = BUFFER_SIZE/nb_new_comms * (new_rank    ) / new_world_size;
  long long end   = BUFFER_SIZE/nb_new_comms * (new_rank + 1) / new_world_size;
  long long nsize = end - begin;

  long long nsize_global[new_world_size];
  int iret = MPI_Allgather(&nsize, 1, MPI_INT, nsize_global, 1, MPI_INT, new_comm);

  long long displs[new_world_size];
  displs[0] = 0;
  for (int i = 1; i < new_world_size; i++){
     displs[i] = displs[i-1] + nsize_global[i-1];
  }
#else
  long nsize = BUFFER_SIZE/world_size;
  long ncount = nsize / sizeof(double); 
  long displs[new_world_size];
  displs[new_rank] = ncount * new_rank;
#endif
  //if (rank == 0) std::cout << " hello! " << std::endl;
  if (DEBUG && rank == 0) std::cout << "create dataset " << ncount << " " << nb_new_comms << std::endl;
  return_val = create_hdf5_dataset(hdf5_file_id, "test_dataset", ncount/nb_new_comms*world_size, mpi_hdf5_comm);
  assert(return_val == true);

  return_val = close_hdf5_file(hdf5_file_id, mpi_hdf5_comm);
  assert(return_val == true);

  // Initialize data buffer
  double *buffer_checkpoint = NULL;
  buffer_checkpoint = (double*)malloc(nsize);
  assert(buffer_checkpoint != NULL);

  double *buffer_recover = NULL;
  buffer_recover = (double*)malloc(nsize);
  assert(buffer_recover != NULL);

  for (long i = 0; i < ncount; i++) {
    buffer_checkpoint[i] = 1.0 + (double)i + ncount*rank;
    buffer_recover[i] = 0.0;
  }

  double t1, t2, t3;
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    t1 = MPI_Wtime();
  }

  // checkpoint
  if (rank == 0) std::cout << "Writing checkpoint" << std::endl;
  return_val = open_hdf5_file(hdf5_file_id, file_name, mpi_hdf5_comm);
  assert(return_val == true);
  return_val = write_data_to_hdf5(hdf5_file_id, "test_dataset", buffer_checkpoint, ncount, displs[new_rank], mpi_hdf5_comm);
  assert(return_val == true);
  return_val = close_hdf5_file(hdf5_file_id, mpi_hdf5_comm);
  assert(return_val == true);

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    t2 = MPI_Wtime();
    std::cout << "Time to write checkpoint is " << t2 - t1 << " secs. "
              << "Bytes written " << BUFFER_SIZE << " bytes per file. "
              << "Write rate is " << BUFFER_SIZE/1024.0/1024.0/(t2 - t1) << " MiBs/sec per file" 
              << std::endl;
  }

  // recover
  if (rank == 0) std::cout << "Recovering  checkpoint" << std::endl;
  return_val = open_hdf5_file(hdf5_file_id, file_name, mpi_hdf5_comm);
  assert(return_val == true);
  return_val = read_data_from_hdf5(hdf5_file_id, "test_dataset", buffer_recover, ncount, displs[new_rank], mpi_hdf5_comm);
  assert(return_val == true);
  return_val = close_hdf5_file(hdf5_file_id, mpi_hdf5_comm);
  assert(return_val == true);
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    t3 = MPI_Wtime();
    std::cout << "Time to recover checkpoint is " << t3 - t2 << " secs. "
              << "Bytes read " << BUFFER_SIZE << " bytes per file. "
              << "Read rate is " << BUFFER_SIZE/1024.0/1024.0/(t3 - t2) << " MiBs/sec per file"
              << std::endl;
  }

  if (rank == 0) {
    t3 = MPI_Wtime();
    std::cout << "Elapsed time is " << t3 - t1 << std::endl;
  }

  if (rank == 0) std::cout << "Verifying  checkpoint" << std::endl;
  int ierr = 0;
  // verification
  for (long i = 0; i < ncount; i++) {
    assert(buffer_checkpoint[i] == buffer_recover[i]);
    if (buffer_checkpoint[i] != buffer_recover[i]) ierr++;
  }
  int ierr_global = 0;
  MPI_Allreduce(&ierr, &ierr_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0 && ierr_global == 0) std::cout << "Checkpoint has been verified" << std::endl;

  free(buffer_checkpoint);
  free(buffer_recover);

  MPI_Finalize();
}
