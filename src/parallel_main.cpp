#include <iostream>
#include <fstream>
#include <sstream>
#include "Smooth.h"
#include <mpi.h>

int main(int argc, char **argv){
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  Smooth smooth(500,400,21,rank,size);
  smooth.SeedRing(300,200);
  smooth.SeedDisk(250,150);
  smooth.SeedDisk(225,125);
  std::ostringstream fname;
  fname << "test" << rank << ".dat" << std::flush;
  std::ofstream outfile(fname.str().c_str());
  outfile << smooth.LocalXSize() << ", " << smooth.Sizey() << ", " << rank << ", " << size << std::endl;
  std::cout << "Rank " << rank << "ready" << std::endl;
  while (smooth.Frame()<10) {
    smooth.Write(outfile);
    smooth.CommunicateMPI();
    smooth.QuickUpdate();
    std::cout << "Rank " << rank << "Completed frame: " << smooth.Frame() << std::endl;
  }
  MPI_Finalize();
}
