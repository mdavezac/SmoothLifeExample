#include <sstream>
#include <fstream>
#include "SmoothWriter.h"

void SmoothWriter::Write() {

  char * start_to_write=reinterpret_cast<char*>(smooth.StartOfWritingBlock());
  unsigned int element_count=smooth.LocalXSize()*smooth.Sizey();
  xdr_vector(&xdrfile,start_to_write,element_count,sizeof(double),reinterpret_cast<xdrproc_t>(xdr_double));
}

SmoothWriter::~SmoothWriter(){
}

SmoothWriter::SmoothWriter(Smooth & smooth, int rank, int size)
    :smooth(smooth),rank(rank),size(size)
{
     std::ostringstream fname;
     fname << "frames" << rank << ".dat" << std::flush;
     std::string mode("w");
     std::FILE * myFile = std::fopen(fname.str().c_str(),mode.c_str());
     xdrstdio_create(&xdrfile, myFile, XDR_ENCODE);
}

void SmoothWriter::Header(int frames){
  int sizey=smooth.Sizey();
  int sizex=smooth.LocalXSize();
  xdr_int(&xdrfile,&sizex);
  xdr_int(&xdrfile,&sizey);
  xdr_int(&xdrfile,&rank);
  xdr_int(&xdrfile,&size);
  xdr_int(&xdrfile,&frames);

}
