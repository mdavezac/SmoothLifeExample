#include <fstream>
#include <iostream>
#include "Smooth.h"

int main(int argc, char **argv) {
  Smooth smooth(200, 500, 21);
  smooth.AddRing();
  std::ofstream outfile("test.dat");
  outfile << smooth.Sizex() << ", " << smooth.Sizey() << std::endl;
  while(smooth.Frame() < 10) {
    smooth.Write(outfile);
    smooth.Update();
    std::cout << smooth.Frame() << std::endl;
  }
}
