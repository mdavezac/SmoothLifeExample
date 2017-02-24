#include <string>
#include <vector>

typedef double density;
typedef double distance;
typedef double filling;

class Smooth {
public:
  Smooth(int sizex = 100, int sizey = 100, distance inner = 21.0, filling birth_1 = 0.278,
         filling birth_2 = 0.365, filling death_1 = 0.267, filling death_2 = 0.445,
         filling smoothing_disk = 0.147, filling smoothing_ring = 0.028);
  int Size() const;
  int Sizex() const;
  int Sizey() const;
  int Range() const;
  int Frame() const;
  const std::vector<density> &Field() const;

  //! \brief Piecewise linear function defining the disk
  //! \details
  //! - 1 inside the disk
  //! - 0 outside the disk
  //! - 0 < x < 1 in the smoothing region
  double Disk(distance radius) const;
  //! \brief Piecewise linear function defining the ring
  //! \details
  //! - 1 inside the ring
  //! - 0 outside the ring
  //! - 0 < x < 1 in the smoothing region
  double Ring(distance radius) const;
  //! Smooth step function: 0 at -∞, 1 at +∞
  static double Sigmoid(double variable, double center, double width);
  //! $e^{-4x / width}$: 0 at -∞, 1 at +∞
  static double Sigmoid(double x, double width);
  density Transition(filling disk, filling ring) const;
  int TorusDifference(int x1, int x2, int size) const;
  double Radius(int x1, int y1, int x2, int y2) const;
  // Value of the integral over a single ring
  double NormalisationRing() const;
  // Value of the integral over a single disk
  double NormalisationDisk() const;
  //! Sets the playing field to random values
  void SeedRandom();
  //! Sets the playing field to constant values
  void SeedConstant(density constant = 0);
  //! Adds a disk to the playing field
  void AddDisk(int x0 = 0, int y0 = 0);
  //! Adds a ring to the playing field
  void AddRing(int x0 = 0, int y0 = 0);
  //! Sets a single pixel in the field
  void AddPixel(int x0, int y0, density value);
  //! Moves to next step
  void Update();
  //! Prints current field to standard output
  void Write(std::ostream &out);

  std::pair<density, density> Integrals(int x, int y) const;

  //! index of element at position i, j
  int Index(int i, int j) const;

private:
  int sizex;
  int sizey;
  std::vector<density> field;
  std::vector<density> work_field;
  distance inner;
  filling birth_1;
  filling birth_2;
  filling death_1;
  filling death_2;
  filling smoothing_disk;
  filling smoothing_ring;
  distance outer;
  distance smoothing;
  int frame;
  double normalisation_disk;
  double normalisation_ring;
};
