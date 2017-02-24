#include <cmath>
#include <cstdlib>
#include <iostream>

#include "smooth.h"

Smooth::Smooth(int sizex, int sizey, distance inner, filling birth_1, filling birth_2,
               filling death_1, filling death_2, filling smoothing_disk, filling smoothing_ring)
    : sizex(sizex), sizey(sizey), field(sizex * sizey), work_field(sizex * sizey), inner(inner),
      birth_1(birth_1), birth_2(birth_2), death_1(death_1), death_2(death_2),
      smoothing_disk(smoothing_disk), smoothing_ring(smoothing_ring), outer(inner * 3),
      smoothing(1.0) {
  normalisation_disk = NormalisationDisk();
  normalisation_ring = NormalisationRing();
}

int Smooth::Range() const { return outer + smoothing / 2; }

int Smooth::Sizex() const { return sizex; }
int Smooth::Sizey() const { return sizey; }
int Smooth::Size() const { return sizex * sizey; }

double Smooth::Disk(distance radius) const {
  if(radius > inner + smoothing / 2)
    return 0.0;
  if(radius < inner - smoothing / 2)
    return 1.0;
  return (inner + smoothing / 2 - radius) / smoothing;
}

double Smooth::Ring(distance radius) const {
  if(radius < inner - smoothing / 2)
    return 0.0;
  if(radius < inner + smoothing / 2)
    return (radius + smoothing / 2 - inner) / smoothing;
  if(radius < outer - smoothing / 2)
    return 1.0;
  if(radius < outer + smoothing / 2)
    return (outer + smoothing / 2 - radius) / smoothing;
  return 0.0;
}

double Smooth::Sigmoid(double variable, double center, double width) {
  return Sigmoid(variable - center, width);
}
double Smooth::Sigmoid(double x, double width) { return 1.0 / (1.0 + std::exp(-4.0 * x / width)); }

density Smooth::Transition(filling disk, filling ring) const {
  auto const sdisk = Sigmoid(disk - 0.5, smoothing_disk);
  auto const t1 = birth_1 * (1.0 - sdisk) + death_1 * sdisk;
  auto const t2 = birth_2 * (1.0 - sdisk) + death_2 * sdisk;
  return Sigmoid(ring - t1, smoothing_ring) * Sigmoid(t2 - ring, smoothing_ring);
}

const std::vector<density> &Smooth::Field() const { return field; };
int Smooth::Index(int i, int j) const { return i * Sizex() + j; }

int Smooth::TorusDifference(int x1, int x2, int size) const {
  auto const remainder = std::abs(x1 - x2) % size;
  return std::min(remainder, std::abs(remainder - size));
}

double Smooth::Radius(int x1, int y1, int x2, int y2) const {
  int xdiff = TorusDifference(x1, x2, sizex);
  int ydiff = TorusDifference(y1, y2, sizey);
  return std::sqrt(xdiff * xdiff + ydiff * ydiff);
}

double Smooth::NormalisationDisk() const {
  double total = 0.0;
  for(int x = 0; x < sizex; x++)
    for(int y = 0; y < sizey; y++)
      total += Disk(Radius(0, 0, x, y));
  return total;
}

double Smooth::NormalisationRing() const {
  double total = 0.0;
  for(int x = 0; x < sizex; x++)
    for(int y = 0; y < sizey; y++)
      total += Ring(Radius(0, 0, x, y));
  return total;
}

void Smooth::Update() {
  for(int x = 0; x < sizex; x++)
    for(int y = 0; y < sizey; y++) {
      auto const integrals = Integrals(x, y);
      work_field[Index(x, y)] = Transition(integrals.first, integrals.second);
    }
  std::swap(field, work_field);
  frame++;
}

std::pair<density, density> Smooth::Integrals(int x, int y) const {
  density ring_total(0), disk_total(0);
  for(int x1 = 0; x1 < sizex; x1++) {
    int deltax = TorusDifference(x, x1, sizex);
    if(deltax > outer + smoothing / 2)
      continue;

    for(int y1 = 0; y1 < sizey; y1++) {
      int deltay = TorusDifference(y, y1, sizey);
      if(deltay > outer + smoothing / 2)
        continue;

      double radius = std::sqrt(deltax * deltax + deltay * deltay);
      double fieldv = field[Index(x1, y1)];
      ring_total += fieldv * Ring(radius);
      disk_total += fieldv * Disk(radius);
    }
  }
  return {disk_total / NormalisationDisk(), ring_total / NormalisationRing()};
}

void Smooth::SeedRandom() {
  for(int x = 0; x < sizex; x++)
    for(int y = 0; y < sizey; y++)
      field[Index(x, y)] += (static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
}

void Smooth::SeedConstant(density constant) { std::fill(field.begin(), field.end(), constant); }
void Smooth::AddDisk(int x0, int y0) {
  for(int x = 0; x < sizex; x++)
    for(int y = 0; y < sizey; y++)
      field[Index(x, y)] += Disk(Radius(x0, y0, x, y));
}

void Smooth::AddRing(int x0, int y0) {
  for(int x = 0; x < sizex; x++)
    for(int y = 0; y < sizey; y++)
      field[Index(x, y)] += Ring(Radius(x0, y0, x, y));
}

void Smooth::AddPixel(int x0, int y0, density value) { field[Index(x0, y0)] = value; }

void Smooth::Write(std::ostream &out) {
  for(int x = 0; x < sizex; x++) {
    for(int y = 0; y < sizey; y++)
      out << field[Index(x, y)] << " , ";
    out << std::endl;
  }
  out << std::endl;
}

int Smooth::Frame() const { return frame; }
