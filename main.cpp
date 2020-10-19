#include <iostream>
#include "Particle.hpp"

int main()
{
  Particle::AddParticleType("Proton", 10., 2);
  Particle::AddParticleType("Boson", 5., 1., 10.);
  Particle p1{"Proton", 1, 2, 3};
  Particle p2{"Boson", 4, 5, 6};
  std::cout << p2.GetPx() << '\n';
  std::cout << p2.Energy() << p2.InvMass(p1) << '\n';
  p1.Print();
  p2.PrintParticleTypes();
  p1.SetParticleType("Boson");
  p1.Print();
}