#include "ParticleType.hpp"
#include <iostream>

ParticleType::ParticleType(std::string name, double mass, int charge)
    : fName{name}, fMass{mass}, fCharge{charge}
{
}

void ParticleType::Print() const
{
  std::cout << "Name:" << fName << " Mass:" << fMass << " Charge:" << fCharge
            << '\n';
}

double ParticleType::GetWidth() const { return 0.; }