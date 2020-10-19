#include "ResonanceType.hpp"
#include <iostream>

ResonanceType::ResonanceType(std::string name,
                             double mass,
                             int charge,
                             double width)
    : ParticleType{name, mass, charge}, fWidth{width}
{
}

void ResonanceType::Print() const
{
  std::cout << "Name:" << GetName() << " Mass:" << GetMass()
            << " Charge:" << GetCharge() << " Width:" << fWidth << '\n';
}