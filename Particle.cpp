#include "Particle.hpp"
#include <iostream>
#include "TMath.h"

int Particle::fNParticleType = 0;

ParticleType* Particle::fParticleTypes[fMaxNumParticleType] = {};

Particle::Particle(std::string name, double Px, double Py, double Pz)
    : fPx{Px}, fPy{Py}, fPz{Pz}
{
  fIParticle = FindParticle(name);
  if (fIParticle == -1) {
    std::cout << "Uknown ParticleType passed to Particle constructor" << '\n';
  }
}

void Particle::SetParticleType(std::string name)
{
  fIParticle = FindParticle(name);
  if (fIParticle == -1) {
    std::cout << "Uknown ParticleType passed to SetParticleType" << '\n';
  }
}

double Particle::Energy() const
{
  ParticleType* particalType = fParticleTypes[fIParticle];
  double mass = particalType->GetMass();
  return TMath::Sqrt(mass * mass + fPx * fPx + fPy * fPy + fPz * fPz);
}

double Particle::InvMass(Particle& particle2) const
{
  double E1 = Energy();
  double E2 = particle2.Energy();
  double Px2 = particle2.GetPx();
  double Py2 = particle2.GetPy();
  double Pz2 = particle2.GetPz();

  return TMath::Sqrt(
      // clang-format off
      (E1 + E2) * (E1 + E2) -
      ((fPx + Px2) * (fPx + Px2) + 
      (fPx + Px2) * (fPx + Px2) +
      (fPx + Px2) * (fPx + Px2))
      // clang-format on
  );
}

void Particle::Print() const
{
  ParticleType* particalType = fParticleTypes[fIParticle];
  std::cout << "IParticle:" << fIParticle << " Name:" << particalType->GetName()
            << " Pulse:(" << fPx << ',' << fPy << ',' << fPz << ")\n";
}

void Particle::AddParticleType(std::string name,
                               double mass,
                               int charge,
                               double width)
{
  if (fMaxNumParticleType != fNParticleType && FindParticle(name) == -1) {
    if (width == 0.) {
      fParticleTypes[fNParticleType] = new ParticleType(name, mass, charge);
    } else {
      fParticleTypes[fNParticleType] =
          new ResonanceType(name, mass, charge, width);
    }
    ++fNParticleType;
  }
}

void Particle::PrintParticleTypes()
{
  for (int i = 0; i != fNParticleType; ++i) {
    fParticleTypes[i]->Print();
  }
}

int Particle::FindParticle(std::string name)
{
  for (int i = 0; i != fNParticleType; ++i) {
    if (name == fParticleTypes[i]->GetName()) {
      return i;
    }
  }
  return -1;
}