#include "Particle.hpp"
#include <cmath>    // for M_PI
#include <cstdlib>  //for RAND_MAX
#include <iostream>

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

double Particle::GetEnergy() const
{
  ParticleType* particalType = fParticleTypes[fIParticle];
  double mass = particalType->GetMass();
  return std::sqrt(mass * mass + fPx * fPx + fPy * fPy + fPz * fPz);
}

double Particle::InvMass(Particle& particle2) const
{
  double E1 = GetEnergy();
  double E2 = particle2.GetEnergy();
  double Px2 = particle2.GetPx();
  double Py2 = particle2.GetPy();
  double Pz2 = particle2.GetPz();

  return std::sqrt(
      // clang-format off
      (E1 + E2) * (E1 + E2) -
      ((fPx + Px2) * (fPx + Px2) + 
      (fPy + Py2) * (fPy + Py2) +
      (fPz + Pz2) * (fPz + Pz2))
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

int Particle::Decay2body(Particle& dau1, Particle& dau2) const
{
  if (GetMass() == 0.0) {
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }

  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if (fIParticle > -1) {  // add width effect

    // gaussian random numbers

    float x1, x2, w, y1, y2;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;

    massMot += fParticleTypes[fIParticle]->GetWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    printf(
        "Decayment cannot be preformed because mass is too low in this "
        "channel\n");
    return 2;
  }

  double pout =
      sqrt(
          (massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) *
          (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) /
      massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;
  dau1.SetP(pout * sin(theta) * cos(phi),
            pout * sin(theta) * sin(phi),
            pout * cos(theta));
  dau2.SetP(-pout * sin(theta) * cos(phi),
            -pout * sin(theta) * sin(phi),
            -pout * cos(theta));

  double energy = sqrt(fPx * fPx + fPy * fPy + fPz * fPz + massMot * massMot);

  double bx = fPx / energy;
  double by = fPy / energy;
  double bz = fPz / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);

  return 0;
}
void Particle::Boost(double bx, double by, double bz)
{
  double energy = GetEnergy();

  // Boost this Lorentz vector
  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx * fPx + by * fPy + bz * fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  fPx += gamma2 * bp * bx + gamma * bx * energy;
  fPy += gamma2 * bp * by + gamma * by * energy;
  fPz += gamma2 * bp * bz + gamma * bz * energy;
}