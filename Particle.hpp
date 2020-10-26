#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <string>
#include "ResonanceType.hpp"

class Particle
{
 public:
  Particle(std::string name = {},
           double Px = 0.,
           double Py = 0.,
           double Pz = 0.);

  int GetIParticle() const { return fIParticle; }
  void SetParticleType(int IParticle) { fIParticle = IParticle; }
  void SetParticleType(std::string name);

  double GetPx() const { return fPx; }
  double GetPy() const { return fPy; }
  double GetPz() const { return fPz; }

  void SetP(double x, double y, double z)
  {
    fPx = x;
    fPy = y;
    fPz = z;
  }

  std::string GetName() const { return fParticleTypes[fIParticle]->GetName(); }
  double GetCharge() const { return fParticleTypes[fIParticle]->GetCharge(); }
  double GetMass() const { return fParticleTypes[fIParticle]->GetMass(); }
  double GetEnergy() const;
  double InvMass(Particle& particle2) const;

  void Print() const;

  int Decay2body(Particle& dau1, Particle& dau2) const;

  static void AddParticleType(std::string name,
                              double mass,
                              int charge,
                              double width = 0.);
  static void PrintParticleTypes();

 private:
  int fIParticle;

  double fPx;
  double fPy;
  double fPz;

  static const int fMaxNumParticleType = 10;
  static int fNParticleType;
  static ParticleType* fParticleTypes[];
  static int FindParticle(std::string name);

  void Boost(double bx, double by, double bz);
};

#endif  // PARTICLE_HPP