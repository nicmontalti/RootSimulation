#ifndef PARTICLE_TYPE_HPP
#define PARTICLE_TYPE_HPP

#include <string>

class ParticleType
{
 public:
  ParticleType(std::string name = std::string{},
               double mass = 0.,
               int charge = 0);
  std::string GetName() const { return fName; }
  double GetMass() const { return fMass; }
  double GetCharge() const { return fCharge; }
  virtual void Print() const;
  virtual double GetWidth() const;

 private:
  const std::string fName;
  const double fMass;
  const int fCharge;
};

#endif  // PARTICLE_TYPE_HPP