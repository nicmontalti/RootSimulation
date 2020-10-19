#ifndef RESONANCE_TYPE_HPP
#define RESONANCE_TYPE_HPP

#include <string>
#include "ParticleType.hpp"

class ResonanceType : public ParticleType
{
 public:
  ResonanceType(std::string name = std::string{},
                double mass = 0.,
                int charge = 0,
                double width = 0.);
  double GetWidth() const { return fWidth; }
  void Print() const override;

 private:
  const double fWidth;
};

#endif  // RESONANCE_TYPE_HPP