#include <array>
#include "Particle.hpp"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"

int constexpr maxNumParticle = 100;
int constexpr maxResonanceDecay = 20;

void AddParticleTypes()
{
  Particle::AddParticleType("pion+", 0.13957, 1);
  Particle::AddParticleType("pion-", 0.13957, -1);
  Particle::AddParticleType("kaon+", 0.49367, 1);
  Particle::AddParticleType("kaon-", 0.49367, -1);
  Particle::AddParticleType("proton+", 0.93827, 1);
  Particle::AddParticleType("proton-", 0.93827, -1);
  Particle::AddParticleType("kaon*", 0.89166, 0, 0.050);
}

int GenerateParticles(
    std::array<Particle, maxNumParticle + maxResonanceDecay>& particles)
{
  int nResonanceDecay = 0;

  for (int i = 0; i != maxNumParticle; ++i) {
    double phi = gRandom->Uniform(2 * TMath::Pi());
    double theta = gRandom->Uniform(TMath::Pi());
    double P = gRandom->Exp(1);

    double Px = P * TMath::Sin(theta) * TMath::Cos(phi);
    double Py = P * TMath::Sin(theta) * TMath::Sin(phi);
    double Pz = P * TMath::Cos(theta);

    particles[i].SetP(Px, Py, Pz);

    double x = gRandom->Rndm();
    if (x < 0.4) {
      particles[i].SetParticleType("pion+");
    } else if (x < 0.8) {
      particles[i].SetParticleType("pion-");
    } else if (x < 0.85) {
      particles[i].SetParticleType("kaon+");
    } else if (x < 0.9) {
      particles[i].SetParticleType("kaon-");
    } else if (x < 0.945) {
      particles[i].SetParticleType("proton+");
    } else if (x < 0.99) {
      particles[i].SetParticleType("proton-");
    } else {
      particles[i].SetParticleType("kaon*");

      double y = gRandom->Rndm();
      if (y < 0.5) {
        particles[maxNumParticle + nResonanceDecay].SetParticleType("pion+");
        particles[maxNumParticle + nResonanceDecay + 1].SetParticleType(
            "kaon-");
      } else {
        particles[maxNumParticle + nResonanceDecay].SetParticleType("pion-");
        particles[maxNumParticle + nResonanceDecay + 1].SetParticleType(
            "kaon+");
      }

      particles[i].Decay2body(particles[maxNumParticle + nResonanceDecay],
                              particles[maxNumParticle + nResonanceDecay + 1]);

      nResonanceDecay += 2;
    }
  }
  return nResonanceDecay;
}

int main()
{
  gRandom->SetSeed();
  int nEvents = 1e5;

  AddParticleTypes();

  auto hParticleTypes = new TH1F("hParticleTypes", "Particle types", 7, 0, 7);
  auto hAzimutalAngles = new TH1F(
      "hAzimutalAngles", "Azimutal angles distribution", 100, 0, TMath::Pi());
  auto hPolarAngles = new TH1F(
      "hPolarAngles", "Polar angles distribution", 100, 0, 2 * TMath::Pi());
  auto hAngles = new TH2F("hAngles",
                          "Angles distribution",
                          100,
                          100,
                          0,
                          TMath::Pi(),
                          0,
                          2 * TMath::Pi());
  auto hPulse = new TH1F("hPulse", "Pulse", 100, 0, 5);
  auto hTransversePulse =
      new TH1F("hTransversePulse", "Transverse pulse", 100, 0, 5);
  auto hEnergy = new TH1F("hEnergy", "Energy", 100, 0, 5);
  auto hInvMass = new TH1F("hInvMass", "Invariant mass", 100, 0, 5);
  auto hConcordantInvMass =
      new TH1F("hConcordantInvMass",
               "Invariant mass of particles with concordant charge sign",
               100,
               0,
               5);
  auto hDiscordantInvMass =
      new TH1F("hDiscordantInvMass",
               "Invariant mass of particles with discordant charge sign",
               100,
               0,
               5);
  auto hConcordantPionKaonInvMass =
      new TH1F("hConcordantPionKaonInvMass",
               "Invariant mass of kaons and pions with concordant charge sign",
               100,
               0,
               5);
  auto hDiscordantPionKaonInvMass =
      new TH1F("hDiscordantPionKaonInvMass",
               "Invariant mass of kaons and pions with discordant charge sign",
               100,
               0,
               5);
  auto hResonanceCoupleInvMass =
      new TH1F("hResonanceCoupleInvMass",
               "Invariant mass of decayed particles couples",
               100,
               0.5,
               1.5);

  for (int i = 0; i != nEvents; ++i) {
    std::array<Particle, maxNumParticle + maxResonanceDecay> particles;
    int nResonanceDecay = GenerateParticles(particles);

    for (int i = 0; i != nResonanceDecay + maxNumParticle; ++i) {
      auto& particle = particles[i];

      auto phi = TMath::ATan(particle.GetPy() / particle.GetPx());
      auto theta =
          TMath::ATan(TMath::Sqrt(particle.GetPx() * particle.GetPx() +
                                  particle.GetPy() * particle.GetPy()) /
                      particle.GetPz());
      auto P = TMath::Sqrt(particle.GetPx() * particle.GetPx() +
                           particle.GetPy() * particle.GetPy() +
                           particle.GetPz() * particle.GetPz());
      auto transP = TMath::Sqrt(particle.GetPx() * particle.GetPx() +
                                particle.GetPy() * particle.GetPy());

      hParticleTypes->Fill(particle.GetIParticle());
      hPolarAngles->Fill(phi);
      hAzimutalAngles->Fill(theta);
      hAngles->Fill(theta, phi);
      hPulse->Fill(P);
      hTransversePulse->Fill(transP);
      hEnergy->Fill(particle.GetEnergy());

      for (int j = i + 1; j != maxNumParticle + nResonanceDecay; ++j) {
        auto& particle2 = particles[j];
        double invMass = particle.InvMass(particle2);
        hInvMass->Fill(invMass);
        if (particle.GetCharge() * particle2.GetCharge() > 0) {
          hConcordantInvMass->Fill(invMass);
          if ((particle.GetName() == "pion+" &&
               particle2.GetName() == "kaon+") ||
              (particle.GetName() == "pion-" &&
               particle2.GetName() == "kaon-")) {
            hConcordantPionKaonInvMass->Fill(invMass);
          }
        } else {
          hDiscordantInvMass->Fill(invMass);
          if ((particle.GetName() == "pion+" &&
               particle2.GetName() == "kaon-") ||
              (particle.GetName() == "pion-" &&
               particle2.GetName() == "kaon+")) {
            hDiscordantPionKaonInvMass->Fill(invMass);
          }
        }
      }
      for (int j = maxNumParticle; j != maxNumParticle + nResonanceDecay;
           j += 2) {
        double invMass = particles[j].InvMass(particles[j + 1]);
        hResonanceCoupleInvMass->Fill(invMass);
      }
    }
  }
  auto file = new TFile("Histograms.root", "RECREATE");
  hParticleTypes->Write();
  hAzimutalAngles->Write();
  hPolarAngles->Write();
  hAngles->Write();
  hPulse->Write();
  hTransversePulse->Write();
  hEnergy->Write();
  hInvMass->Write();
  hConcordantInvMass->Write();
  hDiscordantInvMass->Write();
  hConcordantPionKaonInvMass->Write();
  hDiscordantPionKaonInvMass->Write();
  hResonanceCoupleInvMass->Write();
  file->Close();
}