#include "Particle.hpp"
#include "TBenchmark.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include <array>

int constexpr maxNumParticle = 100;
int constexpr maxNumDecay = 20;

void AddParticleTypes() {
  Particle::AddParticleType("pion+", 0.13957, 1);
  Particle::AddParticleType("pion-", 0.13957, -1);
  Particle::AddParticleType("kaon+", 0.49367, 1);
  Particle::AddParticleType("kaon-", 0.49367, -1);
  Particle::AddParticleType("proton+", 0.93827, 1);
  Particle::AddParticleType("proton-", 0.93827, -1);
  Particle::AddParticleType("kaon*", 0.89166, 0, 0.050);
}

int GenerateParticles(
    std::array<Particle, maxNumParticle + maxNumDecay> &particles) {
  int nDecay = 0;

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
        particles[maxNumParticle + nDecay].SetParticleType("pion+");
        particles[maxNumParticle + nDecay + 1].SetParticleType("kaon-");
      } else {
        particles[maxNumParticle + nDecay].SetParticleType("pion-");
        particles[maxNumParticle + nDecay + 1].SetParticleType("kaon+");
      }

      particles[i].Decay2body(particles[maxNumParticle + nDecay],
                              particles[maxNumParticle + nDecay + 1]);

      nDecay += 2;
    }
  }
  return nDecay;
}

TList *FillHistos(int nEvents = 1e5) {
  auto hParticleTypes = new TH1F("hParticleTypes", "Particle types", 7, 0, 7);
  auto hAzimutalAngles = new TH1F(
      "hAzimutalAngles", "Azimutal angles distribution", 10, 0, TMath::Pi());
  auto hPolarAngles = new TH1F("hPolarAngles", "Polar angles distribution", 10,
                               0, 2 * TMath::Pi());
  auto hAngles = new TH2F("hAngles", "Angles distribution", 10, 10, 0,
                          TMath::Pi(), 0, 2 * TMath::Pi());
  auto hPulse = new TH1F("hPulse", "Pulse", 20, 0, 5);
  auto hTransversePulse =
      new TH1F("hTransversePulse", "Transverse pulse", 20, 0, 5);
  auto hEnergy = new TH1F("hEnergy", "Energy", 200, 0, 5);
  auto hInvMass = new TH1F("hInvMass", "Invariant mass", 200, 0, 5);
  auto hConcordantInvMass = new TH1F(
      "hConcordantInvMass",
      "Invariant mass of particles with concordant charge sign", 200, 0, 5);
  auto hDiscordantInvMass = new TH1F(
      "hDiscordantInvMass",
      "Invariant mass of particles with discordant charge sign", 200, 0, 5);
  auto hConcordantPionKaonInvMass =
      new TH1F("hConcordantPionKaonInvMass",
               "Invariant mass of kaons and pions with concordant charge sign",
               200, 0, 5);
  auto hDiscordantPionKaonInvMass =
      new TH1F("hDiscordantPionKaonInvMass",
               "Invariant mass of kaons and pions with discordant charge sign",
               200, 0, 5);
  auto hResonanceCoupleInvMass = new TH1F(
      "hResonanceCoupleInvMass", "Invariant mass of decayed particles couples",
      24, 0.89166 - 0.2, 0.89166 + 0.2);

  for (int i = 0; i != nEvents; ++i) {
    std::array<Particle, maxNumParticle + maxNumDecay> particles;
    int nResonanceDecay = GenerateParticles(particles);

    for (int i = 0; i != nResonanceDecay + maxNumParticle; ++i) {
      auto &particle = particles[i];

      if (i < maxNumParticle) {
        double Px = particle.GetPx();
        double Py = particle.GetPy();
        double Pz = particle.GetPz();

        auto phi = TMath::ATan2(Py, Px) + TMath::Pi();
        auto theta = TMath::ATan2(TMath::Sqrt(Px * Px + Py * Py), Pz);
        auto P = TMath::Sqrt(Px * Px + Py * Py + Pz * Pz);
        auto transP = TMath::Sqrt(Px * Px + Py * Py);

        hParticleTypes->Fill(particle.GetIParticle());
        hPolarAngles->Fill(phi);
        hAzimutalAngles->Fill(theta);
        hAngles->Fill(theta, phi);
        hPulse->Fill(P);
        hTransversePulse->Fill(transP);
        hEnergy->Fill(particle.GetEnergy());
      }

      for (int j = i + 1; j != maxNumParticle + nResonanceDecay; ++j) {
        auto &particle2 = particles[j];
        auto name = particle.GetName();
        auto name2 = particle2.GetName();
        double invMass = particle.InvMass(particle2);

        hInvMass->Fill(invMass);

        if (particle.GetCharge() * particle2.GetCharge() > 0) {
          hConcordantInvMass->Fill(invMass);
          if ((name == "pion+" && name2 == "kaon+") ||
              (name == "pion-" && name2 == "kaon-")) {
            hConcordantPionKaonInvMass->Fill(invMass);
          }
        } else if (particle.GetCharge() * particle2.GetCharge() < 0) {
          hDiscordantInvMass->Fill(invMass);
          if ((name == "pion+" && name2 == "kaon-") ||
              (name == "pion-" && name2 == "kaon+")) {
            hDiscordantPionKaonInvMass->Fill(invMass);
          }
        }
      }
    }
    for (int j = maxNumParticle; j != maxNumParticle + nResonanceDecay;
         j += 2) {
      double invMass = particles[j].InvMass(particles[j + 1]);
      hResonanceCoupleInvMass->Fill(invMass);
    }
  }

  auto listHistos = new TList{};
  listHistos->SetOwner();

  listHistos->Add(hParticleTypes);
  listHistos->Add(hAzimutalAngles);
  listHistos->Add(hPolarAngles);
  listHistos->Add(hAngles);
  listHistos->Add(hPulse);
  listHistos->Add(hTransversePulse);
  listHistos->Add(hEnergy);
  listHistos->Add(hInvMass);
  listHistos->Add(hConcordantInvMass);
  listHistos->Add(hDiscordantInvMass);
  listHistos->Add(hConcordantPionKaonInvMass);
  listHistos->Add(hDiscordantPionKaonInvMass);
  listHistos->Add(hResonanceCoupleInvMass);

  return listHistos;
}

int main() {
  gRandom->SetSeed();
  AddParticleTypes();

  auto listHistos = FillHistos(1e5);

  auto file = new TFile("Histograms.root", "RECREATE");
  for (auto const histo : *listHistos) {
    histo->Write();
  }
  file->Close();

  delete listHistos;
}