#include "Particle.hpp"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include <array>
#include <iostream>

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

TList *FillHistos(int nEvents = 1e6) {
  auto hParticleTypes = new TH1D("hParticleTypes", "Particle types", 7, 0, 7);
  auto hAzimutalAngles = new TH1F(
      "hAzimutalAngles", "Azimutal angles distribution", 10, 0, TMath::Pi());
  auto hPolarAngles = new TH1F("hPolarAngles", "Polar angles distribution", 10,
                               0, 2 * TMath::Pi());
  auto hAngles = new TH2F("hAngles", "Angles distribution", 10, 10, 0,
                          TMath::Pi(), 0, 2 * TMath::Pi());
  auto hPulse = new TH1D("hPulse", "Pulse", 20, 0, 5);
  auto hTransversePulse =
      new TH1D("hTransversePulse", "Transverse pulse", 20, 0, 5);
  auto hEnergy = new TH1D("hEnergy", "Energy", 20, 0, 5);
  auto hInvMass = new TH1D("hInvMass", "Invariant mass", 20, 0, 5);
  auto hConcordantInvMass = new TH1D(
      "hConcordantInvMass",
      "Invariant mass of particles with concordant charge sign", 200, 0, 5);
  auto hDiscordantInvMass = new TH1D(
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

  for (int k = 0; k != nEvents; ++k) {
    std::array<Particle, maxNumParticle + maxNumDecay> particles;
    int nDecay = 0;

    for (int i = 0; i != maxNumParticle; ++i) {
      double phi = gRandom->Uniform(2 * TMath::Pi());
      double theta = gRandom->Uniform(TMath::Pi());
      double P = gRandom->Exp(1);

      double Px = P * TMath::Sin(theta) * TMath::Cos(phi);
      double Py = P * TMath::Sin(theta) * TMath::Sin(phi);
      double Pz = P * TMath::Cos(theta);
      double transP = TMath::Sqrt(Px * Px + Py * Py);

      particles[i].SetP(Px, Py, Pz);

      hPolarAngles->Fill(phi);
      hAzimutalAngles->Fill(theta);
      hAngles->Fill(theta, phi);
      hPulse->Fill(P);
      hTransversePulse->Fill(transP);

      double x = gRandom->Rndm();
      if (x < 0.4) {
        particles[i].SetParticleType("pion+");
        hParticleTypes->Fill(0);
      } else if (x < 0.8) {
        particles[i].SetParticleType("pion-");
        hParticleTypes->Fill(1);
      } else if (x < 0.85) {
        particles[i].SetParticleType("kaon+");
        hParticleTypes->Fill(2);
      } else if (x < 0.9) {
        particles[i].SetParticleType("kaon-");
        hParticleTypes->Fill(3);
      } else if (x < 0.945) {
        particles[i].SetParticleType("proton+");
        hParticleTypes->Fill(4);
      } else if (x < 0.99) {
        particles[i].SetParticleType("proton-");
        hParticleTypes->Fill(5);
      } else {
        particles[i].SetParticleType("kaon*");
        hParticleTypes->Fill(6);

        double y = gRandom->Rndm();
        auto &childParticle1 = particles[maxNumParticle + nDecay];
        auto &childParticle2 = particles[maxNumParticle + nDecay + 1];

        if (y < 0.5) {
          childParticle1.SetParticleType("pion+");
          childParticle2.SetParticleType("kaon-");
        } else {
          childParticle1.SetParticleType("pion-");
          childParticle2.SetParticleType("kaon+");
        }

        particles[i].Decay2body(childParticle1, childParticle2);

        double invMass = childParticle1.InvMass(childParticle2);
        hResonanceCoupleInvMass->Fill(invMass);

        nDecay += 2;
      }
      hEnergy->Fill(particles[i].GetEnergy());
    }

    for (int i = 0; i != maxNumParticle + nDecay; ++i) {
      auto &particle = particles[i];

      for (int j = i + 1; j != maxNumParticle + nDecay; ++j) {
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
    if ((k % 10000) == 0) {
      std::cout << k / 10000 << "%\n";
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

  auto listHistos = FillHistos(1e6);

  auto file = new TFile("Histograms.root", "RECREATE");
  for (auto const histo : *listHistos) {
    histo->Write();
  }
  file->Close();

  delete listHistos;
}