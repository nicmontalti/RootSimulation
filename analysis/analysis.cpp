#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TList.h"
#include "TROOT.h"
#include "TStyle.h"

#include "myStyle.hpp"

void analysis(char* filePath)
{
  setMyStyle();
  TGaxis::SetMaxDigits(3);
  gROOT->SetStyle("myStyle");

  auto fileHistograms = new TFile(filePath);
  auto cGeneration = new TCanvas("cGeneration", "Particle Generation");
  auto cInvMass = new TCanvas("cInvMass", "Invariant mass");
  cInvMass->Divide(2, 2);
  cGeneration->Divide(2, 2);

  {
    cGeneration->cd(1);

    auto hParticleTypes = (TH1F*)fileHistograms->Get("hParticleTypes");
    int entries = hParticleTypes->GetEntries();
    std::cout << "PARTICLE TYPES";
    for (int i = 1; i != 8; ++i) {
      std::cout << "IParticle: " << i << '\t'
                << "Ratio: " << hParticleTypes->GetBinContent(i) / entries
                << " +/- " << hParticleTypes->GetBinError(i) / entries << '\n';
    }
    std::cout << '\n';

    hParticleTypes->UseCurrentStyle();
    hParticleTypes->SetTitle("Generated particles");
    hParticleTypes->SetXTitle("Particle");
    hParticleTypes->SetYTitle("Occurrences");
    hParticleTypes->GetXaxis()->SetBinLabel(1, "Pion+");
    hParticleTypes->GetXaxis()->SetBinLabel(2, "Pion-");
    hParticleTypes->GetXaxis()->SetBinLabel(3, "Kaon+");
    hParticleTypes->GetXaxis()->SetBinLabel(4, "Kaon-");
    hParticleTypes->GetXaxis()->SetBinLabel(5, "Proton+");
    hParticleTypes->GetXaxis()->SetBinLabel(6, "Proton-");
    hParticleTypes->GetXaxis()->SetBinLabel(7, "Kaon*");
    hParticleTypes->SetFillColor(kBlue - 8);
    hParticleTypes->SetLabelSize(0.06, "x");
    hParticleTypes->SetStats(0);

    hParticleTypes->Draw();

    auto hExpectedParticleTypes = new TH1F(*hParticleTypes);
    hExpectedParticleTypes->SetBinContent(1, 0.4 * entries);
    hExpectedParticleTypes->SetBinContent(2, 0.4 * entries);
    hExpectedParticleTypes->SetBinContent(3, 0.05 * entries);
    hExpectedParticleTypes->SetBinContent(4, 0.05 * entries);
    hExpectedParticleTypes->SetBinContent(5, 0.045 * entries);
    hExpectedParticleTypes->SetBinContent(6, 0.045 * entries);
    hExpectedParticleTypes->SetBinContent(7, 0.01 * entries);

    hExpectedParticleTypes->SetMarkerStyle(kCircle);
    hExpectedParticleTypes->SetMarkerSize(0.7);
    hExpectedParticleTypes->Draw("P,same");

    auto legend = new TLegend(0.7, 0.8, 0.95, 0.95);
    legend->AddEntry(hParticleTypes, "Generated particles", "f");
    legend->AddEntry(hExpectedParticleTypes, "Expected particles", "p");
    legend->Draw();
  }

  {
    cGeneration->cd(2);

    auto hPulse = (TH1F*)fileHistograms->Get("hPulse");
    hPulse->Fit("expo", "Q");
    auto fitFunc = hPulse->GetFunction("expo");

    std::cout << "PULSE FIT\n"
              << "Tau: " << -fitFunc->GetParameter(1) << " +/- "
              << fitFunc->GetParError(1) << '\n'
              << "Chi^2/dof: "
              << fitFunc->GetChisquare() /
                     (hPulse->GetNbinsX() - fitFunc->GetNpar())
              << "\n\n";

    hPulse->UseCurrentStyle();
    hPulse->GetXaxis()->SetTitle("Pulse (GeV)");
    hPulse->GetYaxis()->SetTitle("Occurences");
    hPulse->SetFillColor(42);
    hPulse->SetLineColor(kBlack);
    fitFunc->SetLineColor(kBlack);
    fitFunc->SetLineWidth(2);

    hPulse->Draw();
  }

  {
    cGeneration->cd(3);

    auto hAzimutalAngles = (TH1F*)fileHistograms->Get("hAzimutalAngles");
    auto hPolarAngles = (TH1F*)fileHistograms->Get("hPolarAngles");

    auto listAngles = new TList();
    listAngles->Add(hAzimutalAngles);
    listAngles->Add(hPolarAngles);

    for (int i = 0; i != 2; ++i) {
      cGeneration->cd(i + 3);
      auto h = (TH1F*)listAngles->At(i);

      h->Fit("pol0", "Q");
      auto fitFunc = h->GetFunction("pol0");

      std::cout << "AZIMUTAL ANGLES FIT\n"
                << "Parameter: " << fitFunc->GetParameter(0) << " +/- "
                << fitFunc->GetParError(0) << '\n'
                << "Chi^2/dof: "
                << fitFunc->GetChisquare() /
                       (h->GetNbinsX() - fitFunc->GetNpar())
                << "\n\n";

      h->UseCurrentStyle();
      h->SetXTitle("Angle (rad)");
      h->SetYTitle("Occurences");
      h->GetYaxis()->SetTitleOffset(1.26);
      h->SetMarkerStyle(kFullCircle);
      h->SetMarkerSize(0.6);
      h->SetMarkerColor(kBlack);
      h->SetLineColor(kBlack);
      fitFunc->SetLineColor(kRed + 1);

      h->Draw("E1P");
    }
  }

  {
    auto hConcordantInvMass = (TH1F*)fileHistograms->Get("hConcordantInvMass");
    auto hDiscordantInvMass = (TH1F*)fileHistograms->Get("hDiscordantInvMass");
    auto hConcordantPionKaonInvMass =
        (TH1F*)fileHistograms->Get("hConcordantPionKaonInvMass");
    auto hDiscordantPionKaonInvMass =
        (TH1F*)fileHistograms->Get("hDiscordantPionKaonInvMass");
    auto hResonanceCoupleInvMass =
        (TH1F*)fileHistograms->Get("hResonanceCoupleInvMass");

    auto hDifferencePionKaonInvMass = new TH1F(*hDiscordantPionKaonInvMass);
    hDifferencePionKaonInvMass->Add(
        hDiscordantPionKaonInvMass, hConcordantPionKaonInvMass, 1, -1);
    hDifferencePionKaonInvMass->GetXaxis()->SetRangeUser(0.89166 - 0.2,
                                                         0.89166 + 0.2);

    auto hDifferenceInvMass = new TH1F(*hDiscordantInvMass);
    hDifferenceInvMass->Add(hDiscordantInvMass, hConcordantInvMass, 1, -1);
    hDifferenceInvMass->GetXaxis()->SetRangeUser(0.89166 - 0.2, 0.89166 + 0.2);

    auto listInvMass = new TList();
    listInvMass->Add(hDifferencePionKaonInvMass);
    listInvMass->Add(hDifferenceInvMass);
    listInvMass->Add(hResonanceCoupleInvMass);

    int colors[3] = {45, 31, 38};

    for (int i = 0; i != 3; ++i) {
      cInvMass->cd(i + 1);
      auto h = (TH1F*)listInvMass->At(i);
      h->UseCurrentStyle();
      h->Fit("gaus", "Q");
      auto func = h->GetFunction("gaus");

      h->SetXTitle("Invariant mass (GeV/c^2)");
      h->SetYTitle("Occurences");
      h->SetFillColorAlpha(colors[i], 0.5);
      h->SetLineColor(kBlack);
      func->SetLineColor(kBlack);
      func->SetLineWidth(2);

      h->Draw();
    }
  }

  cGeneration->SaveAs("Generation.pdf");
  cGeneration->SaveAs("Generation.C");
  cGeneration->SaveAs("Generation.root");

  cInvMass->SaveAs("InvMass.pdf");
  cInvMass->SaveAs("InvMass.C");
  cInvMass->SaveAs("InvMass.root");
}

int main(int argc, char* argv[]) { analysis(argv[1]); }