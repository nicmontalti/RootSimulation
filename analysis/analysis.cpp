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

void analysis(char *filePath) {
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

    auto hParticleTypes = (TH1D *)fileHistograms->Get("hParticleTypes");

    std::cout << "PARTICLE TYPES\n";
    for (int i = 1; i != 8; ++i) {
      std::cout << "IParticle: " << i << '\t'
                << "Occurrences: " << hParticleTypes->GetBinContent(i)
                << " +/- " << hParticleTypes->GetBinError(i) << '\n';
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

    auto hExpectedParticleTypes = new TH1D(*hParticleTypes);
    int entries = hParticleTypes->GetEntries();

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

    auto hPulse = (TH1D *)fileHistograms->Get("hPulse");
    hPulse->Fit("expo", "Q");
    auto fitFunc = hPulse->GetFunction("expo");
    double Chi = fitFunc->GetChisquare();
    int dof = fitFunc->GetNDF();

    std::cout << "PULSE FIT\n"
              << "Tau: " << -fitFunc->GetParameter(1) << " +/- "
              << fitFunc->GetParError(1) << '\n'
              << "Chi^2: " << Chi << '\n'
              << "dof: " << dof << '\n'
              << "Chi^2/dof: " << Chi / dof << "\n\n";

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

    auto hAzimutalAngles = (TH1F *)fileHistograms->Get("hAzimutalAngles");
    auto hPolarAngles = (TH1F *)fileHistograms->Get("hPolarAngles");

    auto listAngles = new TList();
    listAngles->Add(hAzimutalAngles);
    listAngles->Add(hPolarAngles);

    for (int i = 0; i != 2; ++i) {
      cGeneration->cd(i + 3);
      auto h = (TH1F *)listAngles->At(i);

      h->Fit("pol0", "Q");
      auto fitFunc = h->GetFunction("pol0");
      double Chi = fitFunc->GetChisquare();
      int dof = fitFunc->GetNDF();

      std::cout << h->GetTitle() << " FIT\n"
                << "Parameter: " << fitFunc->GetParameter(0) << " +/- "
                << fitFunc->GetParError(0) << '\n'
                << "Chi^2: " << Chi << '\n'
                << "dof: " << dof << '\n'
                << "Chi^2/dof: " << Chi / dof << "\n\n";

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
    auto hConcordantInvMass = (TH1D *)fileHistograms->Get("hConcordantInvMass");
    auto hDiscordantInvMass = (TH1D *)fileHistograms->Get("hDiscordantInvMass");
    auto hConcordantPionKaonInvMass =
        (TH1F *)fileHistograms->Get("hConcordantPionKaonInvMass");
    auto hDiscordantPionKaonInvMass =
        (TH1F *)fileHistograms->Get("hDiscordantPionKaonInvMass");
    auto hResonanceCoupleInvMass =
        (TH1F *)fileHistograms->Get("hResonanceCoupleInvMass");
    hResonanceCoupleInvMass->Rebin(2);

    auto hDifferencePionKaonInvMass =
        (TH1F *)hDiscordantPionKaonInvMass->Clone("hDifferencePionKaonInvMass");
    hDifferencePionKaonInvMass->Sumw2();
    hDifferencePionKaonInvMass->Add(hDiscordantPionKaonInvMass,
                                    hConcordantPionKaonInvMass, 1, -1);
    hDifferencePionKaonInvMass->GetXaxis()->SetRangeUser(0.89166 - 0.25,
                                                         0.89166 + 0.25);
    hDiscordantInvMass->Sumw2();
    auto hDifferenceInvMass = (TH1D *)hDiscordantInvMass->Clone("hDiff"
                                                                "erenc"
                                                                "eInvM"
                                                                "ass");
    hDifferenceInvMass->Add(hDiscordantInvMass, hConcordantInvMass, 1, -1);
    hDifferenceInvMass->Rebin(2);
    hDifferenceInvMass->GetXaxis()->SetRangeUser(0.89166 - 0.25,
                                                 0.89166 + 0.25);

    auto listInvMass = new TList();
    listInvMass->Add(hDifferenceInvMass);
    listInvMass->Add(hDifferencePionKaonInvMass);
    listInvMass->Add(hResonanceCoupleInvMass);

    int colors[3] = {45, 31, 38};

    for (int i = 0; i != 3; ++i) {
      cInvMass->cd(i + 1);
      auto h = (TH1D *)listInvMass->At(i);

      h->Fit("gaus", "Q");
      auto fitFunc = h->GetFunction("gaus");
      double Chi = fitFunc->GetChisquare();
      int dof = fitFunc->GetNDF();

      std::cout << "GAUS FIT " << i << '\n'
                << "Amplitude : " << fitFunc->GetParameter(0) << " +/- "
                << fitFunc->GetParError(0) << '\n'
                << "Mean : " << fitFunc->GetParameter(1) << " +/- "
                << fitFunc->GetParError(1) << '\n'
                << "Sigma : " << fitFunc->GetParameter(2) << " +/- "
                << fitFunc->GetParError(2) << '\n'
                << "Chi^2: " << Chi << '\n'
                << "dof: " << dof << '\n'
                << "Chi^2/dof: " << Chi / dof << "\n\n";

      h->UseCurrentStyle();
      h->SetXTitle("Invariant mass "
                   "(GeV/c^2)");
      h->SetYTitle("Occurences");
      h->SetFillColorAlpha(colors[i], 0.5);
      h->SetLineColor(kBlack);
      fitFunc->SetLineColor(kBlack);
      fitFunc->SetLineWidth(2);
      h->Draw("E1, X0");
      h->Draw("HIST, same");
      h->Draw("E1, X0, same");
    }
  }

  cGeneration->SaveAs("Generation.pdf");
  cGeneration->SaveAs("Generation.C");
  cGeneration->SaveAs("Generation.root");

  cInvMass->SaveAs("InvMass.pdf");
  cInvMass->SaveAs("InvMass.C");
  cInvMass->SaveAs("InvMass.root");
}

int main(int argc, char *argv[]) { analysis(argv[1]); }