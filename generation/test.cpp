{
  auto h = new TH1F("test", "test", 24, -4, 4);
  for (int i = 0; i != 1e6; ++i) {
    double y1 = gRandom->Gaus(0, 1);

    h->Fill(y1);
  }
  h->Fit("gaus");
  cout << h->GetRMS() * 1e6 / (1e6 - 1);
}
