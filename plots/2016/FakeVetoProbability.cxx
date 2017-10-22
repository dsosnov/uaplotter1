// #include <string>
void FakeVetoProbability(std::string ifilename = "uaplot.root")
{
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open(ifilename.c_str(),"r");
  auto outF = TFile::Open("fake_veto_probability.root","recreate");
  // TH1F *hcal_h, *hf_p_h, *hf_m_h;
  TString hN;
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  vector<pair<int,string>> typesOfBptx = { make_pair(0,"bptxxor_"),make_pair(1,"bptxplus_"),make_pair(2,"bptxminus_")};
  const double PERCENT_OF_ENTRIES = 99.7 / 100.0;
  //! For caloTower threshold for plus and minus side 
  std::vector<std::tuple<TString, std::string, std::string>> histogramData = {
    std::make_tuple("CMScalo/hf_max_towerE_minus_h[", "Calo Towers threshold for HF minus side", "hf_max_towerE_minus_h.eps"),
    std::make_tuple("CMScalo/hf_max_towerE_plus_h[", "Calo Towers threshold for HF plus side", "hf_max_towerE_plus_h.eps"),
    // std::make_tuple("CMScalo/hf_max_towerE_minus_ill_h[", "Calo Towers threshold for ill tower on HF minus side", "hf_max_towerE_minus_ill_h.eps"),
    // std::make_tuple("CMScalo/hf_max_towerE_plus_ill_h[", "Calo Towers threshold for ill tower on HF plus side", "hf_max_towerE_plus_ill_h.eps"),
  };
  // auto plothfm = TH1F("FakeVetoProbabilityHF-","FakeVetoProbability\\_HF-", 190, 1, 20);
  // auto plothfp = TH1F("FakeVetoProbabilityHF+","FakeVetoProbability\\_HF+", 190, 1, 20);
  vector<TH1F*> plots = {
    new TH1F("FakeVetoProbabilityHF-","P_{Fake\\_HF-}", 190, 1, 20),
    new TH1F("FakeVetoProbabilityHF+","P_{Fake\\_HF+}", 190, 1, 20),
    new TH1F("ThresholdEfficiencyHF-","ThresholdEfficiency\\_HF-", 190, 1, 20),
    new TH1F("ThresholdEfficiencyHF+","ThresholdEfficiency\\_HF+", 190, 1, 20)
  };
    // plothfm, plothfp};
  auto currentPlot = 0;
  // for(auto p: typesOfBptx)
  for(auto t: histogramData){
    TString histoName=std::get<0>(t);
    hN = histoName; hN+=0; hN+="]";
    std::string name = std::get<1>(t);
    // TCanvas *c4 = new TCanvas(name.c_str(), "", 800, 200); // Canvas for printing of projections by bin
    // c4->SetLogy(); c4->SetLogx();
    TH1F *histo = (TH1F*) f->Get(hN.Data());
    // histo->SetTitle((name).c_str());
    histo->Sumw2();
    int sumI = histo->Integral(1, histo->GetNbinsX());
    // sumI=1;
    // histo->Scale(1.0/sumI,"nosw2");
    for(double th = 1; th < 20; th+=0.5){
      int binN = histo->FindBin(th);
      double errEff = 0, errFvp = 0;
      double eff = histo->IntegralAndError(1, binN, errEff);
      eff/=sumI; errEff/=sumI;
      double fvp = histo->IntegralAndError(binN, histo->GetNbinsX(), errFvp);
      fvp/=sumI; errFvp/=sumI;
      plots.at(currentPlot)->Fill(th,fvp);
      plots.at(currentPlot)->SetBinError(plots.at(currentPlot)->FindBin(th),errFvp);

      plots.at(currentPlot+2)->Fill(th,eff);
      plots.at(currentPlot+2)->SetBinError(plots.at(currentPlot+2)->FindBin(th),errEff);
      // purityMinus->SetBinError(purityMinus->FindBin(param.thMinus), param.purity[2].second);
      cout << "For histogram \"" << name << "\" (" << hN << ")" << " fake veto probability at " << th << "(bin " << binN << ")" << ": " << fvp << "\t" << eff << endl;
      // cout << "\tSum is: " << (fvp + eff)/sumI << " " << errEff/sumI << " " << errFvp/sumI << endl;
    }
    // for(uint i = 1010; i < 1300; i+=1)
    //   cout << "\t\t: " << i << "\t" << histo->GetBinContent(i) << "\t" << histo->GetBinError(i) << endl;
    currentPlot++;
  }
  f->Close();
  outF->Write();
  outF->Close();
}
