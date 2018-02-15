// #include <string>
void FakeVetoProbabilityStability(int type = 0, double thMinus = 8, double thPlus = 8)
{
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);

  pair<string,vector<string>> runs8TeVPbp = {"C",{"285480","285505","285517","285530","285537","285538","285539","285549","285684","285718","285726","285739","285750","285759","285832"}};
  pair<string,vector<string>> runs8TeVpPb = {"C",{"285956", "285975", "285993", "285994", "285995", "286009", "286010", "286023", "286031", "286033", "286034", "286051", "286054", "286069", "286070", "286178", "286200", "286201", "286288", "286301",
                                                  "286302",
                                                  "286309", "286314", "286327", "286329", "286365", "286420", "286422", "286425", "286441", "286442", "286450", "286471", "286496"}};
  pair<string,vector<string>> runs5TeVPbp = {"B",{"285090", "285216", "285244", "285368", "285369", "285371", "285374", "285383"}};
  pair<string,vector<string>> runs5TeVpPb = {"D",{"286513", "286514", "286515", "286516", "286517", "286518", "286519", "286520"}};
  vector<pair<string,vector<string>>> runs = {runs8TeVPbp, runs8TeVpPb, runs5TeVPbp, runs5TeVpPb};

  struct runData{
    string name;
    pair<double,double> p_fake_minus;
    pair<double,double> p_fake_plus;
    pair<double,double> veto_efficiency_minus;
    pair<double,double> veto_efficiency_plus;
   };
  vector <runData> runsData;
  
  auto currentRun = runs.at(type);
  for(auto r: currentRun.second){
    string filename = string("store_hidata_PARun2016") + currentRun.first + string("_PAEmptyBX_AOD_PromptReco-v1_000_") +
      r.substr(0,3) + string("_") + r.substr(3,3) + string("/uaplot_bptxxor_without_tracks.root");
    auto f = TFile::Open(filename.c_str(),"r");
    auto hf_minus = (TH1F*) f->Get("CMScalo/hf_max_towerE_minus_h[0]");
    auto hf_plus = (TH1F*) f->Get("CMScalo/hf_max_towerE_plus_h[0]");
    runData d = {r};

    auto sumIMinus = hf_minus->Integral(1, hf_minus->GetNbinsX());
    auto sumIPlus = hf_plus->Integral(1, hf_plus->GetNbinsX());
    auto binNMinus = hf_minus->FindBin(thMinus);
    auto binNPlus = hf_plus->FindBin(thPlus);
    
    double eff = 0, errEff = 0, fvp = 0, errFvp = 0;

    eff = hf_minus->IntegralAndError(1, binNMinus, errEff);
    fvp = hf_minus->IntegralAndError(binNMinus, hf_minus->GetNbinsX(), errFvp);
    d.p_fake_minus = make_pair(fvp/sumIMinus, errFvp/sumIMinus);
    d.veto_efficiency_minus = make_pair(eff/sumIMinus, errEff/sumIMinus);

    eff = hf_plus->IntegralAndError(1, binNPlus, errEff);
    fvp = hf_plus->IntegralAndError(binNPlus, hf_plus->GetNbinsX(), errFvp);
    d.p_fake_plus = make_pair(fvp/sumIPlus, errFvp/sumIPlus);
    d.veto_efficiency_plus = make_pair(eff/sumIPlus, errEff/sumIPlus);

    f->Close();
    runsData.push_back(d);    
    cout << "For run " << r << " efficiencis: " << d.veto_efficiency_minus.first << "\t" << d.veto_efficiency_plus.first << endl;
  }

  auto outF = TFile::Open("FakeVetoProbabilityStability.root","recreate");

  auto hf_fvp_minus = new TH1F("FakeVetoProbabilityHF-","P_{Fake\\_HF-}", runsData.size(), 0, runsData.size());
  auto hf_fvp_plus = new TH1F("FakeVetoProbabilityHF+","P_{Fake\\_HF+}", runsData.size(), 0, runsData.size());
  auto hf_eff_minus = new TH1F("ThresholdEfficiencyHF-","ThresholdEfficiency\\_HF-", runsData.size(), 0, runsData.size());
  auto hf_eff_plus = new TH1F("ThresholdEfficiencyHF+","ThresholdEfficiency\\_HF+", runsData.size(), 0, runsData.size());

  hf_fvp_minus->SetLineColor(46);   hf_fvp_plus->SetLineColor(36);
  hf_fvp_minus->SetMarkerColor(46); hf_fvp_plus->SetMarkerColor(36);
  hf_eff_minus->SetLineColor(46);   hf_eff_plus->SetLineColor(36);
  hf_eff_minus->SetMarkerColor(46); hf_eff_plus->SetMarkerColor(36);

  for(auto i = 0; i < runsData.size(); i++){
    auto run = runsData.at(i);
    // cout << "Adding data for run " << run.name << endl;
    for(auto h: {hf_fvp_minus,hf_fvp_plus,hf_eff_minus,hf_eff_plus})
      h->GetXaxis()->SetBinLabel(i+1,run.name.c_str());
    hf_fvp_minus->SetBinContent(i+1,run.p_fake_minus.first);
    hf_fvp_minus->SetBinError(i+1,run.p_fake_minus.second);
    hf_fvp_plus->SetBinContent(i+1,run.p_fake_plus.first);
    hf_fvp_plus->SetBinError(i+1,run.p_fake_plus.second);
    hf_eff_minus->SetBinContent(i+1,run.veto_efficiency_minus.first);
    hf_eff_minus->SetBinError(i+1,run.veto_efficiency_minus.second);
    hf_eff_plus->SetBinContent(i+1,run.veto_efficiency_plus.first);
    hf_eff_plus->SetBinError(i+1,run.veto_efficiency_plus.second);
  }
  // for(auto h: {hf_fvp_minus,hf_fvp_plus,hf_eff_minus,hf_eff_plus})
  //   h->GetXaxis()->LabelsOption("v");

  auto c = new TCanvas("c1", "c1", 800, 600);
  gPad->SetMargin(0.07, 0.05, 0.10, 0.02);
  {
    char thM[10], thP[10]; sprintf(thM,"%.1f",thMinus); sprintf(thP,"%.1f",thPlus);
    string name = string("P_{Fake} for thresholds equal ") + thM + string(" and ") + thP + string(" for HF- and HF+");
    auto hs = new THStack (name.c_str(),name.c_str());
    auto leg = new TLegend(0.84, 0.89, 0.99, 0.99); //0.14, 0.14, 0.54, 0.34);
    hs->Add(hf_fvp_minus);
    leg->AddEntry(hf_fvp_minus);
    hs->Add(hf_fvp_plus);
    leg->AddEntry(hf_fvp_plus);

    hs->Draw("nostack,e1p");
    leg->Draw();
    c->SaveAs("p_fake.eps");
  };
  c->Clear();
  {
    char thM[10], thP[10]; sprintf(thM,"%.1f",thMinus); sprintf(thP,"%.1f",thPlus);
    string name = string("Threshold efficiency for thresholds equal ") + thM + string(" and ") + thP + string(" for HF- and HF+");
    auto hs = new THStack (name.c_str(),name.c_str());
    auto leg = new TLegend(0.84, 0.89, 0.99, 0.99); //0.14, 0.14, 0.54, 0.34);
    hs->Add(hf_eff_minus);
    leg->AddEntry(hf_eff_minus);
    hs->Add(hf_eff_plus);
    leg->AddEntry(hf_eff_plus);

    hs->Draw("nostack,e1p");
    leg->Draw();
    c->SaveAs("th_eff.eps");
  };
  c->Clear();
  outF->Write();
  outF->Close();

}
