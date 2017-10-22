// #include <string>
void plotEnergyDistribution(std::string additionalText = "", std::string ifilename = "uaplot.root")
{
// string additionalText=""; string ifilename = "uaplot.root";

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  // gStyle->SetPalette(55); //palettes from https://root.cern.ch/doc/master/classTColor.html
  TFile *f = TFile::Open(ifilename.c_str(),"r");
  TH1F *hcal_h, *hf_p_h, *hf_m_h;
  TString hN;
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  vector<pair<int,string>> typesOfBptx = {make_pair(0,"bptxxor_"),make_pair(1,"bptxplus_"),make_pair(2,"bptxminus_")};
  { //! For ordinary plots for CaloTower and ParticleFlow
    std::vector<std::tuple<TString, std::string, std::string>> histogramData = {
      std::make_tuple("CMScalo/calotower_e_eta_phi_minus_h[", "Energy at HF per tower, minus side", "calotower_e_eta_phi_minus_h.eps"),
      std::make_tuple("CMScalo/calotower_e_eta_phi_plus_h[", "Energy at HF per tower, plus side", "calotower_e_eta_phi_plus_h.eps"),
      std::make_tuple("CMScalo/calotower_e_eta_phi_h[", "Energy distribution per tower", "calotower_e_eta_phi_h.eps")
    };
    for(auto p: typesOfBptx){
      for(auto t: histogramData){
        TString histoName=std::get<0>(t);
        std::string name = std::get<1>(t);
        hN = histoName; hN+=p.first; hN+="]";
        TH2F *e_eta_phi_h = (TH2F*) f->Get(hN.Data());
        e_eta_phi_h->SetTitle((name+additionalText).c_str());

        // double sum = e_eta_phi_h->Integral(1,e_eta_phi_h->GetNbinsX(),1,e_eta_phi_h->GetNbinsY());
        // e_eta_phi_h->Scale(1.0/sum);
      
        // c1->SetLogy();
        // c1->SetLogz();
        c1->SetRightMargin(0.1);
        e_eta_phi_h->Draw("COLZ");
        c1->SaveAs((p.second+std::get<2>(t)).c_str());
        c1->Clear();
      }
    }
  }
}
