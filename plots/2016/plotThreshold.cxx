// #include <string>
void plotThreshold(std::string additionalText = "", std::string ifilename = "uaplot.root")
{
  // string additionalText=""; string ifilename = "uaplot.root";

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open(ifilename.c_str(),"r");
  TH1F *hcal_h, *hf_p_h, *hf_m_h;
  TString hN;
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  vector<pair<int,string>> typesOfBptx = { make_pair(0,"bptxxor_"),make_pair(1,"bptxplus_"),make_pair(2,"bptxminus_")};
  { //! For ordinary plots for CaloTower and ParticleFlow
    std::vector<std::tuple<TString, std::string, std::string>> histogramData = {
      std::make_tuple("CMScalo/calotower_e_eta_h[", "Calo Towers", "summary_calotower_e_eta.eps"),
      std::make_tuple("CMScalo/calotower_eMaxTower_eta_h[", "Calo Towers", "perTower_calotower_e_eta.eps"),
      std::make_tuple("CMSpf/pf_e_eta_h[", "Particle Flow", "summary_pf_e_eta.eps"),
      std::make_tuple("CMSpf/pf_eMaxTower_eta_h[", "Particle Flow", "perTower_pf_e_eta.eps")
    };
    for(auto p: typesOfBptx){
      for(auto t: histogramData){
        TString histoName=std::get<0>(t);
        std::string name = std::get<1>(t);
        hN = histoName; hN+=p.first; hN+="]";
        TH2F *pf_e_eta_h = (TH2F*) f->Get(hN.Data());
        pf_e_eta_h->SetTitle((name+additionalText).c_str());
        c1->SetLogy(); c1->SetLogz();
        pf_e_eta_h->Draw("COLZ");
        c1->SaveAs((p.second+std::get<2>(t)).c_str());
        c1->Clear();
      }
    }
  }
  { //! For ParticleFlow summary
    std::vector<std::tuple<std::vector<TString>, std::vector<TString>, std::string>> histogramData = {
      std::make_tuple<std::vector<TString>, std::vector<TString>, std::string>(
    {"CMSpf/pf_e_eta_h[", "CMSpf/pf_e_hch_h[", "CMSpf/pf_e_h0_h[", "CMSpf/pf_e_em0_h["},
    {"Particle Flow", "Charged hadrons", "h0", "gamma"}, "summary_pf_e_summary.eps"
        ),
      std::make_tuple<std::vector<TString>, std::vector<TString>, std::string>(
    {"CMSpf/pf_eMaxTower_eta_h[", "CMSpf/pf_eMaxTower_hch_h[", "CMSpf/pf_eMaxTower_h0_h[", "CMSpf/pf_eMaxTower_em0_h["},
    {"Particle Flow", "Charged hadrons", "h0", "gamma"}, "perTower_pf_e_summary.eps"
        )
    };
    for(auto p: typesOfBptx){
      for(auto t: histogramData){
        std::vector<TString> histoNames=std::get<0>(t);
        std::vector<TString> names=std::get<1>(t);
        c1->SetLogy(); c1->SetLogz();
        c1->Divide(2,2);
        std::vector<TH2F*> hs; hs.resize(histoNames.size());
        for(auto i = 0; i<histoNames.size();++i){
          hN = histoNames.at(i); hN+=p.first; hN+="]";
          hs.at(i) = (TH2F*) f->Get(hN.Data());
          hs.at(i)->SetTitle((names.at(i)+additionalText).Data());
          auto pad = c1->cd(i+1);
          pad->SetLogy(); pad->SetLogz();
          hs.at(i)->Draw("COLZ");
          // hs.at(i)->Print();
        }
        c1->cd();
        c1->SaveAs((p.second+std::get<2>(t)).c_str());
        c1->Clear();
        histoNames.clear(); hs.clear();
      }
    }
  }
  const double PERCENT_OF_ENTRIES = 99.7 / 100.0;
  { //! For definition of thresholds
    std::vector<std::tuple<std::vector<TString>, std::vector<TString>, std::vector<int>, std::string, double>> histogramData = {
      std::make_tuple<std::vector<TString>, std::vector<TString>, std::vector<int>, std::string, double>(
    {"CMSpf/pf_e_eta_h", "CMScalo/calotower_e_eta_h"},
    {"Particle Flow", "Calo Towers"},
    {2, 4},
    "summary_thresolds_comparing.eps", 45
        ),
      std::make_tuple<std::vector<TString>, std::vector<TString>, std::vector<int>, std::string, double>(
    {"CMSpf/pf_eMaxTower_eta_h", "CMScalo/calotower_eMaxTower_eta_h"},
    {"Particle Flow", "Calo Towers"},
    {2, 4},
    "perTower_thresolds_comparing.eps", 30
        )
    };
    for(auto p: typesOfBptx){
      for(auto t: histogramData){
        std::vector<TString> histoNames = std::get<0>(t);
        std::vector<TString> namesForHisto = std::get<1>(t);
        std::vector<int>     colors =        std::get<2>(t);
        std::vector<std::vector<double>> thresholds; // Matrix of thresholds for compared print
        c1->SetLogy();
        { //! Calculating of thresholds and printing it
          for (int i = 0; i<histoNames.size(); ++i){
            TString h = histoNames.at(i);
            hN = h; hN+="["; hN+=p.first; hN+="]";
            auto hs = (TH2F*) f->Get(hN.Data());
            std::vector<double> thresholds_current;
            auto nbins = hs->GetNbinsX();
            TString h_name = h.Remove(0,h.Index("/")+1); // Name without folder
            for(int b = 1; b<=nbins; b++){ // Finding threshold for bin
              gStyle->SetOptTitle(0);
              TString title_main=h_name; title_main+="_bin_"; title_main+=b<10?"0":""; title_main+=b;
              TCanvas *c2 = new TCanvas(title_main.Data(), "", 800, 200); // Canvas for printing of projections by bin
              TString ptText = "#eta="; ptText+=hs->GetXaxis()->GetBinLowEdge(b);
              TPaveText *pt = new TPaveText(.7,.6,.85,.8,"NDC"); pt->AddText(ptText.Data()); pt->SetFillColor(kWhite);
              TString title_full=title_main; title_full+="; E_{towers}";
              TH1D* projection = hs->ProjectionY(title_full.Data(),b,b);
              int sumI = projection->Integral(1, projection->GetNbinsX());
              double threshold; auto currentBin = 0;
              double sum = 0;
              if(sumI == 0){
                threshold=0;
              } else {
                projection->Scale(1.0/sumI);
                while ( sum < PERCENT_OF_ENTRIES ){
                  sum = projection->Integral(1, ++currentBin);
                }
                // std::cout << h_name.Data() << ":" << b << ", bin " << currentBin << ": " << sum << std::endl;
                // std::cout << "::4\n";
                threshold = projection->GetBinLowEdge(currentBin) + projection->GetBinWidth(currentBin);
              }
              std::cout << "For histogram: " << hN << ", bin: " << b << " threshold at " << threshold << std::endl;
              //std::cout <<" [ " << sum << " at bin " << currentBin << "; "<< sumI <<" events at bin ]" << std::endl;
              // std::printf("For histogram: %s, bin: %i threshold at %f [ %f at bin %i]\n", h.Data(), b, threshold, sum, currentBin);
              // std::printf("{ %f = %f + %i * %f }\n", threshold, projection->GetBinLowEdge(1), currentBin, projection->GetBinWidth(currentBin));
              thresholds_current.push_back(threshold);
              c2->SetLogy(); c2->SetLogx();
              projection->SetAxisRange(0.1, 1000);
              // projection->SetAxisRange(0, 1, "Y");
              projection->SetTitle(title_main);
              // projection->Print();
              projection->Draw();
              pt->Draw("SAME");
              if (sumI!=0) {
                double widthArray[3] = {0.1,threshold,2000};
                TH1F* th = new TH1F("","",2,widthArray); th->SetBinContent(1,1);
                // th->SetFillColorAlpha(kYellow,.35);
                th->SetFillColor(kYellow);
                th->SetFillStyle(3002);
                th->SetLineWidth(0);
                th->Draw("SAME");
              }
              projection->Draw("SAME");
              // projection->GetYaxis()->DrawClone("SAME");
              TString filename = p.second; filename+=title_main; filename+=".eps";
              c2->SaveAs(filename.Data());
              c2->Clear();
              gStyle->SetOptTitle(1);
            }
            thresholds.push_back(thresholds_current);
            thresholds_current.clear();
          }
        }
        { //! printing comparing CaloTowers and PF
          TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
          THStack *stack = new THStack("Comparing thresholds",
                                       "Comparing thresholds from CaloTower and ParticleFlow"); //THStack is as 'SAME' draw option, but better
          TLegend* leg = new TLegend(0.13,0.79,0.33,0.89);
          for (int i = 0; i<histoNames.size(); ++i){
            hN = histoNames.at(i); hN+="["; hN+=p.first; hN+="]";
            auto hs = (TH2F*) f->Get(hN.Data());
            auto nbins = hs->GetNbinsX();
            TH1F* threshold_h = new TH1F("threshold_h",
                                         (namesForHisto.at(i)+";#eta;E_{thr} [GeV]").Data(),
                                         nbins,
                                         hs->GetXaxis()->GetBinLowEdge(1),
                                         hs->GetXaxis()->GetBinLowEdge(nbins) + hs->GetXaxis()->GetBinWidth(nbins));
            for(int bin = 0; bin<nbins; ++bin) threshold_h->SetBinContent(bin+1, thresholds.at(i).at(bin));
            std::cout << "Thresholds for " << namesForHisto.at(i).Data() <<": ";
            for(int bin = 0; bin<nbins; ++bin) std::cout<< threshold_h->GetBinContent(bin+1) << " ";std::cout << std::endl;
            threshold_h->SetName((namesForHisto.at(i)+additionalText).Data());
            threshold_h->SetXTitle("#eta");
            threshold_h->SetYTitle("E_{thr} [GeV]");
            threshold_h->SetLineColor(colors.at(i));
            threshold_h->SetLineWidth(2);
            threshold_h->SetAxisRange(0, std::get<4>(t),"Y");
            // threshold_h->Draw("H SAME");
            stack->Add(threshold_h);
            leg->AddEntry(threshold_h, (namesForHisto.at(i)+additionalText).Data());
            // TPaveText *pt = new TPaveText(.05,.1,.95,.8);
          }
          string name="Comparing thresholds from CaloTowers and ParticleFlow"; name+=additionalText;
          stack->SetTitle(name.c_str()); c3->SetTitle(name.c_str()); c3->SetName(name.c_str());
          stack->Draw("nostack");
          stack->GetXaxis()->SetTitle("#eta");
          stack->GetYaxis()->SetTitle("E_{thr} [GeV]");
          leg->Draw();
          c3->Update();
          c3->SaveAs((p.second+std::get<3>(t)).c_str());
          c3->Clear();
        }
      }
    }
  }
  { //! For caloTower threshold for plus and minus side 
    std::vector<std::tuple<TString, std::string, std::string>> histogramData = {
      std::make_tuple("CMScalo/hf_max_towerE_minus_h[", "Calo Towers threshold for HF minus side", "hf_max_towerE_minus_h.eps"),
      std::make_tuple("CMScalo/hf_max_towerE_plus_h[", "Calo Towers threshold for HF plus side", "hf_max_towerE_plus_h.eps"),
      // std::make_tuple("CMScalo/hf_max_towerE_minus_2_h[", "Calo Towers threshold for HF minus side", "hf_max_towerE_minus_2_h.eps"),
      // std::make_tuple("CMScalo/hf_max_towerE_plus_2_h[", "Calo Towers threshold for HF plus side", "hf_max_towerE_plus_2_h.eps"),
      std::make_tuple("CMScalo/hf_max_towerE_minus_ill_h[", "Calo Towers threshold for ill tower on HF minus side", "hf_max_towerE_minus_ill_h.eps"),
      std::make_tuple("CMScalo/hf_max_towerE_plus_ill_h[", "Calo Towers threshold for ill tower on HF plus side", "hf_max_towerE_plus_ill_h.eps"),
      // std::make_tuple("CMScalo/hf_max_towerE_minus_2_full_h[", "Calo Towers threshold for HF minus side", "hf_max_towerE_minus_2_full_h.eps"),
      // std::make_tuple("CMScalo/hf_max_towerE_plus_2_full_h[", "Calo Towers threshold for HF plus side", "hf_max_towerE_plus_2_full_h.eps"),
    };
    for(auto p: typesOfBptx){
      for(auto t: histogramData){
        TString histoName=std::get<0>(t);
        hN = histoName; hN+=p.first; hN+="]";
        std::string name = std::get<1>(t);
        TCanvas *c4 = new TCanvas(name.c_str(), "", 800, 200); // Canvas for printing of projections by bin
        c4->SetLogy(); c4->SetLogx();
        TH1F *histo = (TH1F*) f->Get(hN.Data());
        histo->SetTitle((name+additionalText).c_str());
        int sumI = histo->Integral(1, histo->GetNbinsX());
        double threshold; auto currentBin = 0;
        double sum = 0;
        if(sumI == 0){
          threshold=0;
        } else {
          histo->Scale(1.0/sumI);
          while ( sum < PERCENT_OF_ENTRIES ){
            sum = histo->Integral(1, ++currentBin);
          }
          threshold = histo->GetBinLowEdge(currentBin) + histo->GetBinWidth(currentBin);
        }
        std::cout << "For histogram: " << name << " (" << hN << ")" << " threshold at " << threshold << std::endl;
        //std::cout <<" [ " << sum << " at bin " << currentBin << "; "<< sumI <<" events at bin ]" << std::endl;

        histo->SetAxisRange(0.5, 1000);
        histo->Draw("H");
        if (sumI!=0) {
          double widthArray[3] = {0.1,threshold,7000};
          TH1F* th = new TH1F("","",2,widthArray); th->SetBinContent(1,1);
          // th->SetFillColorAlpha(kYellow,.35);
          th->SetFillColor(kYellow);
          th->SetFillStyle(3002);
          th->SetLineWidth(0);
          th->Draw("SAME");
        }
        c4->SaveAs((p.second+std::get<2>(t)).c_str());
        c4->Clear();
      }
    }
  }
  f->Close();
}
