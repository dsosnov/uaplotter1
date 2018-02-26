#define COUT_VALUE(X) std::cout << #X << ": " << X << std::endl;

void plotnRGs()
{
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);

  auto epos_ppb = TFile::Open("uaplot_output_histos_EPOS_pPb.root","r");
  auto epos_pbp = TFile::Open("uaplot_output_histos_EPOS_Pbp.root","r");
  auto data_ppb = TFile::Open("uaplot_output_histos_data_PARun2016_pPb.root","r");
  auto data_pbp = TFile::Open("uaplot_output_histos_data_PARun2016_Pbp.root","r");
  auto hijing_ppb = TFile::Open("uaplot_output_histos_Hijing_pPb.root","r");
  auto hijing_pbp = TFile::Open("uaplot_output_histos_Hijing_Pbp.root","r");
  auto qgsjet_ppb = TFile::Open("uaplot_output_histos_QGSJET_pPb.root","r");
  auto qgsjet_pbp = TFile::Open("uaplot_output_histos_QGSJET_Pbp.root","r");

  std::vector<std::tuple<TFile*, std::string, int, int, int>> roots_ppb = {
    std::make_tuple(epos_ppb,  "EPOS, pPb", 3, 3006, 22),
    std::make_tuple(hijing_ppb,"Hijing, pPb", 4, 3007, 23),
    std::make_tuple(qgsjet_ppb,  "QGSJet-II, pPb", 8, 3006, 33),
    std::make_tuple(data_ppb,  "Data, pPb", 2, 3004, 20)
  };
  std::vector<std::tuple<TFile*, std::string, int, int, int>> roots_pbp = {
    std::make_tuple(epos_pbp,  "EPOS, Pbp", 3, 3006, 22),
    std::make_tuple(hijing_pbp,"Hijing, Pbp", 4, 3007, 23),
    std::make_tuple(qgsjet_pbp,  "QGSJet-II, Pbp", 8, 3006, 33),
    std::make_tuple(data_pbp,  "Data, Pbp", 2, 3004, 20)
  };
  std::vector<std::tuple< std::vector<std::tuple<TFile*, std::string, int, int, int>>, std::string  >> roots = {
    std::make_tuple(roots_ppb, "ppb"),
    std::make_tuple(roots_pbp, "pbp")
  };

  if(true){ // Empty bins distribituon
    vector<tuple<string, string>> histograms = {
      make_tuple("CMSpf/activity_loose_h", "pf_emptyBins_loose"),
      make_tuple("CMSpf/activity_tight_h", "pf_emptyBins_tight"),
      make_tuple("CMStracking/activity_tight_h", "tracking_emptyBins_tight"),
      make_tuple("CMStracking/activity_loose_h", "tracking_emptyBins_loose"),
      make_tuple("central_activity_h", "emptyBins")
    };
    for(auto h: histograms){
      for(auto rr: roots){
        string outFile = get<1>(h) + string("_") + std::get<1>(rr);
        TCanvas *c = new TCanvas(outFile.c_str(), "", 800, 600);
        auto th = new THStack(outFile.c_str(), "; #eta");
        auto leg = new TLegend(0.69, 0.69, 0.99, 0.92);
        for (auto r: std::get<0>(rr)){
          auto hist = (TH1F*)(std::get<0>(r))->Get( (get<0>(h) + string("[0]")).c_str());
          // hist->Print();
          hist->SetTitle(std::get<1>(r).c_str());
          hist->SetLineColor(std::get<2>(r));
          hist->SetLineWidth(3);
          hist->SetMarkerColor(std::get<2>(r));
          hist->SetMarkerStyle(std::get<4>(r));
          // hist->SetFillColor(std::get<2>(r));
          // hist->SetFillStyle(std::get<3>(r));
          auto countOfEvents = ((TH1F*)(std::get<0>(r))->Get("n_sd_minus_bins_h[0]"))->GetEntries();
          // auto countOfEvents = ((TH1F*)(std::get<0>(r))->Get("CMSCalo/hf_etotal_minus_h[0]"))->GetEntries();
          // cout << "count of events: " << countOfEvents << endl;
          for(int i = 1; i <= hist->GetNbinsX(); i++){ //Setting x axis to eta range
            int n = hist->GetBinContent(i);
            if(n!=0) hist->SetBinContent(i, countOfEvents-n);
          }
          hist->Sumw2(); hist->Scale(1.0/countOfEvents);
          // for(int i = 1; i <= hist->GetNbinsX(); i++) COUT_VALUE(hist->GetBinError(i));
          th->Add(hist); leg->AddEntry(hist);
        }
        th->Draw("nostack,e");
        leg->Draw();
        c->SetLogy(false);
        c->SaveAs((outFile+string("_lin")+string(".eps")).c_str());
        // c->SetLogy();
        // c->SaveAs((outFile+string("_log")+string(".eps")).c_str());
        c->Clear();
      }
    }
  }
  if(true){ // Distribution of \Delta_{\eta^F} minus/plus.
    for(auto rr: roots){
      std::vector<std::tuple<std::string, std::string, int>> histograms = {
        std::make_tuple("n_sd_minus_bins_h", "minusBins", 10),
        std::make_tuple("n_sd_plus_bins_h",  "plusBins" , 9 )
      };
      for(auto h: histograms){
        string outFile = std::get<1>(h) + string("_") + std::get<1>(rr);
        TCanvas *c = new TCanvas(outFile.c_str(), "", 600, 600);
        auto th = new THStack(std::get<1>(h).c_str(), "; |#Delta_{#eta^F}|");
        auto leg = new TLegend(0.69, 0.69, 0.99, 0.92);
        for (auto r: std::get<0>(rr)){
          auto histname = std::get<0>(h) + string("[")+to_string(std::get<2>(h))+string("]");
          auto h0 = (TH1F*)((std::get<0>(r))->Get(histname.c_str())->Clone());
          auto hist = new TH1F(std::get<1>(r).c_str(),std::get<1>(r).c_str(),13, 0, 6.50);
          // h0->Sumw2(); h0->Scale(1.0/h0->GetEntries());
          h0->Sumw2();
          for(int i = 2; i < 15; i++){ //Setting x axis to eta range
            hist->SetBinContent(i-1, h0->GetBinContent(i));
            hist->SetBinError(i-1, h0->GetBinError(i));
            // cout << "hist " << i << " " << h0->GetBinContent(i-1) << " " << h0->GetBinError(i-1) << endl;
          }
          // auto hist = (TH1F*)(std::get<0>(r))->Get((std::get<0>(h) + string("[")+to_string(std::get<2>(h))+string("]")).c_str());
          hist->Print();
          hist->SetTitle(std::get<1>(r).c_str());
          hist->SetLineColor(std::get<2>(r));
          hist->SetLineWidth(3);
          hist->SetMarkerColor(std::get<2>(r));
          hist->SetMarkerStyle(std::get<4>(r));
          // hist->SetFillColor(std::get<2>(r));
          // hist->SetFillStyle(std::get<3>(r));
          // for(int i = 1; i <= hist->GetNbinsX(); i++){ //Setting x axis to eta range
          //   hist->GetXaxis()->SetBinLabel(i, ((i%2)?string(""):to_string(i/2-1)).c_str());
          // }
          // hist->Sumw2(); hist->Scale(1.0/hist->GetEntries());
          // hist->Sumw2();
          hist->Scale(1.0/hist->Integral(1, hist->GetNbinsX()));
          hist->Print();
          th->Add(hist);
          leg->AddEntry(hist);
        }
        th->Draw("nostack,e");
        // th->GetXaxis()->SetLimits(0,13);
        th->SetMaximum(1);
        leg->Draw();
        c->SetLogy();
        c->SaveAs((outFile+string("_log")+string(".eps")).c_str());
        c->SetLogy(false);
        th->SetMaximum(0.03);
        th->Draw("nostack,e");
        c->SaveAs((outFile+string("_lin")+string(".eps")).c_str());
        c->Clear();
      }
    }
  }
  if(true){ // Distribution of \Delta_{\eta^F} minus/plus per PID. //TODO
    std::vector<std::tuple<std::string, std::string, int>> histograms = {
      std::make_tuple("n_sd_minus_bins_h", "minusBins_byPids", 25), //25..29 is "on right side > 4"
      std::make_tuple("n_sd_plus_bins_h",  "plusBins_byPids" , 20)  //20..25 is "on left side > 4"
    };
    std::vector<std::tuple<std::string, int, std::string, int>> pids = {
      std::make_tuple("ND",         3, "ND", 1),
      // std::make_tuple("CD",         7, "ND && CD", 1),
      // std::make_tuple("SD-[RG@-Z]", 4, "ND && CD && SD-", 1),
      // std::make_tuple("SD+[RG@+Z]", 2, "ND && CD && SD- && SD+", 1),
      std::make_tuple("SD & CD",   4, "CD && SD- && SD+", 3),
      std::make_tuple("DD",         6, "ND && CD && SD- && SD+ && DD", 1)
    };
    // std::vector<std::tuple<std::string, int, std::string>> pids = {
    //   std::make_tuple("ND",         3, "ND"),
    //   std::make_tuple("CD",         6, "ND && CD"),
    //   std::make_tuple("SD-[RG@-Z]", 4, "ND && CD && SD-"),
    //   std::make_tuple("SD+[RG@+Z]", 2, "ND && CD && SD- && SD+"),
    //   std::make_tuple("DD",         7, "ND && CD && SD- && SD+ && DD")
    // };
    std::vector<std::tuple< TFile*, std::string>> files_with_pids = {
      std::make_tuple(epos_ppb, "epos_ppb"),
      std::make_tuple(epos_pbp, "epos_pbp")
    };
    for(auto r: files_with_pids){
      for(auto h: histograms){
        string outFile = std::get<1>(r) + string("_") + std::get<1>(h);
        auto *c = new TCanvas(outFile.c_str(), "", 600, 600);
        auto th = new THStack(std::get<1>(h).c_str(), "; |#Delta_{#eta^F}|");
        auto leg = new TLegend(0.69, 0.69, 0.99, 0.92);
        auto countOfEvents = ((TH1F*)(std::get<0>(r))->Get("CMScalo/hf_etotal_minus_h[0]"))->GetEntries();
        auto currentHistN = 0;
        for (int i = 0; i < pids.size(); i++){
          TH1F* h0;
          for(auto j = 0; j < std::get<3>(pids.at(i)); j++){
            auto histname = std::get<0>(h) + string("[") + to_string(std::get<2>(h)+currentHistN) + string("]");
          COUT_VALUE(outFile);
          COUT_VALUE(histname);
            currentHistN++;
            auto h00  = (TH1F*)((std::get<0>(r))->Get(histname.c_str())->Clone());
            if(j==0){
              h0 = h00;
            }else{
              h0->Add(h00);
            }
          }
          // auto h0 = (TH1F*)((std::get<0>(r))->Get(histname.c_str())->Clone());
          auto hist = new TH1F(std::get<1>(r).c_str(),std::get<1>(r).c_str(),13, 0, 6.50);
          // h0->Sumw2(); h0->Scale(1.0/h0->GetEntries());
          h0->Sumw2();
          for(int i = 2; i < 15; i++){ //Setting x axis to eta range
            hist->SetBinContent(i-1, h0->GetBinContent(i));
            hist->SetBinError(i-1, h0->GetBinError(i));
            // cout << "hist " << i << " " << hist->GetBinContent(i-1) << " " << hist->GetBinError(i-1) << endl;
          }
          hist->Print();
          hist->SetTitle(std::get<0>(pids.at(i)).c_str());
          // hist->SetTitle(std::get<2>(pids.at(i)).c_str());
          hist->SetLineColor(std::get<1>(pids.at(i)));
          hist->SetLineWidth(3);
          hist->SetMarkerColor(std::get<1>(pids.at(i)));
          // hist->SetMarkerStyle(std::get<4>(pids.at(i)));
          // hist->SetFillColor(std::get<1>(pids.at(i)));
          // hist->SetFillStyle(3004);
          // hist->Sumw2(); hist->Scale(1.0/countOfEvents);
          hist->Scale(1.0/countOfEvents);
          // for(int i = 1; i <= hist->GetNbinsX(); i++){ //Setting x axis to eta range
          //  hist->GetXaxis()->SetBinLabel(i, ((i%2)?string(""):to_string(i/2-1)).c_str());
          // }
          hist->Print();
          th->Add(hist); leg->AddEntry(hist);
        }
        th->Draw("nostack,e");
        leg->Draw();
        th->SetMinimum(1E-6);
        th->SetMaximum(1);
        c->SetLogy();
        c->SaveAs((outFile+string("_log")+string(".eps")).c_str());
        c->SetLogy(false);
        th->SetMaximum(0.3);
        c->SaveAs((outFile+string("_lin")+string(".eps")).c_str());
        c->Clear();
      }
    }
  }
  if(true){ // Distribution of \Delta_{\eta^F} minus/plus per PID. -- QGSJet //TODO
    std::vector<std::tuple<std::string, std::string, int>> histograms = {
      std::make_tuple("n_sd_minus_bins_h", "minusBins_byPids", 25), //25..29 is "on right side > 4"
      std::make_tuple("n_sd_plus_bins_h",  "plusBins_byPids" , 20)  //20..25 is "on left side > 4"
    };
    std::vector<std::tuple<std::string, int, std::string, int>> pids = {
      std::make_tuple("ND",         3, "ND", 1),
      // std::make_tuple("CD",         7, "ND && CD", 1),
      // std::make_tuple("SD-[RG@-Z]", 4, "ND && CD && SD-", 1),
      // std::make_tuple("SD+[RG@+Z]", 2, "ND && CD && SD- && SD+", 1),
      std::make_tuple("SD & CD",   4, "CD && SD- && SD+", 3),
      std::make_tuple("DD",         6, "ND && CD && SD- && SD+ && DD", 1)
    };
    // std::vector<std::tuple<std::string, int, std::string>> pids = {
    //   std::make_tuple("ND",         3, "ND"),
    //   std::make_tuple("CD",         6, "ND && CD"),
    //   std::make_tuple("SD-[RG@-Z]", 4, "ND && CD && SD-"),
    //   std::make_tuple("SD+[RG@+Z]", 2, "ND && CD && SD- && SD+"),
    //   std::make_tuple("DD",         7, "ND && CD && SD- && SD+ && DD")
    // };
    std::vector<std::tuple< TFile*, std::string>> files_with_pids = {
      std::make_tuple(qgsjet_ppb, "qgsjet_ppb"),
      std::make_tuple(qgsjet_pbp, "qgsjet_pbp")
    };
    for(auto r: files_with_pids){
      for(auto h: histograms){
        string outFile = std::get<1>(r) + string("_") + std::get<1>(h);
        auto *c = new TCanvas(outFile.c_str(), "", 600, 600);
        auto th = new THStack(std::get<1>(h).c_str(), "; |#Delta_{#eta^F}|");
        auto leg = new TLegend(0.69, 0.69, 0.99, 0.92);
        auto countOfEvents = ((TH1F*)(std::get<0>(r))->Get("CMScalo/hf_etotal_minus_h[0]"))->GetEntries();
        auto currentHistN = 0;
        for (int i = 0; i < pids.size(); i++){
          TH1F* h0;
          for(auto j = 0; j < std::get<3>(pids.at(i)); j++){
            auto histname = std::get<0>(h) + string("[") + to_string(std::get<2>(h)+currentHistN) + string("]");
          COUT_VALUE(outFile);
          COUT_VALUE(histname);
            currentHistN++;
            auto h00  = (TH1F*)((std::get<0>(r))->Get(histname.c_str())->Clone());
            if(j==0){
              h0 = h00;
            }else{
              h0->Add(h00);
            }
          }
          // auto h0 = (TH1F*)((std::get<0>(r))->Get(histname.c_str())->Clone());
          auto hist = new TH1F(std::get<1>(r).c_str(),std::get<1>(r).c_str(),13, 0, 6.50);
          // h0->Sumw2(); h0->Scale(1.0/h0->GetEntries());
          h0->Sumw2();
          for(int i = 2; i < 15; i++){ //Setting x axis to eta range
            hist->SetBinContent(i-1, h0->GetBinContent(i));
            hist->SetBinError(i-1, h0->GetBinError(i));
            // cout << "hist " << i << " " << hist->GetBinContent(i-1) << " " << hist->GetBinError(i-1) << endl;
          }
          hist->Print();
          hist->SetTitle(std::get<0>(pids.at(i)).c_str());
          // hist->SetTitle(std::get<2>(pids.at(i)).c_str());
          hist->SetLineColor(std::get<1>(pids.at(i)));
          hist->SetLineWidth(3);
          hist->SetMarkerColor(std::get<1>(pids.at(i)));
          // hist->SetMarkerStyle(std::get<4>(pids.at(i)));
          // hist->SetFillColor(std::get<1>(pids.at(i)));
          // hist->SetFillStyle(3004);
          // hist->Sumw2(); hist->Scale(1.0/countOfEvents);
          hist->Scale(1.0/countOfEvents);
          // for(int i = 1; i <= hist->GetNbinsX(); i++){ //Setting x axis to eta range
          //  hist->GetXaxis()->SetBinLabel(i, ((i%2)?string(""):to_string(i/2-1)).c_str());
          // }
          hist->Print();
          th->Add(hist); leg->AddEntry(hist);
        }
        th->Draw("nostack,e");
        leg->Draw();
        th->SetMinimum(1E-6);
        th->SetMaximum(1);
        c->SetLogy();
        c->SaveAs((outFile+string("_log")+string(".eps")).c_str());
        c->SetLogy(false);
        th->SetMaximum(0.3);
        c->SaveAs((outFile+string("_lin")+string(".eps")).c_str());
        c->Clear();
      }
    }
  }
  if(true){ // HF for SDs
    for(auto rr: roots){
      vector<tuple<tuple<string, string>, tuple<string, string>, string>> hist_types = {
        make_tuple(make_tuple("CMScalo/hf_etotal_minus_h","hfm_sume"), make_tuple("CMScalo/hf_etotal_plus_h","hfp_sume"), "; E_{Sum}"),
        make_tuple(make_tuple("CMScalo/hf_max_towerE_minus_h","hfm_twre"), make_tuple("CMScalo/hf_max_towerE_plus_h","hfp_twre"), "; E_{Max_Tower}")
      };
      int countOfBins=7;
      int start_defm = 30+countOfBins*2; int start_defp = 30+countOfBins*3;
      vector<tuple<tuple<string,int>, tuple<string,int>, int, tuple<tuple<int, int>, tuple<int, int>>>> sd_types= {
        make_tuple(make_tuple("deFm00", 18            ), make_tuple("deFp00", 19            ), 1, make_tuple(make_tuple(5000, 50), make_tuple(200, 25))),
        make_tuple(make_tuple("deFm01", start_defm + 0), make_tuple("deFp01", start_defp + 0), 1, make_tuple(make_tuple(5000, 50), make_tuple(200, 25))),
        make_tuple(make_tuple("deFm12", start_defm + 1), make_tuple("deFp12", start_defp + 1), 1, make_tuple(make_tuple(500, 10), make_tuple(80, 20))),
        make_tuple(make_tuple("deFm23", start_defm + 2), make_tuple("deFp23", start_defp + 2), 1, make_tuple(make_tuple(350, 10), make_tuple(40, 20))),
        make_tuple(make_tuple("deFm34", start_defm + 3), make_tuple("deFp34", start_defp + 3), 1, make_tuple(make_tuple(350, 10), make_tuple(40, 20))),
        make_tuple(make_tuple("deFm45", start_defm + 4), make_tuple("deFp45", start_defp + 4), 1, make_tuple(make_tuple(200, 10), make_tuple(25, 10))),
        make_tuple(make_tuple("deFm56", start_defm + 5), make_tuple("deFp56", start_defp + 5), 1, make_tuple(make_tuple(200, 10), make_tuple(10, 10))),
        make_tuple(make_tuple("deFm67", start_defm + 6), make_tuple("deFp67", start_defp + 6), 1, make_tuple(make_tuple(200, 10), make_tuple( 6, 2))),
        make_tuple(make_tuple("deFm37", start_defm + 3), make_tuple("deFp37", start_defp + 3), 4, make_tuple(make_tuple(200, 5), make_tuple( 6, 2)))
      }; // ( (minus_name, minus_cut), (plus_name, plus_cut), count_plots_merged, ((maxx_sume, rebin_sume), (maxx_twre, rebin_twre)) )
      for(auto hn = 0; hn < hist_types.size(); hn++){
        auto h = hist_types.at(hn);
        for(auto s: sd_types){
          for(auto j = 0; j < 2; j++){
            auto current_cut  = (j==1) ? get<1>(s) : get<0>(s);
            auto current_hist = (j==1) ? get<1>(h) : get<0>(h);
            auto hist_params = get<3>(s);
            auto hist_param_pair = (hn==1) ? get<1>(hist_params) : get<0>(hist_params);
            auto maxx = get<0>(hist_param_pair);
            auto rebin = get<1>(hist_param_pair);
            string outFile = std::get<0>(current_cut) + string("_") + std::get<1>(current_hist) + string("_") + std::get<1>(rr);
            TCanvas *c = new TCanvas(outFile.c_str(), "", 600, 600);
            auto th = new THStack(outFile.c_str(), get<2>(h).c_str());
            auto leg = new TLegend(0.69, 0.69, 0.99, 0.92);
            for (auto r: std::get<0>(rr)){
              TH1F* hist;
              for(auto i = 0; i < std::get<2>(s); i++){
                // cout << "outFile: " << outFile << endl; 
                string histname = std::get<0>(current_hist) + string("[")+std::to_string(i+get<1>(current_cut))+string("]");
                auto h0  = (TH1F*)((std::get<0>(r))->Get(histname.c_str())->Clone());
                if(i==0){
                  hist = h0;
                }else{
                  hist->Add(h0);
                }
              }

              hist->SetTitle(std::get<1>(r).c_str());
              hist->SetLineColor(std::get<2>(r));
              hist->SetLineWidth(3);
              hist->SetMarkerColor(std::get<2>(r));
              hist->SetMarkerStyle(std::get<4>(r));
              // hist->SetFillColor(std::get<2>(r));
              // hist->SetFillStyle(std::get<3>(r));
              if(rebin>1) hist->Rebin(rebin);
              // hist->Sumw2(); hist->Scale(1.0/hist->Integral(1,hist->GetNbinsX()));
              hist->Sumw2(); hist->Scale(1.0/hist->Integral(0,hist->GetNbinsX()+1));
              // hist->Sumw2(); hist->Scale(1.0/hist->GetEntries());
              th->Add(hist); leg->AddEntry(hist);
            }
            th->Draw("nostack,e");
            leg->Draw();
            th->GetXaxis()->SetLimits(0.8, 7000);
            th->SetMinimum(5E-5);
            c->SetLogx(); c->SetLogy();
            c->SaveAs((outFile+string("_log")+string(".eps")).c_str());
            c->SetLogx(false); c->SetLogy(false);
            if(maxx > 0) th->GetXaxis()->SetLimits(0, maxx);
            c->SaveAs((outFile+string("_lin")+string(".eps")).c_str());
            c->Clear();
          }
        }
      }
    }
  }
  if(true){ // PIDs
    for(auto rr: roots){
      string outFile = string("pids_") + std::get<1>(rr);
      cout << "outFile: " << outFile << endl;
      TCanvas *c = new TCanvas(outFile.c_str(), "", 800, 600);
      auto th = new THStack("PIDs", "");
      auto leg = new TLegend(0.74, 0.74, 0.99, 0.92);
      for (auto r: std::get<0>(rr)){
        string histname = "sd_flag_central_reco_h[0]";
        auto hist = (TH1F*)(std::get<0>(r))->Get(histname.c_str());
        // hist->Print();
        hist->SetTitle(std::get<1>(r).c_str());
        hist->SetLineColor(std::get<2>(r));
        hist->SetLineWidth(3);
        hist->SetMarkerColor(std::get<2>(r));
        hist->SetMarkerStyle(std::get<4>(r));
        // hist->SetFillColor(std::get<2>(r));
        // hist->SetFillStyle(std::get<3>(r));
        auto countOfEvents = ((TH1F*)(std::get<0>(r))->Get("n_sd_minus_bins_h[0]"))->GetEntries();
        // auto countOfEvents = ((TH1F*)(std::get<0>(r))->Get("CMSCalo/hf_etotal_minus_h[0]"))->GetEntries();
        // cout << "count of events: " << countOfEvents << endl;
        auto a = hist->GetXaxis();
        vector<string> l = {"","SD-", "ND", "SD+", "DD", "CD", "Elastic", ""};
        for(auto i = 0; i < l.size(); i++) a->SetBinLabel(i+1, l.at(i).c_str());
        // ChangeLabel(1,-1,-1,-1,2);
        hist->Sumw2(); hist->Scale(1.0/countOfEvents);
        th->Add(hist); leg->AddEntry(hist);
      }
      c->SetLogy();
      th->Draw("nostack,e");
      leg->Draw();
      c->SaveAs((outFile+string("")+string(".eps")).c_str());
      c->Clear();
    }
  }
  if(false){ // E_Max_Tower on sides
    vector<tuple<string, int>> types = {make_tuple("all_", 17), make_tuple("inelastic_", 0)};
    for(auto t: types){
      for(auto rr: roots){
        string outFile = string("towers_mp_") + get<0>(t) + std::get<1>(rr) + string(".eps");
        cout << "outFile: " << outFile << endl;
        TCanvas *c = new TCanvas(outFile.c_str(), "", 800, 600);
        auto th = new THStack("towers_mp_", "; E_{TowerHF-}; E_{TowerHF+}");
        auto leg = new TLegend(0.74, 0.74, 0.99, 0.92);
        for (auto r: std::get<0>(rr)){
          string histname = string("CMScalo/hf_maxTowerE_minus_plus_h[")+to_string(get<1>(t)) + string("]");
          auto hist = (TH2F*)(std::get<0>(r))->Get(histname.c_str());
          // hist->Print();
          hist->SetTitle(std::get<1>(r).c_str());
          hist->SetLineColor(std::get<2>(r));
          hist->SetLineWidth(3);
          // hist->SetMarkerColor(std::get<2>(r));
          // hist->SetMarkerStyle(std::get<4>(pids.at(i)));
          // hist->SetFillColor(std::get<2>(r));
          // hist->SetFillStyle(std::get<3>(r));
          auto countOfEvents = hist->GetEntries();
          COUT_VALUE(histname);
          COUT_VALUE(std::get<1>(r));
          COUT_VALUE(countOfEvents);
          COUT_VALUE(hist->Integral(1, 8, 1, 8));
          hist->Sumw2(); hist->Scale(1.0/countOfEvents);
          COUT_VALUE(hist->Integral(1, 8, 1, 8));
          th->Add(hist); leg->AddEntry(hist);
        }
        th->Draw("nostack,box");
        leg->Draw();
        c->SaveAs(outFile.c_str());
        c->Clear();
      }
    }
  }
  if(true){ // HF for elastic with comparing with BPTX_XOR 
    vector<tuple<tuple<string, string>, tuple<string, string>, string, int, int>> hist_types = {
      make_tuple(make_tuple("CMScalo/hf_etotal_minus_h","hfm_sume"), make_tuple("CMScalo/hf_etotal_plus_h","hfp_sume"), "; E_{Sum}", 200, 2),
      make_tuple(make_tuple("CMScalo/hf_max_towerE_minus_h","hfm_twre"), make_tuple("CMScalo/hf_max_towerE_plus_h","hfp_twre"), "; E_{Max_Tower}", 6, 2)
    };
    auto corename = "elastic";
    auto noise_pbp = TFile::Open("../bptxor/bptxxor_pbp.root", "r");
    auto noise_ppb = TFile::Open("../bptxor/bptxxor_ppb.root", "r");
    vector<tuple<TFile*, string, int, int, int, int>> hist_files_ppb = {
      make_tuple(epos_ppb, "EPOS, pPb", 3, 2, 3, 22),
      make_tuple(hijing_ppb, "HIJING, pPb", 4, 2, 3, 23),
      make_tuple(data_ppb, "DATA, pPb", 2, 2, 3, 20),
      make_tuple(noise_ppb, "Noise, pPb", 6, 5, 6, 21)
    };
    vector<tuple<TFile*, string, int, int, int, int>> hist_files_pbp = {
      make_tuple(epos_pbp, "EPOS, Pbp", 3, 2, 3, 22),
      make_tuple(hijing_pbp, "HIJING, Pbp", 4, 2, 3, 23),
      make_tuple(data_pbp, "DATA, Pbp", 2, 2, 3, 20),
      make_tuple(noise_pbp, "Noise, Pbp", 6, 5, 6, 21)
    };
    vector<tuple<vector<tuple<TFile*, string, int, int, int, int>>, string>> hist_files = {
      make_tuple(hist_files_ppb, "ppb"),
      make_tuple(hist_files_pbp, "pbp")
    };
    for(auto rr: hist_files){
      for(auto hn = 0; hn < hist_types.size(); hn++){
        auto h = hist_types.at(hn);
        for(auto j = 0; j < 2; j++){
            auto current_cut  = corename;
            auto current_hist = (j==1) ? get<1>(h) : get<0>(h);
            auto maxx = get<3>(h);
            auto rebin = get<4>(h);
            string outFile = current_cut + string("_") + std::get<1>(current_hist) + string("_") + std::get<1>(rr);
            TCanvas *c = new TCanvas(outFile.c_str(), "", 600, 600);
            auto th = new THStack(outFile.c_str(), get<2>(h).c_str());
            auto leg = new TLegend(0.69, 0.69, 0.99, 0.92);
            for (auto r: std::get<0>(rr)){
              auto histn = (j==1) ? get<4>(r) : get<3>(r);
              string histname = std::get<0>(current_hist) + string("[")+std::to_string(histn)+string("]");
              auto hist  = (TH1F*)((std::get<0>(r))->Get(histname.c_str())->Clone());
              hist->Print();
              hist->SetTitle(std::get<1>(r).c_str());
              hist->SetLineColor(std::get<2>(r));
              hist->SetLineWidth(3);
              hist->SetMarkerColor(std::get<2>(r));
              hist->SetMarkerStyle(std::get<5>(r));
              // hist->SetFillColor(std::get<2>(r));
              // hist->SetFillStyle(std::get<3>(r));
              if(rebin>1) hist->Rebin(rebin);
              // hist->Sumw2(); hist->Scale(1.0/hist->Integral(1,hist->GetNbinsX()));
              hist->Sumw2(); hist->Scale(1.0/hist->Integral(0,hist->GetNbinsX()+1));
              // hist->Sumw2(); hist->Scale(1.0/hist->GetEntries());
              th->Add(hist); leg->AddEntry(hist);
            }
            th->Draw("nostack,e");
            leg->Draw();
            c->SetLogx(); c->SetLogy();
            th->GetXaxis()->SetLimits(0.8, 7000);
            th->SetMinimum(5E-5);
            c->SaveAs((outFile+string("_log")+string(".eps")).c_str());
            c->SetLogx(false); c->SetLogy(false);
            if(maxx > 0) th->GetXaxis()->SetLimits(0, maxx);
            c->SaveAs((outFile+string("_lin")+string(".eps")).c_str());
            c->Clear();
        }
      }
    }
  }

}
