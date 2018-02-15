#include <string>

#define DEBUG

double GetErrorOfDivision(double A, double dA, double B, double dB){
  return (dA - A/B*dB)/B;
}

void PrintCounts(string file="", bool recalculate = true){
  struct th_param{
    double thMinus, thPlus;
    vector<pair<double,double>> purity;
    vector<pair<double,double>> efficiency;
    vector<pair<double,double>> NCorrect_of_all_this_type;
  };
  vector <th_param> threshold_table;
  if(file==string("") && recalculate){
    cout << "Empty file string! Set recalculate to false\n";
    recalculate=false;
  }
  if(recalculate){
    TChain *t = new TChain();
    auto countFiles = t->Add((file+string("/hf_by_processID")).c_str());
    std::cout << "Count of files: " << countFiles << std::endl;
  
    vector<pair<double,double>> thresholds = {};
    for(double t = 1; t < 20; t+=0.5) thresholds.push_back(make_pair(t,t));
    // for(double tm = 7; tm < 10.5; tm+=0.5) for(double tp = 7; tp < 13.5; tp+=0.5) thresholds.push_back(make_pair(tm,tp));
    
    std::vector<string> processes = {"","&&processID==101","&&processID==102","&&processID==103","&&processID==104","&&processID==105"};
    std::vector<string> processTypes = {"All","ND ","CD ","SDm","SDp","DD "};
    std::cout << "Summary by all processID's" << std::endl;
    TString defaultCut = "l1Triggers[3]==1";
    std::vector<int> countTotal = {};
    countTotal.push_back(t->GetEntries(defaultCut+""));                 std::cout << "All:               " << countTotal.back() << std::endl;
    countTotal.push_back(t->GetEntries(defaultCut+"&&processID==101")); std::cout << "101 (ND):          " << countTotal.back() << std::endl;
    countTotal.push_back(t->GetEntries(defaultCut+"&&processID==102")); std::cout << "102 (CD):          " << countTotal.back() << std::endl;
    countTotal.push_back(t->GetEntries(defaultCut+"&&processID==103")); std::cout << "103 (SDm [RG@-Z]): " << countTotal.back() << std::endl;
    countTotal.push_back(t->GetEntries(defaultCut+"&&processID==104")); std::cout << "104 (SDp [RG@+Z]): " << countTotal.back() << std::endl;
    countTotal.push_back(t->GetEntries(defaultCut+"&&processID==105")); std::cout << "105 (DD):          " << countTotal.back() << std::endl;
    std::vector<int> countTotal_errors = {}; for (auto c: countTotal) countTotal_errors.push_back(TMath::Sqrt(c));
    vector<pair<char,char>> operations = {make_pair('>','>'),make_pair('<','<'),make_pair('<','>'),make_pair('>','<')};
    for(auto thr: thresholds){
      th_param param = {thr.first, thr.second};
      std::cout << "===================================================================================================\n";
      std::cout << "For thresholds on minus and plus side equals " << thr.first << " and " << thr.second << ":\n";
      for(int type = 1; type < 5; type++){
        std::cout << "  For events counted as " << processTypes[type] << ":\n";
        std::vector<int> countByType = {};
        for (int p = 0; p < 6; p++){
          TString cut = defaultCut;
          cut += "&&hfMinus"; cut += operations.at(type-1).first ; cut += thr.first;
          cut += "&&hfPlus";  cut += operations.at(type-1).second; cut += thr.second;
          cut += processes[p];
          countByType.push_back(t->GetEntries(cut));
          // std::cout << "// cut: " << cut << std::endl;
          std::cout << "    " << processTypes[p] << " processes [" << 100+p << "]: " << t->GetEntries(cut) << std::endl;
        }
        std::vector<int> countByType_errors = {}; for(auto c: countByType) countByType_errors.push_back(TMath::Sqrt(c));
        double purity = double(countByType[type])/double(countByType[0]);
        double purityError = GetErrorOfDivision (countByType[type],countByType_errors[type],countByType[0],countByType_errors[0]);
        param.purity.push_back(make_pair(purity, purityError));
        double efficiency = double(countByType[0])/double(countTotal[0]);
        double efficiencyError = GetErrorOfDivision(countByType[0], countByType_errors[0], countTotal[0], countTotal_errors[0]);
        param.efficiency.push_back(make_pair(efficiency, efficiencyError));
        double efficiency_L = double(countByType[type])/double(countTotal[type]);
        double efficiencyError_L = GetErrorOfDivision(countByType[type], countByType_errors[type], countTotal[type], countTotal_errors[type]);
        param.NCorrect_of_all_this_type.push_back(make_pair(efficiency_L, efficiencyError_L));
        std::cout << "    => purity       = " << purity        << " ± " << purityError       << std::endl;
        std::cout << "       efficiency   = " << efficiency    << " ± " << efficiencyError   << std::endl;
        std::cout << "       efficiency_L = " << efficiency_L  << " ± " << efficiencyError_L << std::endl;
      }
      threshold_table.push_back(param);
    }
  //   {
  //     ofstream texTable ("summary_table.tex");
  //     texTable.precision(3);
  //     texTable << "%===================================================================================================\n";
  //     texTable << "%    SUMMARY TABLE \n";
  //     texTable << "%===================================================================================================\n";
  //     texTable << "\\begin{frame} \\frametitle{}" << endl;
  //     texTable << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
  //     texTable << "  \\hline" << endl;
  //     texTable << "  \\multicolumn{2}{|c|}{Thresholds} & \\multicolumn{2}{c|}{SDm [RG@-Z]}  & \\multicolumn{2}{c|}{SDp [RG@+Z]} \\\\" << endl;
  //     texTable << "  \\hline" << endl;
  //     texTable << "  HF- & HF+ & Purity & Efficiency & Purity & Efficiency \\\\" << endl;
  //     texTable << "  \\hline" << endl;
  //     // std::cout << "\\hline" << std::endl;
  //     // std::cout << "\\hline" << std::endl;
  //     int count = 0;
  //     for(auto param: threshold_table){
  //       texTable << "  " << param.thMinus << " & " << param.thPlus  << " & " << endl;
  //       texTable << "    \\begin{tabular}{@{}c@{}}" << param.purity[2].first << " $\\pm$\\\\" << param.purity[2].second << "\\end{tabular}" << " & " << endl;
  //       texTable << "    \\begin{tabular}{@{}c@{}}" << param.efficiency[2].first << " $\\pm$\\\\" << param.efficiency[2].second << "\\end{tabular}" << " & " << endl;
  //       texTable << "    \\begin{tabular}{@{}c@{}}" << param.purity[3].first << " $\\pm$\\\\" << param.purity[3].second << "\\end{tabular}" << " & " << endl;
  //       texTable << "    \\begin{tabular}{@{}c@{}}" << param.efficiency[3].first << " $\\pm$\\\\" << param.efficiency[3].second << "\\end{tabular}" << " \\\\ " << endl;
  //       texTable << "  \\hline" << endl;
  //       if((++count)%8==0){
  //         texTable << "\\end{tabular}" << endl;
  //         texTable << "\\end{frame}\n\\begin{frame}" << endl;
  //         texTable << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
  //         texTable << "  \\hline" << endl;
  //         texTable << "  \\multicolumn{2}{|c|}{Thresholds} & \\multicolumn{2}{c|}{SDm [RG@-Z]}  & \\multicolumn{2}{c|}{SDp [RG@+Z]} \\\\" << endl;
  //         texTable << "  \\hline" << endl;
  //         texTable << "  HF- & HF+ & Purity & Efficiency & Purity & Efficiency \\\\" << endl;
  //         texTable << "  \\hline" << endl;
  //       }
        
  //     }
  //     texTable << "  \\hline" << endl;
  //     texTable << "\\end{tabular}" << endl;
  //     texTable << "\\end{frame}" << endl;    
  //     texTable.close();
  //   };
  }
  {
    auto fvpFile = TFile::Open("fake_veto_probability.root", "r"); 
    auto fvpm = (TH1F*) fvpFile->Get("FakeVetoProbabilityHF-");
    auto fvpp = (TH1F*) fvpFile->Get("FakeVetoProbabilityHF+");
    auto theffm = (TH1F*) fvpFile->Get("ThresholdEfficiencyHF-");
    auto theffp = (TH1F*) fvpFile->Get("ThresholdEfficiencyHF+");
    fvpm->SetLineColor(46); fvpp->SetLineColor(36);
    fvpm->SetMarkerColor(46); fvpp->SetMarkerColor(36);
    theffm->SetLineColor(46); theffp->SetLineColor(36);
    theffm->SetMarkerColor(46); theffp->SetMarkerColor(36);
    TH1F *purityMinus, *efficiencyMinus, *purityPlus, *efficiencyPlus,
      *purityND, *efficiencyND,
      *efficiencyFullMinus, *efficiencyFullPlus,
      *efficiencyOutMinus,* efficiencyOutPlus,
      *NCorrect_of_all_this_type_minus, *NCorrect_of_all_this_type_plus;
    TFile* outfile;
    if(recalculate){
      outfile = TFile::Open("purity-efficiency.root", "RECREATE"); 
      purityMinus = new TH1F("PurityMinus","Purity Veto--;Threshold", 190, 1, 20);
      purityPlus  = new TH1F("PurityPlus" ,"Purity Veto+;Threshold" , 190, 1, 20);
      efficiencyMinus = new TH1F("EfficiencyMinus","E Veto--;Threshold", 190, 1, 20);
      efficiencyPlus = new TH1F("EfficiencyPlus"  ,"E Veto+;Threshold", 190, 1, 20);
      purityND = new TH1F("PurityND","Purity for ND;Threshold", 190, 1, 20);
      efficiencyND = new TH1F("EfficiencyND","Efficiency for ND;Threshold", 190, 1, 20);
      NCorrect_of_all_this_type_minus = new TH1F("NCorrect_of_all_this_type_SDm","E_{SDm};Threshold", 190, 1, 20);
      NCorrect_of_all_this_type_plus = new TH1F("NCorrect_of_all_this_type_SDp", "E_{SDp};Threshold", 190, 1, 20);
      purityMinus->SetLineColor(2); efficiencyMinus->SetLineColor(6);
      purityMinus->SetMarkerColor(2); efficiencyMinus->SetMarkerColor(6);
      purityPlus->SetLineColor(4); efficiencyPlus->SetLineColor(8);
      purityPlus->SetMarkerColor(4); efficiencyPlus->SetMarkerColor(8);
      purityND->SetLineColor(3); efficiencyND->SetLineColor(9);
      purityND->SetMarkerColor(3); efficiencyND->SetMarkerColor(9);
      NCorrect_of_all_this_type_minus->SetLineColor(6); NCorrect_of_all_this_type_minus->SetMarkerColor(6);
      NCorrect_of_all_this_type_plus->SetLineColor(8); NCorrect_of_all_this_type_plus->SetMarkerColor(8);
      for(auto param: threshold_table){
        if(param.thMinus != param.thPlus) continue;
        purityMinus->Fill(param.thMinus, param.purity[2].first);
        purityMinus->SetBinError(purityMinus->FindBin(param.thMinus), param.purity[2].second);
        efficiencyMinus->Fill(param.thMinus, param.efficiency[2].first);
        efficiencyMinus->SetBinError(efficiencyMinus->FindBin(param.thMinus), param.efficiency[2].second);
        purityPlus->Fill(param.thMinus, param.purity[3].first);
        purityPlus->SetBinError(purityPlus->FindBin(param.thMinus), param.purity[3].second);
        efficiencyPlus->Fill(param.thMinus, param.efficiency[3].first);
        efficiencyPlus->SetBinError(efficiencyPlus->FindBin(param.thMinus), param.efficiency[3].second);
        purityND->Fill(param.thMinus, param.purity[0].first);
        purityND->SetBinError(purityND->FindBin(param.thMinus), param.purity[0].second);
        efficiencyND->Fill(param.thMinus, param.efficiency[0].first);
        efficiencyND->SetBinError(efficiencyND->FindBin(param.thMinus), param.efficiency[0].second);
        NCorrect_of_all_this_type_minus->Fill(param.thMinus, param.NCorrect_of_all_this_type[2].first);
        NCorrect_of_all_this_type_minus->SetBinError(NCorrect_of_all_this_type_minus->FindBin(param.thMinus), param.NCorrect_of_all_this_type[2].second);
        NCorrect_of_all_this_type_plus->Fill(param.thMinus, param.NCorrect_of_all_this_type[3].first);
        NCorrect_of_all_this_type_plus->SetBinError(NCorrect_of_all_this_type_plus->FindBin(param.thMinus), param.NCorrect_of_all_this_type[2].second);
      }
      efficiencyFullMinus = (TH1F*) efficiencyMinus->Clone("efficiencyFullMinus");
      efficiencyFullMinus->Multiply(theffm);
      efficiencyFullPlus =  (TH1F*) efficiencyPlus->Clone( "efficiencyFullPlus");
      efficiencyFullPlus->Multiply( theffp);
      efficiencyOutMinus = (TH1F*) efficiencyFullMinus->Clone("efficiencyOutMinus");
      efficiencyOutMinus->Multiply(purityMinus);
      efficiencyOutPlus = (TH1F*) efficiencyFullPlus->Clone("efficiencyOutPlus");
      efficiencyOutPlus->Multiply(purityPlus);
    }else{
      auto histfile = TFile::Open("purity-efficiency.root", "r");
      // histfile->ls();
      purityMinus = (TH1F*)histfile->Get("PurityMinus");
      efficiencyMinus = (TH1F*)histfile->Get("EfficiencyMinus");
      purityPlus = (TH1F*)histfile->Get("PurityPlus");
      efficiencyPlus = (TH1F*)histfile->Get("EfficiencyPlus");
      purityND = (TH1F*)histfile->Get("PurityND");
      efficiencyND = (TH1F*)histfile->Get("EfficiencyND");
      efficiencyFullMinus = (TH1F*)histfile->Get("efficiencyFullMinus");
      efficiencyFullPlus = (TH1F*)histfile->Get("efficiencyFullPlus");
      efficiencyOutMinus = (TH1F*)histfile->Get("efficiencyOutMinus");
      efficiencyOutPlus = (TH1F*)histfile->Get("efficiencyOutPlus");
      NCorrect_of_all_this_type_minus = (TH1F*)histfile->Get("NCorrect_of_all_this_type_SDm");
      NCorrect_of_all_this_type_plus = (TH1F*)histfile->Get("NCorrect_of_all_this_type_SDp");
    }
    purityMinus->SetTitle("Purity Veto-;Threshold");
    purityPlus->SetTitle("Purity Veto+;Threshold");
    efficiencyMinus->SetTitle("#epsilon_{Veto-};Threshold");
    efficiencyPlus->SetTitle("#epsilon_{Veto+};Threshold");
    NCorrect_of_all_this_type_minus->SetTitle("#epsilon_{SDm};Threshold");
    NCorrect_of_all_this_type_plus->SetTitle("#epsilon_{SDp};Threshold");
    efficiencyFullMinus->SetTitle("Corrected selection efficiency HF-; Threshold");
    efficiencyFullPlus->SetTitle( "Corrected selection efficiency HF+; Threshold");
    efficiencyOutMinus->SetTitle("Selection acceptance HF-; Threshold");
    efficiencyOutPlus->SetTitle("Selection acceptance HF+; Threshold");
#ifdef DEBUG
    cout << "Minus side:\ni\tbin\tPurity\t\tEfficiency\tFVP\t\t(1-FVP)\t\tEff_full\tEff_out\n";
    cout.precision(7);
    for(auto binN = 1; binN < purityMinus->GetNbinsX()+1; binN++){
      if(!purityMinus->GetBinContent(binN)) continue;
      cout << purityMinus->GetBinCenter(binN) << "\t" << binN << "\t" <<
        purityMinus->GetBinContent(binN) << "\t" <<
        efficiencyMinus->GetBinContent(binN) << "\t" <<
        fvpm->GetBinContent(binN) << "\t" <<
        theffm->GetBinContent(binN) << "\t" <<
        efficiencyFullMinus->GetBinContent(binN) << "\t" <<
        efficiencyOutMinus->GetBinContent(binN) << endl;
    }
    cout << "Plus side:\ni\tbin\tPurity\t\tEfficiency\tFVP\t\t(1-FVP)\t\tEff_full\tEff_out\n";
    cout.precision(7);
    for(auto binN = 1; binN < purityPlus->GetNbinsX()+1; binN++){
      if(!purityPlus->GetBinContent(binN)) continue;
      cout << purityPlus->GetBinCenter(binN) << "\t" << binN << "\t" <<
        purityPlus->GetBinContent(binN) << "\t" <<
        efficiencyPlus->GetBinContent(binN) << "\t" <<
        fvpp->GetBinContent(binN) << "\t" <<
        theffp->GetBinContent(binN) << "\t" <<
        efficiencyFullPlus->GetBinContent(binN) << "\t" <<
        efficiencyOutPlus->GetBinContent(binN) << endl;
    }
#endif
    auto canv = new TCanvas("c","c",800,600);
    gPad->SetMargin(0.07, 0.05, 0.09, 0.02);
    auto leg = new TLegend(0.59, 0.79, 0.99, 0.99);
    {
      auto hs = new THStack ("Purity&Efficiency",";Threshold");
      hs->Add(purityMinus); leg->AddEntry(purityMinus);
      hs->Add(efficiencyMinus); leg->AddEntry(efficiencyMinus);
      hs->Add(purityPlus); leg->AddEntry(purityPlus);
      hs->Add(efficiencyPlus); leg->AddEntry(efficiencyPlus);

      hs->Draw("nostack,e1p");
      leg->Draw();
      canv->SetLogy(1);
      canv->SaveAs("purity_vs_efficiency_log.eps");
      canv->SetLogy(0);
      canv->SaveAs("purity_vs_efficiency.eps");
    };
    leg->Clear();
    canv->Clear();
    {
      auto hs2 = new THStack ("Purity&Efficiency_FVP",";Threshold");
      hs2->Add(purityMinus); leg->AddEntry(purityMinus);
      hs2->Add(efficiencyMinus); leg->AddEntry(efficiencyMinus);
      hs2->Add(purityPlus); leg->AddEntry(purityPlus);
      hs2->Add(efficiencyPlus); leg->AddEntry(efficiencyPlus);
      hs2->Add(fvpm); leg->AddEntry(fvpm);
      hs2->Add(fvpp); leg->AddEntry(fvpp);

      hs2->Draw("nostack,e1p");
      leg->Draw();
      canv->SetLogy(1);
      canv->SaveAs("purity_vs_efficiency_fvp_log.eps");
      canv->SetLogy(0);
      canv->SaveAs("purity_vs_efficiency_fvp.eps");
    };
    leg->Clear();
    canv->Clear();
    {
      auto hs3 = new THStack ("Purity&SD_candidates",";Threshold");
      auto leg3 = new TLegend(0.69, 0.79, 0.99, 0.99);
      hs3->Add(efficiencyFullMinus); leg3->AddEntry(efficiencyFullMinus);
      hs3->Add(purityMinus); leg3->AddEntry(purityMinus);
      hs3->Add(purityPlus); leg3->AddEntry(purityPlus);
      hs3->Add(efficiencyFullPlus); leg3->AddEntry(efficiencyFullPlus);

      hs3->Draw("nostack,e1p");
      leg3->Draw();
      canv->SetLogy(1);
      canv->SaveAs("purity_vs_efficiency_corrected_log.eps");
      canv->SetLogy(0);
      canv->SaveAs("purity_vs_efficiency_corrected.eps");
      leg3->Clear();
    };
    canv->Clear();
    // {
    //   auto hs4 = new THStack ("Purity&Efficiency_out",";Threshold");
    //   auto leg4 = new TLegend(0.49, 0.14, 0.99, 0.29); //0.14, 0.14, 0.54, 0.34);
    //   hs4->Add(efficiencyOutMinus); leg4->AddEntry(efficiencyOutMinus);
    //   hs4->Add(efficiencyOutPlus); leg4->AddEntry(efficiencyOutPlus);

    //   hs4->Draw("nostack,e1p");
    //   leg4->Draw();
    //   canv->SetLogy(1);
    //   canv->SaveAs("purity_vs_out_efficiency_log.eps");
    //   canv->SetLogy(0);
    //   canv->SaveAs("purity_vs_out_efficiency.eps");
    //   leg4->Clear();
    // };
    // canv->Clear();
    {
      auto hs4 = new THStack ("Purity&Efficiency_acceptance",";Threshold");
      auto leg4 = new TLegend(0.49, 0.14, 0.99, 0.29); //0.14, 0.14, 0.54, 0.34);
      ((TGaxis*)efficiencyOutMinus->GetXaxis())->SetMaxDigits(3);
      ((TGaxis*)efficiencyOutPlus->GetXaxis())->SetMaxDigits(3);
      hs4->Add(efficiencyOutMinus); leg4->AddEntry(efficiencyOutMinus);
      hs4->Add(efficiencyOutPlus); leg4->AddEntry(efficiencyOutPlus);
      hs4->Add(purityMinus); leg4->AddEntry(purityMinus);
      hs4->Add(purityPlus); leg4->AddEntry(purityPlus);

      hs4->Draw("nostack,e1p");
      leg4->Draw();
      canv->SetLogy(1);
      canv->SaveAs("purity_vs_efficiency_acceptance_log.eps");
      canv->SetLogy(0);
      canv->SaveAs("purity_vs_efficiency_acceptance.eps");
      leg4->Clear();
    };
    canv->Clear();
    {
      auto hs4 = new THStack ("Efficiency_SD",";Threshold");
      auto leg4 = new TLegend(0.79, 0.14, 0.99, 0.29); //0.14, 0.14, 0.54, 0.34);
      hs4->Add(NCorrect_of_all_this_type_minus); leg4->AddEntry(NCorrect_of_all_this_type_minus);
      hs4->Add(NCorrect_of_all_this_type_plus); leg4->AddEntry(NCorrect_of_all_this_type_plus);

      hs4->Draw("nostack,e1p");
      leg4->Draw();
      canv->SetLogy(1);
      canv->SaveAs("efficiency_SD_log.eps");
      canv->SetLogy(0);
      canv->SaveAs("efficiency_SD.eps");
      leg4->Clear();
    };
    canv->Clear();
    {
      gPad->SetMargin(0.07, 0.05, 0.09, 0.09);
      auto hs4 = new THStack ("Purity&Efficiency_acceptance",";Threshold");
      auto leg4 = new TLegend(0.49, 0.14, 0.99, 0.29); //0.14, 0.14, 0.54, 0.34);
      ((TGaxis*)efficiencyOutMinus->GetXaxis())->SetMaxDigits(3);
      ((TGaxis*)efficiencyOutPlus->GetXaxis())->SetMaxDigits(3);
      hs4->Add(efficiencyOutMinus); leg4->AddEntry(efficiencyOutMinus);
      hs4->Add(efficiencyOutPlus); leg4->AddEntry(efficiencyOutPlus);

      hs4->Draw("nostack,e1p");
      leg4->Draw();
      canv->SetLogy(1);
      canv->SaveAs("efficiency_acceptance_log.eps");
      canv->SetLogy(0);
      canv->SaveAs("efficiency_acceptance.eps");
      leg4->Clear();
    };
    canv->Clear();
    if(recalculate){
      // outfile->ls();
      outfile->Write();
      outfile->Close();
    }
  };
}
