{
  TFile *fileE = TFile::Open("uaplot_output_histos_EPOS_Pbp.root");
  TFile *fileH = TFile::Open("uaplot_output_histos_Hijing_Pbp.root");
  TH1F * hEall  = (TH1F*) fileE->Get("MC/proton_pt2_h[0]");//hEall->Sumw2();  
  TH1F * hHall  = (TH1F*) fileH->Get("MC/proton_pt2_h[0]");//hHall->Sumw2();  
  // SD+
  TH1F * hEdiff  = (TH1F*) fileE->Get("MC/proton_pt2_h[19]");//hEdiff->Sumw2();  
  TH1F * hHdiff  = (TH1F*) fileH->Get("MC/proton_pt2_h[19]");//hHdiff->Sumw2();  
  TH1F * hERG1     = (TH1F*) fileE->Get("MC/proton_pt2_h[20]");//hERG1->Sumw2();  
  TH1F * hHRG1     = (TH1F*) fileH->Get("MC/proton_pt2_h[20]");//hHRG1->Sumw2();  
  TH1F * hERG10    = (TH1F*) fileE->Get("MC/proton_pt2_h[30]");//hERG10->Sumw2();  
  TH1F * hHRG10    = (TH1F*) fileH->Get("MC/proton_pt2_h[30]");//hHRG10->Sumw2();  
  /*
  int statE = hEall->GetEntries();
  int statH = hHall->GetEntries();
  TH1F * hE[4] = {hEall, hEdiff, hERG1, hERG10};
  TH1F * hH[4] = {hHall, hHdiff, hHRG1, hHRG10};
  */  
/*
  for(int i=0;i<4;i++){
    hE[i]->Scale(1./float(statE));
    hH[i]->Scale(1./float(statH));
    hH[i]->SetLineColor(1);
  };
  */
};
