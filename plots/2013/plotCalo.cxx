{
  gStyle->SetOptFit(1);
  TFile *f = TFile::Open("/home/kkuzn/Analysis/uaplotter1/uaplot_output_histos_data_pPb.root","r");
  TH1F *hNall, *hNND, *hNEL, *hNSDm, *hNSDp;
  TH1F *hSall, *hSND, *hSEL, *hSSDm, *hSSDp;
  TH1F *h_all, *h_ND, *h_EL, *h_SDm, *h_SDp;
//  TString histoName1="CMScalo/hf_max_towerE_minus_h[";
//  TString histoName2="CMScalo/hf_max_towerE_plus_h[";  
  TString histoName1="CMScalo/hf_etotal_minus_h[";
  TString histoName2="CMScalo/hf_etotal_plus_h[";  
  //TString histoName3="ZDC/FSCm_N_h[";  
  TString hN;
  
  hN = histoName1; hN+=0; hN+="]";
  hNall=(TH1F*) f->Get(hN.Data());
  hNall->SetLineColor(1);
  hNall->SetLineWidth(2);
  hN = histoName2; hN+=0; hN+="]";
  hSall=(TH1F*) f->Get(hN.Data());
  hSall->SetLineColor(1); hSall->SetMarkerColor(1); hSall->SetMarkerStyle(20);
  hSall->SetLineWidth(2);
  
    hN = histoName1; hN+=17; hN+="]";
  hNND=(TH1F*) f->Get(hN.Data());
  hNND->SetLineColor(2); 
  hNND->SetLineWidth(2);
  hN = histoName2; hN+=17; hN+="]";
  hSND=(TH1F*) f->Get(hN.Data());
  hSND->SetLineColor(2); hSND->SetMarkerColor(2); hSND->SetMarkerStyle(20);
  hSND->SetLineWidth(2);

  hN = histoName1; hN+=3; hN+="]";
  hNEL=(TH1F*) f->Get(hN.Data());
  hNEL->SetLineColor(kGreen+3); 
  hNEL->SetLineWidth(2);
  hN = histoName2; hN+=3; hN+="]";
  hSEL=(TH1F*) f->Get(hN.Data());
  hSEL->SetLineColor(kGreen+3); hSEL->SetMarkerColor(kGreen+3); hSEL->SetMarkerStyle(20);
  hSEL->SetLineWidth(2);

  
  hN = histoName1; hN+=4; hN+="]";
  hNSDm=(TH1F*) f->Get(hN.Data());
  hNSDm->SetLineColor(kBlue);
  hNSDm->SetLineWidth(2);
  hN = histoName2; hN+=4; hN+="]";
  hSSDm=(TH1F*) f->Get(hN.Data());
  hSSDm->SetLineColor(kBlue); hSSDm->SetMarkerColor(kBlue); hSSDm->SetMarkerStyle(20);
  hSSDm->SetLineWidth(2);

  hN = histoName1; hN+=19; hN+="]";
  hNSDp=(TH1F*) f->Get(hN.Data());
  hNSDp->SetLineColor(kMagenta);
  hNSDp->SetLineWidth(2);
  hN = histoName2; hN+=19; hN+="]";
  hSSDp=(TH1F*) f->Get(hN.Data());
  hSSDp->SetLineColor(kMagenta); hSSDp->SetMarkerColor(kMagenta); hSSDp->SetMarkerStyle(20);
  hSSDp->SetLineWidth(2);

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  hSall->Draw();
  hSND->Draw("sames");
  hSEL->Draw("sames");
  hSSDm->Draw("sames");
  hSSDp->Draw("sames");
};