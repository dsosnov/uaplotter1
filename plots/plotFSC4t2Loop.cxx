{
  gStyle->SetOptFit(1);
  TFile *f = TFile::Open("/home/kkuzn/Analysis/uaplotter1/uaplot_output_histos_data_pPb.root","r");
  TH1F *hall[8], *hnall[8], *hsall[8];
  TH1F *hND[8],  *hnND[8],  *hsND[8];
  TH1F *hEL[8],  *hnEL[8],  *hsEL[8];
  //  TH1F *hSDm[8], *hnSDm[8], *hsSDm[8];
  //  TH1F *hSDp[8], *hnSDp[8], *hsSDp[8];
  TH1F *hNall, *hNND, *hNEL, *hNSDm, *hNSDp;
  TH1F *hSall, *hSND, *hSEL, *hSSDm, *hSSDp;
  TH1F *h_all, *h_ND, *h_EL, *h_SDm, *h_SDp;
  TString histoName1="ZDC/FSCm_No_ha1[";//"ZDC/FSCm_No8_h["; //No_ha1
  TString histoName2="ZDC/FSCm_Si_ha1[";//"ZDC/FSCm_Si8_h[";  // Si_ha1 
  //  TString histoName3="ZDC/FSCm_N_h[";  
  TString hN;
  /*
  hN = histoName1; hN+=0; hN+="]";
  hNall=(TH1F*) f->Get(hN.Data());
  hNall->SetLineColor(1);
  hNall->SetLineWidth(1);
  hN = histoName2; hN+=0; hN+="]";
  hSall=(TH1F*) f->Get(hN.Data());
  hSall->SetLineColor(1); hSall->SetMarkerColor(1); hSall->SetMarkerStyle(20);
  hSall->SetLineWidth(2);
  hN = histoName3; hN+=0; hN+="]";
  h_all=(TH1F*) f->Get(hN.Data());
  h_all->SetLineColor(1);
  h_all->SetLineWidth(2);
  */
  hN = histoName1; hN+=1; hN+="]"; std::cout << hN.Data() << std::endl;
  hNND=(TH1F*) f->Get(hN.Data());
  hNND->SetLineColor(2); 
  hNND->SetLineWidth(2);
  hN = histoName2; hN+=1; hN+="]";
  hSND=(TH1F*) f->Get(hN.Data());
  hSND->SetLineColor(2); hSND->SetMarkerColor(2); hSND->SetMarkerStyle(20);
  hSND->SetLineWidth(2);
  // hN = histoName3; hN+=1; hN+="]";
  // h_ND=(TH1F*) f->Get(hN.Data());
  // h_ND->SetLineColor(2);
  // h_ND->SetLineWidth(2);

  hN = histoName1; hN+=2; hN+="]";
  hNEL=(TH1F*) f->Get(hN.Data());
  hNEL->SetLineColor(kGreen+3); 
  hNEL->SetLineWidth(1);
  hN = histoName2; hN+=2; hN+="]";
  hSEL=(TH1F*) f->Get(hN.Data());
  hSEL->SetLineColor(kGreen+3); hSEL->SetMarkerColor(kGreen+3); hSEL->SetMarkerStyle(20);
  hSEL->SetLineWidth(2);
  // hN = histoName3; hN+=2; hN+="]";
  // h_EL=(TH1F*) f->Get(hN.Data());
  // h_EL->SetLineColor(kGreen+3);
  // h_EL->SetLineWidth(2);

  /*
  hN = histoName1; hN+=4; hN+="]";
  hNSDm=(TH1F*) f->Get(hN.Data());
  hNSDm->SetLineColor(kBlue);
  hNSDm->SetLineWidth(1);
  hN = histoName2; hN+=4; hN+="]";
  hSSDm=(TH1F*) f->Get(hN.Data());
  hSSDm->SetLineColor(kBlue); hSSDm->SetMarkerColor(kBlue); hSSDm->SetMarkerStyle(20);
  hSSDm->SetLineWidth(2);
  hN = histoName3; hN+=4; hN+="]";
  h_SDm=(TH1F*) f->Get(hN.Data());
  h_SDm->SetLineColor(kBlue);
  h_SDm->SetLineWidth(2);


  hN = histoName1; hN+=19; hN+="]";
  hNSDp=(TH1F*) f->Get(hN.Data());
  hNSDp->SetLineColor(kMagenta);
  hNSDp->SetLineWidth(1);
  hN = histoName2; hN+=19; hN+="]";
  hSSDp=(TH1F*) f->Get(hN.Data());
  hSSDp->SetLineColor(kMagenta); hSSDp->SetMarkerColor(kMagenta); hSSDp->SetMarkerStyle(20);
  hSSDp->SetLineWidth(2);
  hN = histoName3; hN+=19; hN+="]";
  h_SDp=(TH1F*) f->Get(hN.Data());
  h_SDp->SetLineColor(kMagenta);
  h_SDp->SetLineWidth(2);
  */
  /*
  TCanvas * cSN = new TCanvas("cSN", "Signal-Noise Compare", 1200, 800);
  cSN->Divide(3,2);
  cSN->cd(1); cSN_1->SetLogy();
  hNall->Draw();
  hSall->Draw("sames");
  cSN->cd(2); cSN_2->SetLogy();
  hNND->Draw();
  hSND->Draw("sames");
  cSN->cd(3); cSN_3->SetLogy();
  hNEL->Draw();
  hSEL->Draw("sames");
  cSN->cd(4); cSN_4->SetLogy();
  hNSDm->Draw();
  hSSDm->Draw("sames");
  cSN->cd(5); cSN_5->SetLogy();
  hNSDp->Draw();
  hSSDp->Draw("sames");
  
  cSN->cd(6);
  h_all->Draw(); h_ND->Draw("sames"); h_EL->Draw("sames"); h_SDm->Draw("sames"); h_SDp->Draw("sames");
  
  
  for (short unsigned int ch=0; ch<8; ch++){
    TString histoName="ZDC/FSCm_TS_ha"; histoName+=ch; histoName+="[";
    
    TString hN       =histoName;       hN+=0; hN+="]";
    std::cout << hN.Data() << std::endl;
    hall[ch]  = (TH1F*) f->Get(hN.Data());
    hall[ch]->SetLineColor(1);
    hall[ch]->SetLineWidth(2);
    hN.Replace(9,2, "No");
    hnall[ch]  = (TH1F*) f->Get(hN.Data());
    hnall[ch]->SetLineColor(1);
    hnall[ch]->SetLineWidth(2);
    hN.Replace(9,2, "Si");
    hsall[ch]  = (TH1F*) f->Get(hN.Data());
    hsall[ch]->SetLineColor(1);
    hsall[ch]->SetLineWidth(2);
    
    hN       =histoName;       hN+=17; hN+="]";
    TH1F * hND[ch]  = (TH1F*) f->Get(hN.Data());
    hND[ch]->SetLineColor(2);
    hND[ch]->SetLineWidth(2);
    hN.Replace(9,2, "No");
    TH1F * hnND[ch]  = (TH1F*) f->Get(hN.Data());
    hnND[ch]->SetLineColor(2);
    hnND[ch]->SetLineWidth(2);
    hN.Replace(9,2, "Si");
    TH1F * hsND[ch]  = (TH1F*) f->Get(hN.Data());
    hsND[ch]->SetLineColor(2);
    hsND[ch]->SetLineWidth(2);
    
    hN       =histoName;       hN+=3; hN+="]";
    TH1F * hEL[ch]  = (TH1F*) f->Get(hN.Data());
    hEL[ch]->SetLineColor(kGreen+3);
    hEL[ch]->SetLineWidth(2);
    hN.Replace(9,2, "No");
    TH1F * hnEL[ch]  = (TH1F*) f->Get(hN.Data());
    hnEL[ch]->SetLineColor(kGreen+3);
    hnEL[ch]->SetLineWidth(2);
    hN.Replace(9,2, "Si");
    TH1F * hsEL[ch]  = (TH1F*) f->Get(hN.Data());
    hsEL[ch]->SetLineColor(kGreen+3);
    hsEL[ch]->SetLineWidth(2);
    /*
    hN       =histoName;       hN+=4; hN+="]";
    TH1F * hSDm[ch]  = (TH1F*) f->Get(hN.Data());
    hSDm[ch]->SetLineColor(kBlue);
    hSDm[ch]->SetLineWidth(2);
    hN.Replace(9,2, "No");
    TH1F * hnSDm[ch]  = (TH1F*) f->Get(hN.Data());
    hnSDm[ch]->SetLineColor(kBlue);
    hnSDm[ch]->SetLineWidth(2);
    hN.Replace(9,2, "Si");
    TH1F * hsSDm[ch]  = (TH1F*) f->Get(hN.Data());
    hsSDm[ch]->SetLineColor(kBlue);
    hsSDm[ch]->SetLineWidth(2);

    hN       =histoName;       hN+=19; hN+="]";
    TH1F * hSDp[ch]  = (TH1F*) f->Get(hN.Data());
    hSDp[ch]->SetLineColor(kMagenta);
    hSDp[ch]->SetLineWidth(2);
    hN.Replace(9,2, "No");
    TH1F * hnSDp[ch]  = (TH1F*) f->Get(hN.Data());
    hnSDp[ch]->SetLineColor(kMagenta);
    hnSDp[ch]->SetLineWidth(2);
    hN.Replace(9,2, "Si");
    TH1F * hsSDp[ch]  = (TH1F*) f->Get(hN.Data());
    hsSDp[ch]->SetLineColor(kMagenta);    
    hsSDp[ch]->SetLineWidth(2);
  */
  };
  
  /*
  TCanvas * c11 = new TCanvas("c11", "TS -Z: [0,3] ( chID [26-29] )", 800, 800);
  c11->Divide(2,2);
  TCanvas * c12 = new TCanvas("c12", "TS -Z: [4,7] ( chID [30-33] )", 800, 800);
  c12->Divide(2,2);
  TCanvas *c1[2] = {c11, c12};
  
  for (short unsigned int cv=0; cv<2; cv++){
    for (short unsigned int pd=0; pd<4; pd++){
      ch=cv*4+pd;
      c1[cv]->cd(pd+1);
      hSDp[ch]->Draw();
      hall[ch]->Draw("sames");
      hND[ch]->Draw("sames");
      hEL[ch]->Draw("sames");
      hSDm[ch]->Draw("sames");
      hSDp[ch]->Draw("sames");
    };
  };
  
  TCanvas * c21 = new TCanvas("c21", "Noise -Z: [0,3] ( chID [26-29] )", 800, 800);
  c21->Divide(2,2);
  TCanvas * c22 = new TCanvas("c22", "Noise -Z: [4,7] ( chID [30-33] )", 800, 800);
  c22->Divide(2,2);
  TCanvas *c2[2] = {c21, c22};
  
  for (short unsigned int cv=0; cv<2; cv++){
    for (short unsigned int pd=0; pd<4; pd++){
      ch=cv*4+pd;
      c2[cv]->cd(pd+1);
      hnall[ch]->Draw("");
      hnSDp[ch]->Draw("sames");
      hnND[ch]->Draw("sames");
      hnEL[ch]->Draw("sames");
      hnSDm[ch]->Draw("sames");
      hnSDp[ch]->Draw("sames");
    };
  };

  TCanvas * cN = new TCanvas("cN", "Noise ALL", 1200, 800);
  cN->Divide(3,2);
  TF1  * fg[6];
  float mean[6];
  float sigma[6];
  float thr[6];
  for (short unsigned int cv=0; cv<6; cv++){
    cN->cd(cv+1);
    hnall[cv]->Draw();
    hnall[cv]->GetXaxis()->SetRangeUser(1000, 3000);
    fg[cv] = new TF1("Gaussian","gaus",1500,2500); 
    std::cout << "=============== ch " << ch << std::endl;
    hnall[cv]->Fit(fg[cv],"RS");    
    mean[cv]  = fg[cv]->GetParameter(1);
    sigma[cv] = fg[cv]->GetParameter(2);
    thr[cv]   = mean[cv]+3*sigma[cv];
  };
  for (short unsigned int cv=0; cv<6; cv++){
    std::cout << ch << "\t" << mean[cv] << "\t"<< sigma[cv] << "\t" << thr[cv] << std::endl;
    
  };
  for (short unsigned int cv=0; cv<6; cv++){
    std::cout << thr[cv] << ", ";
  };std::cout << std::endl;
  */
  TCanvas * c31 = new TCanvas("c31", "signal -Z: [0,3] ( chID [26-29] )", 800, 800);
  c31->Divide(2,2);
  TCanvas * c32 = new TCanvas("c32", "signal -Z: [4,7] ( chID [30-33] )", 800, 800);
  c32->Divide(2,2);
  TCanvas *c3[2] = {c31, c32};
  
  for (short unsigned int cv=0; cv<2; cv++){
    for (short unsigned int pd=0; pd<4; pd++){
      ch=cv*4+pd;
      c3[cv]->cd(pd+1);
      // hsall[ch]->Draw();
      //      hsSDp[ch]->Draw("sames");
      hsND[ch]->Draw();
      hsEL[ch]->Draw("sames");
      //      hsSDm[ch]->Draw("sames");
      //      hsSDp[ch]->Draw("sames");
    };
  };
}
