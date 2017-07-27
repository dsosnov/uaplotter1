TH1F plot_xi(TFile *f, int cut, bool norm){
  TString hname = "RP/xi_valid_h["; hname += cut; hname += "]";
  TH1F * horig  = (TH1F*) f->Get(hname.Data());
  TH1F h; horig->Copy(h);
  if(norm){
    int stat = h.GetEntries();
    h.Sumw2();
    h.Scale(1./stat);
  };
  return h;
};

TH1F plot_t(TFile *f, int cut, bool norm){
  TString hname = "RP/t_valid_h["; hname += cut; hname += "]";
  TH1F * horig  = (TH1F*) f->Get(hname.Data());
  TH1F h; horig->Copy(h);
  if(norm){
    int stat = h.GetEntries();
    h.Sumw2();
    h.Scale(1./stat);
  };
  return h;
};


TH1F plot_Yratio_profile(TFile *f, int cut, bool norm){
  TString hname = "RP/friciD_h["; hname += cut; hname += "]";
  TH2F * horig  = (TH2F*) f->Get(hname.Data());
  hname         = "ratio_Y[";     hname += cut; hname += "]";
  TH1D *hp = horig->ProjectionY("",0,horig->GetNbinsX(),"e");
  TH1F hc; hp->Copy(hc);
  if(norm){
    int stat = hc.GetEntries();
    hc.Sumw2();
    hc.Scale(1./stat);
  };
  return hc;  
}

TH1F plot_Xi_reco(TFile *f, int cut, bool norm){
  TString hname = "xi_reco_full_h["; hname += cut; hname += "]";
  TH1F * horig  = (TH1F*) f->Get(hname.Data());
  TH1F h; horig->Copy(h);
  if(norm){
    int stat = h.GetEntries();
    h.Sumw2();
    h.Scale(1./stat);
  };
  return h;
}
