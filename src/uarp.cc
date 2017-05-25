#include "uarp.h"

#include <time.h>
#include <string>
#include "iostream"
#include "TString.h"
#include "Rtypes.h"
#include "TMath.h"
#include "TH2.h"
#include "TBranch.h"

ClassImp(uarp)

float uarp::FriciVar()
{
  float var = 0; 
  if(track_valid[0] && track_valid[1]){
    var = -0.1;
  }else if (track_valid[0]){
    var = y[0][0]/y[1][0];
  }else if (track_valid[1]){
    var = y[0][1]/y[1][1];
  };
  return var;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uarp::uarp(TChain     * tree, 	    //!<tree of ua format
     TDirectory * dir,              //!<directory in the output root file
     const bool   cmstotem,         //!<true for merged, false for CMS only
     const bool   CASTORp,          //!<true proton to -Z (pPb)
     const short          int MC,   //!< -1, 0 - data; >0 MC
     const short unsigned int Ncuts //!< number of cuts
    ):uabase(cmstotem, MC, Ncuts, dir)
{
  RPproton  = 0;
  RPtrackNU = 0;
  RPtrackND = 0;
  RPtrackFU = 0;
  RPtrackFD = 0;
  if(tree_combined_flag){
    if(CASTORp){
      tree->SetBranchAddress("rec_prot_right.", &RPproton); // here (TOTEM confirmed)
       //tree->SetBranchAddress("rec_prot_left.",  &RPproton); // wrong
      tree->SetBranchAddress("track_rp_120.", &RPtrackNU);
      tree->SetBranchAddress("track_rp_121.", &RPtrackND);
      tree->SetBranchAddress("track_rp_124.", &RPtrackFU);
      tree->SetBranchAddress("track_rp_125.", &RPtrackFD);
    }else{
      tree->SetBranchAddress("rec_prot_left.",  &RPproton); // here
      //tree->SetBranchAddress("rec_prot_right.", &RPproton); // wrong (TOTEM confirmed)
      tree->SetBranchAddress("track_rp_20.", &RPtrackNU);
      tree->SetBranchAddress("track_rp_21.", &RPtrackND);
      tree->SetBranchAddress("track_rp_24.", &RPtrackFU);
      tree->SetBranchAddress("track_rp_25.", &RPtrackFD);
    };
  };
  create_histos();
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uarp::~uarp(){
  std::cout << "uarp::~uarp(): deleting " 
	    << h1D->size() << "+" << h2D->size() << " histos" << std::endl;
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
bool uarp::ProceedEvent(short unsigned int cut, const bool fill, const bool info){
  if(fill && (cut>n_cuts) ){
    std::cout << "uarp::ProceedEvent: required cut number is larger that possible, please define larger uaforward::n_cut!\n"; 
    return false;
  }; 

  proton_valid = false;
  proton_t     = 0;
  proton_xi    = 1;
  proton_y     = 0;
  memset(y, 0, sizeof(y));
  memset(track_valid, false, sizeof(track_valid));
  
  if(mc>0 || !tree_combined_flag) 
    return false;
  
  proton_valid = RPproton->valid;
  proton_xi    = RPproton->xi;
  proton_t     = RPproton->t;
  proton_y     = RPproton->y0;
  
  //double y[2][2]; // [F,N][U,D]
  if(RPtrackFU->valid){
    y[0][0] = RPtrackFU->y;
  };
  if(RPtrackFD->valid){
    y[0][1] = RPtrackFD->y;
  };
  if(RPtrackNU->valid){
    y[1][0] = RPtrackNU->y;
  };
  if(RPtrackND->valid){
    y[1][1] = RPtrackND->y;
  };
  
  for (short unsigned int ud = 0; ud<2; ud++){
    if( (y[0][ud]!=0) && (y[1][ud]!=0) ){
	track_valid[ud] = true;
    };
  };

  if(info)
    PrintEventInfo(true);
  if(fill)
    FillLastEvent(cut);

  return true;
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
void uarp::PrintEventInfo(bool detailed){
  std::cout << "uarp::PrintEventInfo:\n";
  if(detailed){
    std::cout << "\tvalid: " << proton_valid << "\tXi : "     << proton_xi << ";\tt :"<< proton_t << ";\ty :"<< proton_y << std::endl;
    std::cout << "\ttrack FU y:" << y[0][0] << "\ttrack FD y:" << y[0][1] << std::endl;
    std::cout << "\ttrack NU y:" << y[1][0] << "\ttrack ND y:" << y[1][1] << std::endl;
    std::cout << "\t yF/yN:" << FriciVar() << std::endl;
  }else{
    if(!proton_valid){
      std::cout << "no protons found\n";
    }else{
      std::cout << "\tXi : "     << proton_xi << ";\tt :"<< proton_t << std::endl;      
    };
  };
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uarp::FillLastEvent(const short unsigned int cut)
{
  if( cut>=n_cuts ){
    std::cout << "uarp::FillLastEvent: required cut number is larger that possible, do nothing. Please define larger uaforward::n_cut!\n"; 
    return false;
  }; 
  if(mc<=0){
    if(proton_valid){
      //double log_xi = TMath::Log10(proton_xi);
      //double log_t  = TMath::Log10(-proton_t);
      xi_valid_h[cut]->Fill(proton_xi);
      t_valid_h[cut]->Fill(proton_t);
      xi_vs_t_valid_h[cut]->Fill(TMath::Log10(-proton_t), proton_xi);
    };
    if(track_valid[0]) {
      friciU_h[cut]->Fill(y[1][0], y[0][0]/y[1][0]);
      yp_vs_yN[cut]->Fill(y[1][0], proton_y);
    };
    if(track_valid[1]) {
      friciD_h[cut]->Fill(y[1][1], y[0][1]/y[1][1]);
      yp_vs_yN[cut]->Fill(y[1][1], proton_y);
    };
  };
  return true;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
void uarp::FillRPWithMCtruth(const short unsigned int cut, const double xi, const double t){
  if( cut>=n_cuts ){
    std::cout << "uarp::FillRPWithMCtruth: required cut number is larger that possible, do nothing. Please define larger uaforward::n_cut!\n"; 
    return;
  }; 
  //double log_xi = TMath::Log10(xi);
  //double log_t  = TMath::Log10(-t);
  xi_valid_h[cut]->Fill(xi);
  t_valid_h[cut]->Fill(t);
  xi_vs_t_valid_h[cut]->Fill(TMath::Log10(-t), xi);
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
void uarp::create_histos(){
  TString title1, title2;
  
  n_each_h1D = n_cuts;
  n_each_h2D = n_cuts;
  xi_valid_h      = new TH1F * [n_each_h1D];
  t_valid_h       = new TH1F * [n_each_h1D];
  xi_vs_t_valid_h = new TH2F * [n_each_h2D];
  
  friciU_h = new TH2F * [n_each_h2D];
  friciD_h = new TH2F * [n_each_h2D];
  yp_vs_yN = new TH2F * [n_each_h2D];

  for(short unsigned int i=0; i<n_each_h1D; i++){
    title1 = "xi_valid_h["; title1+=i; title1+="]";
    title2 = title1; title2+="; log_{10}(#xi)";
    xi_valid_h[i] = new TH1F(title1.Data(), title2.Data(), 400, -2, 2);
    xi_valid_h[i]->SetDirectory(directory); 

    title1 = "t_valid_h["; title1+=i; title1+="]";
    title2 = title1; title2+="; t [(GeV/c)^{2}]";
    t_valid_h[i] = new TH1F(title1.Data(), title2.Data(), 10001, -10000, 1);
    t_valid_h[i]->SetDirectory(directory); 
  };
  h1D->push_back(xi_valid_h);
  h1D->push_back(t_valid_h);

  for(short unsigned int i=0; i<n_each_h2D; i++){
    title1 = "xi_vs_t_valid_h["; title1+=i; title1+="]";
    title2 = title1; title2+=";log_{10}(-t [(GeV/c)^{2}]); #xi";//log_{10}(#xi)";
    xi_vs_t_valid_h[i] = new TH2F(title1.Data(), title2.Data(), 700, -3, 4, 400, -2, 2);
    xi_vs_t_valid_h[i]->SetDirectory(directory); 

    title1 = "friciU_h["; title1+=i; title1+="]";
    title2 = title1; title2+=";y_N [mm]; y_F/y_N";
    friciU_h[i] = new TH2F(title1.Data(), title2.Data(), 250, 5, 30, 400, 0.75, 1.15);
    friciU_h[i]->SetDirectory(directory); 

    title1 = "friciD_h["; title1+=i; title1+="]";
    title2 = title1; title2+=";y_N [mm]; y_F/y_N";
    friciD_h[i] = new TH2F(title1.Data(), title2.Data(), 250, -30, -5, 400, 0.75, 1.15);
    friciD_h[i]->SetDirectory(directory); 

    title1 = "yp_vs_yN["; title1+=i; title1+="]";
    title2 = title1; title2+=";y_N [mm]; y0(p) [mm]";
    yp_vs_yN[i] = new TH2F(title1.Data(), title2.Data(), 600, -30, 30, 600, -30, 30);
    yp_vs_yN[i]->SetDirectory(directory); 
  };  
  h2D->push_back(xi_vs_t_valid_h);
  h2D->push_back(friciU_h);
  h2D->push_back(friciD_h);
  h2D->push_back(yp_vs_yN);
};

