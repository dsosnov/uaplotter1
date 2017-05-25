#include "uat2.h"

#include "T2Event.h"

#include <time.h>
#include <bitset>
#include <string>
#include "iostream"

#include "Rtypes.h"
#include "TChain.h"

ClassImp(uat2)

void uat2::delete_histo(TH1F **h, int N){
  if(h){
    for(int i = 0; i<N; i++) delete h[i];
    delete [] h;
  }
}

void uat2::delete_histo(TH1F ***h, int N, int M){
  if(h){
    std::cout << "deleting T2 2xhistos...\n";
    for(int i = 0; i<N; i++){
      if(h[i]){
	for(int j = 0; j<M; j++) 
	  delete h[i][j];
      }; delete [] h[i];
    };
    delete [] h;
  };
  std::cout << "deleted T2 2xhistos\n";
};

void uat2::delete_histo(TH2F ***h, int N, int M){
  if(h){
    std::cout << "deleting T2 2x2Dhistos...\n";
    for(int i = 0; i<N; i++){
      if(h[i]){
	for(int j = 0; j<M; j++) 
	  delete h[i][j];
      }; delete [] h[i];
    };
    delete [] h;
  };
  std::cout << "deleted T2 2x2Dhistos\n";
};




void uat2::resethistos(unsigned int cut){
  if(cut>n_cuts){
      std::cout << "uat2::resethistos: n_cut<cut!\n"; 
      return;
  };
  for(unsigned int i=0; i<n_t2_h; i++){
    n_tracks_plus_h[i]->Clear(); // cms vertices for different conditions
    n_tracks_minus_h[i]->Clear(); // cms vertices for different conditions
    n_tracks_prim_plus_h[i]->Clear(); // cms vertices for different conditions
    n_tracks_prim_minus_h[i]->Clear(); // cms vertices for different conditions
    z_impact_plus_h[i]->Clear();
    z_impact_minus_h[i]->Clear();
    z_impact_prim_plus_h[i]->Clear();
    z_impact_prim_minus_h[i]->Clear();
    for(unsigned short int j=0;j<NTCASES;j++){
      n_clusters_minus_h[i][j]->Clear();
      n_clusters_plus_h[i][j]->Clear();
      n_clusters_minusF_h[i][j]->Clear();
      n_clusters_plusN_h[i][j]->Clear();
      n_clusters_minusF_h[i][j]->Clear();
      n_clusters_plusN_h[i][j]->Clear();
      n_cluster_2h[i][j]->Clear();
      n_qtracks_2h[i][j]->Clear();
      trigger_h[i][j]->Clear();//TBD to be moved outside
    };
  };
    
};



uat2::uat2(TChain* tree, TDirectory * dir, const bool combined, const int ncuts):n_cuts(ncuts){
  directory = dir;
  T2tracks = 0;
  TOTtrig  = 0;// TBD to be moved outside
  T2primaryPlus   = 0;
  T2primaryMinus   = 0;
  if(combined){
    tree->SetBranchAddress("branchT2EV.",   &T2tracks);
    tree->SetBranchAddress("trigger_data.", &TOTtrig);// TBD to be moved outside
  };

  T2primaryMinus = new T2Event();
  T2primaryPlus  = new T2Event();
  
  // vertices
  n_t2_h          = n_cuts;
  n_tracks_plus_h  = new TH1F * [n_t2_h];
  n_tracks_minus_h = new TH1F * [n_t2_h];
  n_tracks_prim_minus_h = new TH1F * [n_t2_h];
  n_tracks_prim_plus_h  = new TH1F * [n_t2_h];
  z_impact_minus_h    = new TH1F * [n_t2_h];
  z_impact_plus_h     = new TH1F * [n_t2_h];
  z_impact_prim_minus_h    = new TH1F * [n_t2_h];
  z_impact_prim_plus_h     = new TH1F * [n_t2_h];
  
  n_clusters_minus_h  = new TH1F ** [n_t2_h];
  n_clusters_plus_h   = new TH1F ** [n_t2_h];
  n_clusters_minusF_h = new TH1F ** [n_t2_h];
  n_clusters_plusF_h  = new TH1F ** [n_t2_h];
  n_clusters_minusN_h = new TH1F ** [n_t2_h];
  n_clusters_plusN_h  = new TH1F ** [n_t2_h];
  n_cluster_2h        = new TH2F ** [n_t2_h];
  n_qtracks_2h        = new TH2F ** [n_t2_h];
  
  trigger_h          = new TH1F ** [n_t2_h];// TBD to be moved outside, other TH1 type
  
  for(unsigned int i = 0; i<n_cuts; i++){
    TString title1 = "n_tracks_plus_h["; title1+=i; title1+="]";
    TString title2 = "n_tracks_plus_h["; title2+=i; title2+="]; t2 tracks +";
    n_tracks_plus_h[i] = new TH1F(title1.Data(), title2.Data(), 51,  -1, 50);
    n_tracks_plus_h[i]->SetDirectory(directory);
    
    title1 = "n_tracks_minus_h["; title1+=i; title1+="]";
    title2 = "n_tracks_minus_h["; title2+=i; title2+="]; t2 tracks -";
    n_tracks_minus_h[i] = new TH1F(title1.Data(), title2.Data(), 51,  -1, 50);
    n_tracks_minus_h[i]->SetDirectory(directory);
    
    title1 = "n_tracks_prim_plus_h["; title1+=i; title1+="]";
    title2 = "n_tracks_prim_plus_h["; title2+=i; title2+="]; t2 primary tracks +";
    n_tracks_prim_plus_h[i] = new TH1F(title1.Data(), title2.Data(), 51,  -1, 50);
    n_tracks_prim_plus_h[i]->SetDirectory(directory);

    title1 = "n_tracks_prim_minus_h["; title1+=i; title1+="]";
    title2 = "n_tracks_prim_minus_h["; title2+=i; title2+="]; t2 primary tracks -";
    n_tracks_prim_minus_h[i] = new TH1F(title1.Data(), title2.Data(), 51,  -1, 50);
    n_tracks_prim_minus_h[i]->SetDirectory(directory);

    title1 = "z_impact_plus_h"; title1+=i; title1+="]";
    title2 = "z_impact_plus_h["; title2+=i; title2+="]; Z imp +";
    z_impact_plus_h[i] = new TH1F(title1.Data(), title2.Data(), 5000,  -500000, 500000);
    z_impact_plus_h[i]->SetDirectory(directory);

    title1 = "z_impact_minus_h["; title1+=i; title1+="]";
    title2 = "z_impact_plus_h["; title2+=i; title2+="]; Z imp -";
    z_impact_minus_h[i] = new TH1F(title1.Data(), title2.Data(), 5000,  -500000, 500000);
    z_impact_minus_h[i]->SetDirectory(directory);

    title1 = "z_impact_prim_plus_h"; title1+=i; title1+="]";
    title2 = "z_impact_prim_plus_h["; title2+=i; title2+="]; Z imp prim +";
    z_impact_prim_plus_h[i] = new TH1F(title1.Data(), title2.Data(), 5000,  -500000, 500000);
    z_impact_prim_plus_h[i]->SetDirectory(directory);

    title1 = "z_impact_prim_minus_h["; title1+=i; title1+="]";
    title2 = "z_impact_prim_plus_h["; title2+=i; title2+="]; Z imp prim -";
    z_impact_prim_minus_h[i] = new TH1F(title1.Data(), title2.Data(), 5000,  -500000, 500000);
    z_impact_prim_minus_h[i]->SetDirectory(directory);
    
    n_clusters_minus_h[i] =  new TH1F * [NTCASES];
    n_clusters_plus_h[i]  =  new TH1F * [NTCASES];
    n_clusters_minusF_h[i] =  new TH1F * [NTCASES];
    n_clusters_plusF_h[i]  =  new TH1F * [NTCASES];
    n_clusters_minusN_h[i] =  new TH1F * [NTCASES];
    n_clusters_plusN_h[i]  =  new TH1F * [NTCASES];
    n_cluster_2h[i]        =  new TH2F * [NTCASES];
    n_qtracks_2h[i]        =  new TH2F * [NTCASES];
    std::string qlabel[6] = { "NM", "FM", "NP", "FP", "1NM10FM", "1NP10FP"};
    
    trigger_h[i] = new TH1F * [NTCASES];
    //https://twiki.cern.ch/twiki/bin/view/TOTEM/PArun
    std::string label[16] = { "RP220_V",
			"RP220_H",
			"RP220_Cross",
			"RP147_V",
			"CMS & L1SA",
			"T2_Single_Arm",
			"T2",
			"T2_HighMultiplicity",
			"T1",
			"BX",
			"T2_LM",
			"L1SA",
			"RP_220_H & CMS",
			"T2 & CMS",
			"T1 & T2 & CMS",
			"CMS"};   
    
    for (unsigned short int j=0; j<NTCASES; j++){
      title1 = "n_clusters_minus_h["; title1+=i; title1+="]["; title1+=j; title1+="]";
      title2 = "n_clusters_minus_h["; title2+=i; title2+="]["; title2+=j; title2+="]; t2 pad clusters -";
      n_clusters_minus_h[i][j] = new TH1F(title1.Data(), title2.Data(), 51,  -1, 50);
      n_clusters_minus_h[i][j]->SetDirectory(directory);
      
      title1 = "n_clusters_plus_h["; title1+=i; title1+="]["; title1+=j; title1+="]";
      title2 = "n_clusters_plus_h["; title2+=i; title2+="]["; title2+=j; title2+="]; t2 pad clusters +";
      n_clusters_plus_h[i][j] = new TH1F(title1.Data(), title2.Data(), 51,  -1, 50);
      n_clusters_plus_h[i][j]->SetDirectory(directory);

      title1 = "n_clusters_minusF_h["; title1+=i; title1+="]["; title1+=j; title1+="]";
      title2 = "n_clusters_minusF_h["; title2+=i; title2+="]["; title2+=j; title2+="]; t2 pad clusters -F";
      n_clusters_minusF_h[i][j] = new TH1F(title1.Data(), title2.Data(), 51,  -1, 50);
      n_clusters_minusF_h[i][j]->SetDirectory(directory);
      
      title1 = "n_clusters_plusF_h["; title1+=i; title1+="]["; title1+=j; title1+="]";
      title2 = "n_clusters_plusF_h["; title2+=i; title2+="]["; title2+=j; title2+="]; t2 pad clusters +F";
      n_clusters_plusF_h[i][j] = new TH1F(title1.Data(), title2.Data(), 51,  -1, 50);
      n_clusters_plusF_h[i][j]->SetDirectory(directory);

      title1 = "n_clusters_minusN_h["; title1+=i; title1+="]["; title1+=j; title1+="]";
      title2 = "n_clusters_minusN_h["; title2+=i; title2+="]["; title2+=j; title2+="]; t2 pad clusters -N";
      n_clusters_minusN_h[i][j] = new TH1F(title1.Data(), title2.Data(), 51,  -1, 50);
      n_clusters_minusN_h[i][j]->SetDirectory(directory);
      
      title1 = "n_clusters_plusN_h["; title1+=i; title1+="]["; title1+=j; title1+="]";
      title2 = "n_clusters_plusN_h["; title2+=i; title2+="]["; title2+=j; title2+="]; t2 pad clusters +N";
      n_clusters_plusN_h[i][j] = new TH1F(title1.Data(), title2.Data(), 51,  -1, 50);
      n_clusters_plusN_h[i][j]->SetDirectory(directory);

      title1 = "n_cluster_2h["; title1+=i; title1+="]["; title1+=j; title1+="]";
      title2 = "n_cluster_2h["; title2+=i; title2+="]["; title2+=j; title2+="];t2 quarter; N pad clusters";
      n_cluster_2h[i][j] = new TH2F(title1.Data(), title2.Data(), 6, 0,6, 52,  -2, 50);
      for(unsigned int bin=1; bin<7; bin++) n_cluster_2h[i][j]->GetXaxis()->SetBinLabel(bin,qlabel[bin-1].data());
      n_cluster_2h[i][j]->SetDirectory(directory);

      title1 = "n_qtracks_2h["; title1+=i; title1+="]["; title1+=j; title1+="]";
      title2 = "n_qtracks_2h["; title2+=i; title2+="]["; title2+=j; title2+="];t2 quarter; N tracks";
      n_qtracks_2h[i][j] = new TH2F(title1.Data(), title2.Data(), 6, 0,6, 52,  -2, 50);
      for(unsigned int bin=1; bin<7; bin++) n_qtracks_2h[i][j]->GetXaxis()->SetBinLabel(bin,qlabel[bin-1].data());
      n_qtracks_2h[i][j]->SetDirectory(directory);
    
      // TBD to be moved outside---------->
      title1 = "trigger_h["; title1+=i; title1+="]["; title1+=j; title1+="]";
      title2 = "trigger_h["; title2+=i; title2+="]["; title2+=j; title2+="]; input_status_bits";
      trigger_h[i][j] = new TH1F(title1.Data(), title2.Data(), 17,  0, 17);
      for(unsigned int bin=1; bin<17; bin++) trigger_h[i][j]->GetXaxis()->SetBinLabel(bin,label[bin-1].data());
      trigger_h[i][j]->SetDirectory(directory);    
    };
    //<---------------
  };
};


uat2::~uat2(){
  delete_histo(n_tracks_minus_h, n_t2_h);
  delete_histo(n_tracks_plus_h,  n_t2_h);
  delete_histo(n_tracks_prim_minus_h, n_t2_h);
  delete_histo(n_tracks_prim_plus_h,  n_t2_h);
  delete_histo(z_impact_minus_h,      n_t2_h);
  delete_histo(z_impact_plus_h,       n_t2_h);
  delete_histo(n_clusters_minus_h,    n_t2_h,NTCASES);
  delete_histo(n_clusters_plus_h,     n_t2_h,NTCASES);
  delete_histo(n_clusters_minusF_h,    n_t2_h,NTCASES);
  delete_histo(n_clusters_plusF_h,     n_t2_h,NTCASES);
  delete_histo(n_clusters_minusN_h,    n_t2_h,NTCASES);
  delete_histo(n_clusters_plusN_h,     n_t2_h,NTCASES);
  delete_histo(n_cluster_2h,           n_t2_h,NTCASES);
  delete_histo(n_qtracks_2h,           n_t2_h,NTCASES);
  delete_histo(trigger_h,             n_t2_h, NTCASES); //TBD to be moved outside
  if(T2primaryPlus) delete T2primaryPlus;
  if(T2primaryMinus) delete T2primaryMinus;
}



/*
double getZimpLH0( double eta );
double getZimpRH0( double eta );
double getZimpLH1( double eta );
double getZimpRH1( double eta );
double getZimpLH2( double eta );
double getZimpRH2( double eta );
double getZimpLH3( double eta );
double getZimpRH3( double eta );
*/

int uat2::ProceedEvent(unsigned int cut, const bool fill, const bool info){
  if((cut>=n_cuts) && fill){
    std::cout << "uat2::ProceedEvent - too large cat number\n";
  };
  
  // TBD to be moved outside---------->
  //std::cout << TOTtrig->input_status_bits << "\t" << (std::bitset<16>)TOTtrig->input_status_bits << std::endl;
  for (short unsigned int b=0; b<16; b++){
    //std::cout << b << "\t" <<  (std::bitset<16>) (1<<b) << "   " << bool (TOTtrig->input_status_bits & (1<<b)) << std::endl;
    trigger[b] = bool(TOTtrig->input_status_bits & (1<<b));
  };
  T2oneArm = trigger[5];
  T2trig   = trigger[6];
  //<----------------------------------
  
  int NofT2M = 0;
  int NofT2P = 0;
  
  resetT2primary(true); // Plus
  resetT2primary(false); // Minus
  
  memset(n_clust_minus,  0, sizeof(n_clust_minus));
  memset(n_clust_plus,   0, sizeof(n_clust_plus));
  memset(n_qtracks_minus,0, sizeof(n_qtracks_minus));
  memset(n_qtracks_plus, 0, sizeof(n_qtracks_plus));
  //unsigned int NumPadCluH0; //Num pad cluster in the whole PN
  //unsigned int NumPadCluH1; //Num pad cluster in the whole PF
  //unsigned int NumPadCluH2; //Num pad cluster in the whole MN
  //unsigned int NumPadCluH3; //Num pad cluster in the whole MF
  unsigned int Ntrks = (T2tracks->TrkEntryX).size();  
  n_clust_minus[0] = T2tracks->NumPadCluH2+T2tracks->NumPadCluH3;
  n_clust_plus[0]  = T2tracks->NumPadCluH0+T2tracks->NumPadCluH1;
  n_clust_minus[1] = T2tracks->NumPadCluH2; // MN
  n_clust_plus[1]  = T2tracks->NumPadCluH0; // PN
  n_clust_minus[2] = T2tracks->NumPadCluH3; // MF
  n_clust_plus[2]  = T2tracks->NumPadCluH1; // PF

  for(unsigned int i_t2=0; i_t2<Ntrks; i_t2++){
    double Z0impact = T2tracks->TrkZImpact.at(i_t2);    
    
    double phi = T2tracks->TrkPhi_RZFit.at(i_t2);
    short unsigned int sideIndx = 0;
    if( (phi<90) || (phi>270) ){ // Near Side
      sideIndx = 1;
    }else{
      sideIndx = 2;
    };
    
    if(T2tracks->TrkEntryZ.at(i_t2)<0){
      NofT2M++;       
      n_qtracks_minus[sideIndx]++;
      if(fill)
	z_impact_minus_h[cut]->Fill(Z0impact);       
    }else{
      
      NofT2P++;
      n_qtracks_plus[sideIndx]++;
      if(fill) 
	z_impact_plus_h[cut]->Fill(Z0impact);
    };

    if( Z0impact<5000 && Z0impact>-5000 && (T2tracks->TrkChiProb.at(i_t2)>0.01)) {// && eta_t2<-5.3 && eta_t2>-6.5){
      addT2primaryTrack(i_t2, cut, fill, info);
    };  
  }; // end track loop
  
  
  n_prim_tracks_minus = T2primaryMinus->TrkEta2.size();
  n_prim_tracks_plus  = T2primaryPlus->TrkEta2.size();
  n_tracks_minus      = NofT2M;
  n_tracks_plus       = NofT2P;
  if(info)
    std::cout << "T2 : -:" << n_tracks_minus << " (" << n_prim_tracks_minus << ")\t\t+:" << n_tracks_plus << " (" << n_prim_tracks_plus << ")\n";
  if(fill){
    n_tracks_minus_h[cut]->Fill(n_tracks_minus);
    n_tracks_plus_h[cut]->Fill(n_tracks_plus);
    n_tracks_prim_minus_h[cut]->Fill(n_prim_tracks_minus);
    n_tracks_prim_plus_h[cut]->Fill(n_prim_tracks_plus);
  };
  short int j[NTCASES];
  memset(j,-1, sizeof(j));
  j[0]=0;
  if(not(T2trig)) j[1]=1;
  if(T2oneArm)    j[2]=2;
  if(T2oneArm && T2trig)      j[3]=3;
  if(not(T2oneArm)) j[4]=4;
  if(fill){
    for(short unsigned int theCase=0; theCase<NTCASES; theCase++){
      if(j[theCase]>=0){
	n_clusters_minus_h[cut][j[theCase]]->Fill(n_clust_minus[0]);
	n_clusters_plus_h[cut][j[theCase]]->Fill(n_clust_plus[0]);
	n_clusters_minusN_h[cut][j[theCase]]->Fill(n_clust_minus[1]);
	n_clusters_plusN_h[cut][j[theCase]]->Fill(n_clust_plus[1]);  
	n_clusters_minusF_h[cut][j[theCase]]->Fill(n_clust_minus[2]);
	n_clusters_plusF_h[cut][j[theCase]]->Fill(n_clust_plus[2]);
	for(short unsigned side=0; side<2; side++){ // minus, plus
	  // 0 - NM
	  // 1 - FM
	  // 2 - NP
	  // 3 - FP
	  n_cluster_2h[cut][j[theCase]]->Fill(side,   n_clust_minus[side+1]); //  0: 0 - NM, 1: 1 - FM
	  n_cluster_2h[cut][j[theCase]]->Fill(side,   -2); //  stat
	  n_cluster_2h[cut][j[theCase]]->Fill(side+2, n_clust_plus[side+1]); //  0: 2 - NP, 3: 3 - FP
	  n_cluster_2h[cut][j[theCase]]->Fill(side+2, -2); //  stat
	  
	  n_qtracks_2h[cut][j[theCase]]->Fill(side,   n_qtracks_minus[side+1]); //  0: 0 - NM, 1: 1 - FM
	  n_qtracks_2h[cut][j[theCase]]->Fill(side,   -2); //  stat
	  n_qtracks_2h[cut][j[theCase]]->Fill(side+2, n_qtracks_plus[side+1]); //  0: 2 - NP, 3: 3 - FP
	  n_qtracks_2h[cut][j[theCase]]->Fill(side+2, -2); //  stat
	};
	
	// [1] -N ; [2] -F
	int combination[2] = {-1,-1};
	if(n_clust_minus[1]>0 && n_clust_minus[2]==0) combination[0] = 1; // NM only
	if(n_clust_minus[2]>0 && n_clust_minus[1]==0) combination[0] = 10; // FM only
	if(n_clust_minus[2]==0 && n_clust_minus[1]==0) combination[0] = 15; // nothing
	if(n_clust_plus[1]>0 && n_clust_plus[2]==0) combination[1] = 1; // NM only
	if(n_clust_plus[2]>0 && n_clust_plus[1]==0) combination[1] = 10; // FM only
	if(n_clust_plus[2]==0 && n_clust_plus[1]==0) combination[1] = 15; // nothing
	
	n_cluster_2h[cut][j[theCase]]->Fill(4,combination[0]); // minus
	n_cluster_2h[cut][j[theCase]]->Fill(5,combination[1]);   // plus
	
	combination[0] = -1; combination[1] = -1;
	if(n_qtracks_minus[1]>0 && n_qtracks_minus[2]==0) combination[0] = 1; // NM only
	if(n_qtracks_minus[2]>0 && n_qtracks_minus[1]==0) combination[0] = 10; // FM only
	if(n_qtracks_minus[2]==0 && n_qtracks_minus[1]==0) combination[0] = 15; // nothing
	if(n_qtracks_plus[1]>0 && n_qtracks_plus[2]==0) combination[1] = 1; // NM only
	if(n_qtracks_plus[2]>0 && n_qtracks_plus[1]==0) combination[1] = 10; // FM only
	if(n_qtracks_plus[2]==0 && n_qtracks_plus[1]==0) combination[1] = 15; // nothing
	
	n_qtracks_2h[cut][j[theCase]]->Fill(4,combination[0]); // minus
	n_qtracks_2h[cut][j[theCase]]->Fill(5,combination[1]);   // plus
	for(short unsigned int b=0;b<16;b++)
	  if(trigger[b]) trigger_h[cut][j[theCase]]->Fill(b); 
      }; // end if (j[theCase]>=0)
    }; // end trigger case loop
  }; // end if(fill)
  
  return n_prim_tracks_minus + n_prim_tracks_plus;
};

/*
// these cuts are defined for 8 TeV!!!
double getZimpLH0( double eta ) {
  double ZimpL;
  if ( fabs(eta)<6.40 && fabs(eta)>6.35 ) ZimpL=-7200.;
  if ( fabs(eta)<6.35 && fabs(eta)>6.30 ) ZimpL=-8500.;
  if ( fabs(eta)<6.30 && fabs(eta)>6.25 ) ZimpL=-6791.97;
  if ( fabs(eta)<6.25 && fabs(eta)>6.20 ) ZimpL=-7300.;
  if ( fabs(eta)<6.20 && fabs(eta)>6.15 ) ZimpL=-7100.13;
  if ( fabs(eta)<6.15 && fabs(eta)>6.10 ) ZimpL=-6287.03;
  if ( fabs(eta)<6.10 && fabs(eta)>6.05 ) ZimpL=-5240.97;
  if ( fabs(eta)<6.05 && fabs(eta)>6.00 ) ZimpL=-5500.;
  if ( fabs(eta)<6.00 && fabs(eta)>5.95 ) ZimpL=-4900.;
  if ( fabs(eta)<5.95 && fabs(eta)>5.90 ) ZimpL=-4200.;
  if ( fabs(eta)<5.90 && fabs(eta)>5.85 ) ZimpL=-4100.;
  if ( fabs(eta)<5.85 && fabs(eta)>5.80 ) ZimpL=-3900.;
  if ( fabs(eta)<5.80 && fabs(eta)>5.75 ) ZimpL=-3800.;
  if ( fabs(eta)<5.75 && fabs(eta)>5.70 ) ZimpL=-4100.;
  if ( fabs(eta)<5.70 && fabs(eta)>5.65 ) ZimpL=-4900.;
  if ( fabs(eta)<5.65 && fabs(eta)>5.60 ) ZimpL=-4700.;
  if ( fabs(eta)<5.45 && fabs(eta)>5.40 ) ZimpL=-4200.;
  if ( fabs(eta)<5.40 && fabs(eta)>5.35 ) ZimpL=-4000.;
  return (ZimpL);
}
double getZimpRH0( double eta ) {
  double ZimpR;
  if ( fabs(eta)<6.40 && fabs(eta)>6.35 ) ZimpR=8864.64;
  if ( fabs(eta)<6.35 && fabs(eta)>6.30 ) ZimpR=9430.82;
  if ( fabs(eta)<6.30 && fabs(eta)>6.25 ) ZimpR=6600.;
  if ( fabs(eta)<6.25 && fabs(eta)>6.20 ) ZimpR=7500.;
  if ( fabs(eta)<6.20 && fabs(eta)>6.15 ) ZimpR=6600.;
  if ( fabs(eta)<6.15 && fabs(eta)>6.10 ) ZimpR=5800.;
  if ( fabs(eta)<6.10 && fabs(eta)>6.05 ) ZimpR=4900.;
  if ( fabs(eta)<6.05 && fabs(eta)>6.00 ) ZimpR=5527.29;
  if ( fabs(eta)<6.00 && fabs(eta)>5.95 ) ZimpR=5463.97;
  if ( fabs(eta)<5.95 && fabs(eta)>5.90 ) ZimpR=4836.54;
  if ( fabs(eta)<5.90 && fabs(eta)>5.85 ) ZimpR=4895.51;
  if ( fabs(eta)<5.85 && fabs(eta)>5.80 ) ZimpR=4658.03;
  if ( fabs(eta)<5.80 && fabs(eta)>5.75 ) ZimpR=4493.24;
  if ( fabs(eta)<5.75 && fabs(eta)>5.70 ) ZimpR=4640.68;
  if ( fabs(eta)<5.70 && fabs(eta)>5.65 ) ZimpR=5515.86;
  if ( fabs(eta)<5.65 && fabs(eta)>5.60 ) ZimpR=5501.98;
  if ( fabs(eta)<5.45 && fabs(eta)>5.40 ) ZimpR=5051.85;
  if ( fabs(eta)<5.40 && fabs(eta)>5.35 ) ZimpR=4755.5;
  return (ZimpR);
}


double getZimpLH1( double eta ) {
  double ZimpL;
  if ( fabs(eta)<6.40 && fabs(eta)>6.35 ) ZimpL=-7600.;
  if ( fabs(eta)<6.35 && fabs(eta)>6.30 ) ZimpL=-7900.;
  if ( fabs(eta)<6.30 && fabs(eta)>6.25 ) ZimpL=-7405.84;
  if ( fabs(eta)<6.25 && fabs(eta)>6.20 ) ZimpL=-5800.;
  if ( fabs(eta)<6.20 && fabs(eta)>6.15 ) ZimpL=-5710.27;
  if ( fabs(eta)<6.15 && fabs(eta)>6.10 ) ZimpL=-5200.;
  if ( fabs(eta)<6.10 && fabs(eta)>6.05 ) ZimpL=-4900.;
  if ( fabs(eta)<6.05 && fabs(eta)>6.00 ) ZimpL=-5900.;
  if ( fabs(eta)<6.00 && fabs(eta)>5.95 ) ZimpL=-4700.;
  if ( fabs(eta)<5.95 && fabs(eta)>5.90 ) ZimpL=-4700.;
  if ( fabs(eta)<5.90 && fabs(eta)>5.85 ) ZimpL=-4300.;
  if ( fabs(eta)<5.85 && fabs(eta)>5.80 ) ZimpL=-4100.;
  if ( fabs(eta)<5.80 && fabs(eta)>5.75 ) ZimpL=-4200.;
  if ( fabs(eta)<5.75 && fabs(eta)>5.70 ) ZimpL=-4600.;
  if ( fabs(eta)<5.70 && fabs(eta)>5.65 ) ZimpL=-4100.;
  if ( fabs(eta)<5.65 && fabs(eta)>5.60 ) ZimpL=-4200.;
  if ( fabs(eta)<5.45 && fabs(eta)>5.40 ) ZimpL=-4100.;
  if ( fabs(eta)<5.40 && fabs(eta)>5.35 ) ZimpL=-4100.;
  return (ZimpL);
}
double getZimpRH1( double eta ) {
  double ZimpR;
  if ( fabs(eta)<6.40 && fabs(eta)>6.35 ) ZimpR=8634.3;
  if ( fabs(eta)<6.35 && fabs(eta)>6.30 ) ZimpR=7900.;
  if ( fabs(eta)<6.30 && fabs(eta)>6.25 ) ZimpR=6900.;
  if ( fabs(eta)<6.25 && fabs(eta)>6.20 ) ZimpR=6000.;
  if ( fabs(eta)<6.20 && fabs(eta)>6.15 ) ZimpR=5700.;
  if ( fabs(eta)<6.15 && fabs(eta)>6.10 ) ZimpR=5408.34;
  if ( fabs(eta)<6.10 && fabs(eta)>6.05 ) ZimpR=5069.3;
  if ( fabs(eta)<6.05 && fabs(eta)>6.00 ) ZimpR=6165.07;
  if ( fabs(eta)<6.00 && fabs(eta)>5.95 ) ZimpR=4926.65;
  if ( fabs(eta)<5.95 && fabs(eta)>5.90 ) ZimpR=5039.14;
  if ( fabs(eta)<5.90 && fabs(eta)>5.85 ) ZimpR=4679.63;
  if ( fabs(eta)<5.85 && fabs(eta)>5.80 ) ZimpR=4419.66;
  if ( fabs(eta)<5.80 && fabs(eta)>5.75 ) ZimpR=4477.99;
  if ( fabs(eta)<5.75 && fabs(eta)>5.70 ) ZimpR=4952.03;
  if ( fabs(eta)<5.70 && fabs(eta)>5.65 ) ZimpR=4495.14;
  if ( fabs(eta)<5.65 && fabs(eta)>5.60 ) ZimpR=4699.89;
  if ( fabs(eta)<5.45 && fabs(eta)>5.40 ) ZimpR=4675.84;
  if ( fabs(eta)<5.40 && fabs(eta)>5.35 ) ZimpR=4663.83;
  return (ZimpR);
}
double getZimpLH2( double eta ) {
  double ZimpL;
  if ( fabs(eta)<6.40 && fabs(eta)>6.35 ) ZimpL=-10623.7;
  if ( fabs(eta)<6.35 && fabs(eta)>6.30 ) ZimpL=-9762.76;
  if ( fabs(eta)<6.30 && fabs(eta)>6.25 ) ZimpL=-8950.99;
  if ( fabs(eta)<6.25 && fabs(eta)>6.20 ) ZimpL=-9378.61;
  if ( fabs(eta)<6.20 && fabs(eta)>6.15 ) ZimpL=-7986.17;
  if ( fabs(eta)<6.15 && fabs(eta)>6.10 ) ZimpL=-7553.04;
  if ( fabs(eta)<6.10 && fabs(eta)>6.05 ) ZimpL=-7746.9;
  if ( fabs(eta)<6.05 && fabs(eta)>6.00 ) ZimpL=-7110.92;
  if ( fabs(eta)<6.00 && fabs(eta)>5.95 ) ZimpL=-7408.7;
  if ( fabs(eta)<5.95 && fabs(eta)>5.90 ) ZimpL=-6898.26;
  if ( fabs(eta)<5.90 && fabs(eta)>5.85 ) ZimpL=-5854.12;
  if ( fabs(eta)<5.85 && fabs(eta)>5.80 ) ZimpL=-5885.02;
  if ( fabs(eta)<5.80 && fabs(eta)>5.75 ) ZimpL=-5583.3;
  if ( fabs(eta)<5.75 && fabs(eta)>5.70 ) ZimpL=-5780.84;
  if ( fabs(eta)<5.70 && fabs(eta)>5.65 ) ZimpL=-5783.51;
  if ( fabs(eta)<5.65 && fabs(eta)>5.60 ) ZimpL=-5863.39;
  if ( fabs(eta)<5.45 && fabs(eta)>5.40 ) ZimpL=-5599.80;
  if ( fabs(eta)<5.40 && fabs(eta)>5.35 ) ZimpL=-5398.44;
  return (ZimpL);
}
double getZimpRH2( double eta ) {
  double ZimpR;
  if ( fabs(eta)<6.40 && fabs(eta)>6.35 ) ZimpR=7000.;
  if ( fabs(eta)<6.35 && fabs(eta)>6.30 ) ZimpR=7800.;
  if ( fabs(eta)<6.30 && fabs(eta)>6.25 ) ZimpR=7800.;
  if ( fabs(eta)<6.25 && fabs(eta)>6.20 ) ZimpR=8000.;
  if ( fabs(eta)<6.20 && fabs(eta)>6.15 ) ZimpR=6900.;
  if ( fabs(eta)<6.15 && fabs(eta)>6.10 ) ZimpR=6400.;
  if ( fabs(eta)<6.10 && fabs(eta)>6.05 ) ZimpR=6300.;
  if ( fabs(eta)<6.05 && fabs(eta)>6.00 ) ZimpR=5600.;
  if ( fabs(eta)<6.00 && fabs(eta)>5.95 ) ZimpR=5700.;
  if ( fabs(eta)<5.95 && fabs(eta)>5.90 ) ZimpR=4900.;
  if ( fabs(eta)<5.90 && fabs(eta)>5.85 ) ZimpR=4200.;
  if ( fabs(eta)<5.85 && fabs(eta)>5.80 ) ZimpR=4400.;
  if ( fabs(eta)<5.80 && fabs(eta)>5.75 ) ZimpR=4300.;
  if ( fabs(eta)<5.75 && fabs(eta)>5.70 ) ZimpR=4600.;
  if ( fabs(eta)<5.70 && fabs(eta)>5.65 ) ZimpR=4600.;
  if ( fabs(eta)<5.65 && fabs(eta)>5.60 ) ZimpR=4500.;
  if ( fabs(eta)<5.45 && fabs(eta)>5.40 ) ZimpR=4300.;
  if ( fabs(eta)<5.40 && fabs(eta)>5.35 ) ZimpR=4100.;
  return (ZimpR);
}

double getZimpLH3( double eta ) {
  double ZimpL;
  if ( fabs(eta)<6.40 && fabs(eta)>6.35 ) ZimpL=-8681.15;
  if ( fabs(eta)<6.35 && fabs(eta)>6.30 ) ZimpL=-8531.18;
  if ( fabs(eta)<6.30 && fabs(eta)>6.25 ) ZimpL=-7066.61;
  if ( fabs(eta)<6.25 && fabs(eta)>6.20 ) ZimpL=-6464.21;
  if ( fabs(eta)<6.20 && fabs(eta)>6.15 ) ZimpL=-6282.4;
  if ( fabs(eta)<6.15 && fabs(eta)>6.10 ) ZimpL=-6409.5;
  if ( fabs(eta)<6.10 && fabs(eta)>6.05 ) ZimpL=-6168.08;
  if ( fabs(eta)<6.05 && fabs(eta)>6.00 ) ZimpL=-6691.73;
  if ( fabs(eta)<6.00 && fabs(eta)>5.95 ) ZimpL=-6716.16;
  if ( fabs(eta)<5.95 && fabs(eta)>5.90 ) ZimpL=-5594.74;
  if ( fabs(eta)<5.90 && fabs(eta)>5.85 ) ZimpL=-5650.25;
  if ( fabs(eta)<5.85 && fabs(eta)>5.80 ) ZimpL=-5496.15;
  if ( fabs(eta)<5.80 && fabs(eta)>5.75 ) ZimpL=-5686.48;
  if ( fabs(eta)<5.75 && fabs(eta)>5.70 ) ZimpL=-5767.4;
  if ( fabs(eta)<5.70 && fabs(eta)>5.65 ) ZimpL=-5544.53;
  if ( fabs(eta)<5.65 && fabs(eta)>5.60 ) ZimpL=-5890.34;
  if ( fabs(eta)<5.45 && fabs(eta)>5.40 ) ZimpL=-5290.79;
  if ( fabs(eta)<5.40 && fabs(eta)>5.35 ) ZimpL=-5377.84;
  return (ZimpL);
}

double getZimpRH3( double eta ) {
  double ZimpR;
  if ( fabs(eta)<6.40 && fabs(eta)>6.35 ) ZimpR=6700.;
  if ( fabs(eta)<6.35 && fabs(eta)>6.30 ) ZimpR=7200.;
  if ( fabs(eta)<6.30 && fabs(eta)>6.25 ) ZimpR=6500.;
  if ( fabs(eta)<6.25 && fabs(eta)>6.20 ) ZimpR=6000.;
  if ( fabs(eta)<6.20 && fabs(eta)>6.15 ) ZimpR=5900.;
  if ( fabs(eta)<6.15 && fabs(eta)>6.10 ) ZimpR=5600.;
  if ( fabs(eta)<6.10 && fabs(eta)>6.05 ) ZimpR=5500.;
  if ( fabs(eta)<6.05 && fabs(eta)>6.00 ) ZimpR=5600.;
  if ( fabs(eta)<6.00 && fabs(eta)>5.95 ) ZimpR=5500.;
  if ( fabs(eta)<5.95 && fabs(eta)>5.90 ) ZimpR=4300.;
  if ( fabs(eta)<5.90 && fabs(eta)>5.85 ) ZimpR=4200.;
  if ( fabs(eta)<5.85 && fabs(eta)>5.80 ) ZimpR=4100.;
  if ( fabs(eta)<5.80 && fabs(eta)>5.75 ) ZimpR=4400.;
  if ( fabs(eta)<5.75 && fabs(eta)>5.70 ) ZimpR=4500.;
  if ( fabs(eta)<5.70 && fabs(eta)>5.65 ) ZimpR=4300.;
  if ( fabs(eta)<5.65 && fabs(eta)>5.60 ) ZimpR=4600.;
  if ( fabs(eta)<5.45 && fabs(eta)>5.40 ) ZimpR=3800.;
  if ( fabs(eta)<5.40 && fabs(eta)>5.35 ) ZimpR=4000.;
  return (ZimpR);
}

*/
void uat2::PrintPrimTracks(bool minus){
  T2Event*    T2ptr;
  if(minus){
    T2ptr = T2primaryMinus;
    std::cout << "T2- tracks: \n";
  }else{
    T2ptr = T2primaryPlus;
    std::cout << "T2+ tracks: \n";
  };
  unsigned int Ntrks = (T2ptr->TrkEntryX).size();  
  for(unsigned int i_t2=0; i_t2<Ntrks; i_t2++){
      std::cout << "eta: " << T2ptr->TrkEta2.at(i_t2) << ";  phi: " << T2ptr->TrkPhi_RZFit.at(i_t2) << ";  zimp: " << T2ptr->TrkZImpact.at(i_t2) << std::endl;
  }; // end track loop
};

void uat2::PrintAllTracks(bool minus){
  T2Event*    T2ptr = T2tracks;
  if(minus){
    std::cout << "T2- tracks: \n";
  }else{
    std::cout << "T2+ tracks: \n";
  };
  unsigned int Ntrks = (T2ptr->TrkEntryX).size();  
  for(unsigned int i_t2=0; i_t2<Ntrks; i_t2++){
    if( (T2ptr->TrkEntryZ.at(i_t2)<0)==minus)
      std::cout << "eta: " << T2ptr->TrkEta2.at(i_t2) << ";  phi: " << T2ptr->TrkPhi_RZFit.at(i_t2) << ";  zimp: " << T2ptr->TrkZImpact.at(i_t2) << std::endl;
  }; // end track loop
};

void uat2::addT2primaryTrack(unsigned int i, unsigned int cut, const bool fill, const bool info){  
// ONLY TRACK INFO MAKES SENSE HERE DUE TO DIFFERENT SIZE OF TRACK AND HIT COLLECTIONS !    
  T2Event*    T2primptr;
  if(T2tracks->TrkEntryZ.at(i)>0){
    T2primptr = T2primaryPlus;
    if(fill){
      z_impact_prim_plus_h[cut]->Fill(T2tracks->TrkZImpact.at(i));
    };
  }else{
    T2primptr = T2primaryMinus;
    if(fill){
      z_impact_prim_minus_h[cut]->Fill(T2tracks->TrkZImpact.at(i));
    };
  }; 

 
//  T2primptr->Pad_row.push_back(T2tracks->Pad_row.at(i));             //strip row   
//  T2primptr->Pad_col.push_back(T2tracks->Pad_col.at(i));             //strip column   
//  T2primptr->Pad_det.push_back(T2tracks->Pad_det.at(i));             //symbolic id of the detector containind the pad
    
//  (T2primptr->Strip_row).push_back(T2tracks->Strip_row.at(i));           //strip row
//  (T2primptr->Strip_col).push_back(T2tracks->Strip_col.at(i));           //strip column
//  (T2primptr->Strip_det).push_back(T2tracks->Strip_det.at(i));           //symbolic id of the detector containind the strip
  
//(T2primptr->TrkPrimaryGeneratorEta).push_back(T2tracks->TrkPrimaryGeneratorEta);
  (T2primptr->TrkEta_XY ).push_back(T2tracks->TrkEta_XY .at(i));            
  (T2primptr->TrkZmin_XY).push_back(T2tracks->TrkZmin_XY.at(i));  
  (T2primptr->TrkRmin_XY).push_back(T2tracks->TrkRmin_XY.at(i));  
  (T2primptr->TrkZImpact).push_back(T2tracks->TrkZImpact.at(i));
  
  (T2primptr->TrkAx).push_back(T2tracks->TrkAx.at(i));            // slope of the track projection in the XZ plane
  (T2primptr->TrkAy).push_back(T2tracks->TrkAy.at(i));            // slope of the track projection in the YZ plane
  (T2primptr->TrkX0).push_back(T2tracks->TrkX0.at(i));            // X at Z=0 for the XZ projected track  
  (T2primptr->TrkY0 ).push_back(T2tracks->TrkY0.at(i));            // Y at Z=0 for the XZ projected track
  (T2primptr->TrkPhi).push_back(T2tracks->TrkPhi.at(i));           // Trk Phi (XY fit)

  (T2primptr->TrkChi2XProb).push_back(T2tracks->TrkChi2XProb.at(i));          //Chi2-X probability (goodness of the XZ projection fit)
  (T2primptr->TrkChi2YProb).push_back(T2tracks->TrkChi2YProb.at(i));          //Chi2-Y probability (goodness of the YZ projection fit)

  (T2primptr->TrkClass1HitCounter).push_back(T2tracks->TrkClass1HitCounter.at(i));   //Number of class1 Hit in the Trk
  (T2primptr->TrkHitCounter).push_back(T2tracks->TrkHitCounter.at(i));    //Number of class1 + cluster Pad hits in the Trk
  
                    
  (T2primptr->TrkThetaR_RZFit).push_back(T2tracks->TrkThetaR_RZFit.at(i));// Trk Polar angle obtained tracking in the Rz plane
  (T2primptr->TrkEta_RZFit).push_back(T2tracks->TrkEta_RZFit.at(i));   // Trk Eta obtained tracking in the Rz plane 
  (T2primptr->TrkPhi_RZFit).push_back(T2tracks->TrkPhi_RZFit.at(i));   // Trk Phi obtained with a constant fit.
  (T2primptr->TrkZ0_RZFit ).push_back(T2tracks->TrkZ0_RZFit .at(i));    // Crossing Point between Trk and Z Axis obtained tracking in the Rz plane 
  (T2primptr->TrkBX_RZFit ).push_back(T2tracks->TrkBX_RZFit .at(i));    // X0 @ Z=0 obtained tracking in the Rz plane 
  (T2primptr->TrkBY_RZFit ).push_back(T2tracks->TrkBY_RZFit .at(i));    // X0 @ Z=0 obtained tracking in the Rz plane 
  
  T2primptr->NumPadCluH0 = T2tracks->NumPadCluH0;
  T2primptr->NumPadCluH1 = T2tracks->NumPadCluH1;
  T2primptr->NumPadCluH2 = T2tracks->NumPadCluH2;
  T2primptr->NumPadCluH3 = T2tracks->NumPadCluH3;
  
  (T2primptr->TrkEta2).push_back(T2tracks->TrkEta2.at(i));
  
  

  (T2primptr->TrkChiProb).push_back(T2tracks->TrkChiProb.at(i));             //Make sense only with T2TrackProducer3
/*
  (T2primptr->Trknumpadonly).push_back(T2tracks->Trknumpadonly.at(i));   
  (T2primptr->Trknumstriponly).push_back(T2tracks->Trknumstriponly.at(i));   
  (T2primptr->Trknumhitall).push_back(T2tracks->Trknumhitall.at(i));   
  (T2primptr->Trknumhitacl1).push_back(T2tracks->Trknumhitacl1.at(i));   
*/

(T2primptr->ProbChi2R_rz).push_back(T2tracks->ProbChi2R_rz.at(i));
(T2primptr->Chi2Rreduced_rz).push_back(T2tracks->Chi2Rreduced_rz.at(i));

/*
  (T2primptr->HitPhi).push_back(T2tracks->HitPhi.at(i));      // Phi position of all the Hits (deg)
  (T2primptr->HitR).push_back(T2tracks->HitR.at(i));        // R position of all the Hits (mm)
  (T2primptr->HitType).push_back(T2tracks->HitType.at(i));     // 0-> only pad; 1-> only strip 2->Class 1 Hit (superimposition Pad/Strip)
  (T2primptr->HitNumPad).push_back(T2tracks->HitNumPad.at(i));   // Cluster Pad Size 
  (T2primptr->HitNumStrip).push_back(T2tracks->HitNumStrip.at(i)); // Cluster Strip Size
  (T2primptr->HitID).push_back(T2tracks->HitID);
  */

  (T2primptr->TrkEntryX).push_back(T2tracks->TrkEntryX.at(i));   //Trk Entry point X
  (T2primptr->TrkEntryY).push_back(T2tracks->TrkEntryY.at(i));   //Trk Entry point Y 
  (T2primptr->TrkEntryZ).push_back(T2tracks->TrkEntryZ.at(i));   //Trk Entry point Z
  (T2primptr->TrkExitX).push_back(T2tracks->TrkExitX.at(i));    //Trk Exit point X
  (T2primptr->TrkExitY).push_back(T2tracks->TrkExitY.at(i));    //Trk Exit point Y  
  (T2primptr->TrkExitZ).push_back(T2tracks->TrkExitZ.at(i));    //Trk Exit point Z
  
 
   
  //T2primptr->run_no = T2tracks->run_no; 
  //T2primptr->ev_no  = T2tracks->ev_no;
  //T2primptr->timestamp = T2tracks->timestamp;
  
  
  //(T2primptr->Pad_noise).push_back(T2tracks->Pad_noise.at(i));
};



void uat2::resetT2primary(bool plus){
  T2Event*    T2primptr;
  if(plus){
    T2primptr = T2primaryPlus;
  }else{
    T2primptr = T2primaryMinus;
  }; 
  
 (T2primptr->Pad_row).clear();             //strip row   
 (T2primptr->Pad_col).clear();             //strip column   
 (T2primptr->Pad_det).clear();             //symbolic id of the detector containind the pad
 
 (T2primptr->Strip_row).clear();           //strip row
 (T2primptr->Strip_col).clear();           //strip column
 (T2primptr->Strip_det).clear();           //symbolic id of the detector containind the strip

 (T2primptr->TrkPrimaryGeneratorEta).clear();
 (T2primptr->TrkEta_XY ).clear();            
  (T2primptr->TrkZmin_XY).clear();  
  (T2primptr->TrkRmin_XY).clear();  
  (T2primptr->TrkZImpact).clear();

  (T2primptr->TrkAx).clear();            // slope of the track projection in the XZ plane
  (T2primptr->TrkAy).clear();            // slope of the track projection in the YZ plane
  (T2primptr->TrkX0).clear();            // X at Z=0 for the XZ projected track  
  (T2primptr->TrkY0 ).clear();            // Y at Z=0 for the XZ projected track
  (T2primptr->TrkPhi).clear();           // Trk Phi (XY fit)
                           
  (T2primptr->TrkChi2XProb).clear();          //Chi2-X probability (goodness of the XZ projection fit)
  (T2primptr->TrkChi2YProb).clear();          //Chi2-Y probability (goodness of the YZ projection fit)
  (T2primptr->TrkClass1HitCounter).clear();   //Number of class1 Hit in the Trk
  (T2primptr->TrkHitCounter).clear();    //Number of class1 + cluster Pad hits in the Trk
                              
                              
  (T2primptr->TrkThetaR_RZFit).clear();// Trk Polar angle obtained tracking in the Rz plane
  (T2primptr->TrkEta_RZFit).clear();   // Trk Eta obtained tracking in the Rz plane 
  (T2primptr->TrkPhi_RZFit).clear();   // Trk Phi obtained with a constant fit.
  (T2primptr->TrkZ0_RZFit ).clear();    // Crossing Point between Trk and Z Axis obtained tracking in the Rz plane 
  (T2primptr->TrkBX_RZFit ).clear();    // X0 @ Z=0 obtained tracking in the Rz plane 
  (T2primptr->TrkBY_RZFit ).clear();    // X0 @ Z=0 obtained tracking in the Rz plane 
                           
  (T2primptr->NumPadCluH0) = 0;
  (T2primptr->NumPadCluH1) = 0;
  (T2primptr->NumPadCluH2) = 0;
  (T2primptr->NumPadCluH3) = 0;

  (T2primptr->TrkNumHitInH0).clear();
  (T2primptr->TrkNumHitInH1).clear();
  (T2primptr->TrkNumHitInH2).clear();
  (T2primptr->TrkNumHitInH3).clear();
  (T2primptr->TrkEta2).clear();
  
 (T2primptr->TrkChiProb).clear();             //Make sense only with T2TrackProducer3
 (T2primptr->Trknumpadonly).clear();   
 (T2primptr->Trknumstriponly).clear();   
 (T2primptr->Trknumhitall).clear();   
 (T2primptr->Trknumhitacl1).clear();   
 (T2primptr->ProbChi2R_rz).clear(); 
  (T2primptr->Chi2Rreduced_rz).clear();
  
 (T2primptr->HitPhi).clear();      // Phi position of all the Hits (deg)
 (T2primptr->HitR).clear();        // R position of all the Hits (mm)
 (T2primptr->HitType).clear();     // 0-> only pad; 1-> only strip 2->Class 1 Hit (superimposition Pad/Strip)
 (T2primptr->HitNumPad).clear();   // Cluster Pad Size 
 (T2primptr->HitNumStrip).clear(); // Cluster Strip Size
 (T2primptr->HitID).clear();

 (T2primptr->TrkEntryX).clear();   //Trk Entry point X
 (T2primptr->TrkEntryY).clear();   //Trk Entry point Y 
 (T2primptr->TrkEntryZ).clear();   //Trk Entry point Z
 (T2primptr->TrkExitX).clear();    //Trk Exit point X
 (T2primptr->TrkExitY).clear();    //Trk Exit point Y  
 (T2primptr->TrkExitZ).clear();    //Trk Exit point Z



 T2primptr->run_no = 0; 
 T2primptr->ev_no  = 0;
 T2primptr->timestamp = 0;

 
 (T2primptr->Pad_noise).clear();
};


