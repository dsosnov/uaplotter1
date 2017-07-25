/*
 *
 */

#include "uaforward.h"
#include "uathresholds.h"

#include "iostream"
#include "stdlib.h"

ClassImp(uaforward)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uaforward::uaforward(TChain* tree, 
		     TDirectory * dir,
		     const bool cmstotem, 
		     const bool ZDC56,
		     const short          int MC, 
		     const short unsigned int Ncuts):
      uabase(cmstotem, MC, Ncuts, dir), zdc56(ZDC56)
{
  ZDCdigis = 0;
  FSCdigis = 0;
  if(mc<=0){
    if(tree_combined_flag){
	tree->SetBranchAddress("cmszdcDigisUA", &ZDCdigis);
	tree->SetBranchAddress("cmsfscDigisUA", &FSCdigis);
    }else{
	tree->SetBranchAddress("zdcDigis",      &ZDCdigis);
	tree->SetBranchAddress("fscDigis",      &FSCdigis);
    };
  };
  startFSCts=3;
  if(startFSCts>FSC_NNOTS){
    std::cout << "uaforward::uaforward(): Number of FSC noise TS is larger than first signal TS!!!"<< std::endl;
    exit(0);
  };
  create_histos();
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uaforward::~uaforward()
{
  std::cout << "uaforward::~uaforward(): deleting " 
	    << h1D->size() << "+" << h2D->size() << " histos" << std::endl;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uaforward::ProceedEvent(short unsigned int cut, const bool fill, const bool info)
{   
  if(fill && (cut>n_cuts) ){
    std::cout << "uaforward::ProceedEvent: required cut number is larger that possible, please define larger uaforward::n_cut!\n"; 
    return false;
  }; 

  memset(Etotal, 0, sizeof(Etotal));
  memset(EEM,    0, sizeof(EEM));

  if(mc>0)
    return false;

  // ZDC is configured to have firstSample be TS5 in runs 210498-210676, and TS4 in runs 210737-211831 
  for(std::vector<MyZDCDigi>::iterator it=ZDCdigis->begin(); it!=ZDCdigis->end(); ++it){      
    float            sig = 0;
    unsigned int minsize = 6;
    if(zdc56)
      minsize = 7;

    std::vector<float> q = (*it).digifC;
    if(q.size()<minsize){
      std::cout << "uaforward::ProceedEvent: number of ZDC digi timeslices < " << minsize << "; doing nothing\n";
    }else{
      if(zdc56){sig = q.at(5) + q.at(6) - q.at(3) - q.at(4);}
      else     {sig = q.at(4) + q.at(5) - q.at(2) - q.at(3);}; // checked with 210885
//       std::cout << (*it).channelId << "\t\t";
//       for (short unsigned int ts=0; ts<10; ts++)
// 	std::cout << q.at(ts) << "  ";
//       std::cout << std::endl;
    };

    // EM section
    if( (*it).section==1 ) sig*=0.1;

    short unsigned int ind = int((*it).side<0); // side = -1 for +Z! true side=-1=>1 +z=1, -z=0
    Etotal[ind] += sig;
    if((*it).section==1) 
      EEM[ind]  += sig;    
  }; // end zdc digy loop




  memset(FSCm_TS, 0, sizeof(FSCm_TS));
  memset(FSCm_No, 0, sizeof(FSCm_No));
  memset(FSCm_Si, 0, sizeof(FSCm_Si));
  memset(FSCm_flag, false, sizeof(FSCm_flag));
  FSCm_No8=0;
  FSCm_Si8=0;
  goodFSC = true;
  // loop FSC digis
  for(std::vector<MyFSCDigi>::iterator it=FSCdigis->begin(); it!=FSCdigis->end(); ++it){      
    float sig = 0;
    float noi = 0;

    std::vector<float> q = (*it).digifC;

    short unsigned int qsize=10<q.size()?10:q.size();
    unsigned int minsize=10;

    if((qsize<minsize) || (FSC_NSITS+startFSCts)>qsize){
      std::cout << "uaforward::ProceedEvent: number of FSC digi timeslices < " << minsize << "; doing nothing\n";
    }else{
      if((*it).side>0){
	
	unsigned int ch = ((*it).channelId-26);
	for( short unsigned int ts=0; ts<FSC_NNOTS; ts++)
	  noi+=q.at(ts);
	for( short unsigned int ts=startFSCts; ts<FSC_NSITS+startFSCts; ts++)
	  sig+=q.at(ts);
	for(short unsigned int j=0; j<qsize; j++)
	  FSCm_TS[ch][j]  = q.at(j);
	
	if( (1200 < noi) && (noi < 2800) ){ //!!! TBD 
	  FSCm_No[ch] = noi;
	  FSCm_Si[ch] = sig;
	}else if (ch<6) {
	  goodFSC = false;
	  //std::cout << "ch=" << ch << std::endl;
	  //for(short unsigned int ts=0; ts<10; ts++)
	  //  std::cout << FSCm_TS[ch][ts] << "   ";
	  //std::cout << "\n";
	};
	//std::cout << FSCm_No[ch] << std::endl;
	if(sig>FSC_SITHR[ch]) FSCm_flag[ch]=true;
	
	if(ch<6){
	  FSCm_No8+=noi;
	  FSCm_Si8+=sig;
	};
      };
    };

  }; // end FSC digi loop

  if(info) PrintEventInfo(true);
  if(fill) FillLastEvent(cut);
  return true;
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uaforward::PrintEventInfo(bool detailed){
  std::cout << "uaforward::PrintEventInfo:\n";
  if(detailed){
    std::cout << "\tZDC- : "     << Etotal[0] << " EM("<< EEM[0] 
	      << ");\tZDC+: "  << Etotal[1] << " EM("<< EEM[1] << ")"<< std::endl;
    std::cout << "\tFSC-:";
    if(goodFSC){
      std::cout << " GOOD\n";
    }else{
      std::cout << " BAD\n";
    };
    std::cout << "\tnoise\t";
    for(short unsigned int ch=0; ch<8; ch++)
      std::cout << FSCm_No[ch] << "\t";
    std::cout << std::endl;
    std::cout << "\tsignal\t";
    for(short unsigned int ch=0; ch<8; ch++)
      std::cout << FSCm_Si[ch] << "\t";
    std::cout << std::endl;
    std::cout << "\tflags\t";
    for(short unsigned int ch=0; ch<8; ch++)
      std::cout << FSCm_flag[ch] << "\t";
    std::cout << std::endl;
    for(short unsigned int ch=0; ch<8; ch++){
      std::cout << "\t" << ch << "\t";
      for(short unsigned int ts=0; ts<10; ts++){
	std::cout << FSCm_TS[ch][ts] << "\t";
      };
      std::cout << std::endl;
    };
    std::cout << "\t\t=> total (6ch) noise: " << FSCm_No8 << "; signal: " << FSCm_Si8 << std::endl;
  }else{
    std::cout << "\tZDC- : "    << Etotal[0] 
	      << ";\tZDC+: "  << Etotal[1] << std::endl;    
    std::cout << "\tFSC-: "; 
    if(goodFSC){
      std::cout << " GOOD\n";
    }else{
      std::cout << " BAD\n";
    };
    std::cout << "\tnoise\t";
    for(short unsigned int ch=0; ch<8; ch++)
      std::cout << FSCm_No[ch] << "\t";
    std::cout << std::endl;
    std::cout << "\tsignal\t";
    for(short unsigned int ch=0; ch<8; ch++)
      std::cout << FSCm_Si[ch] << "\t";
    std::cout << std::endl;
        std::cout << "\tflags\t";
    for(short unsigned int ch=0; ch<8; ch++)
      std::cout << FSCm_flag[ch] << "\t";
    std::cout << std::endl;
    std::cout << "\t\t=> total (6ch) noise: " << FSCm_No8 << "; signal: " << FSCm_Si8 << std::endl;
  };
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uaforward::FillLastEvent(const short unsigned int cut)
{
  if( cut>=n_cuts ){
    std::cout << "uaforward::FillLastEvent: required cut number is larger that possible, do nothing. Please define larger uaforward::n_cut!\n"; 
    return false;
  }; 
  ZDCd_Etot_minus_h[cut]->Fill(Etotal[0]);
  ZDCd_Etot_plus_h[cut] ->Fill(Etotal[1]);
  ZDCd_EM_minus_h[cut]->Fill(EEM[0]);
  ZDCd_EM_plus_h[cut] ->Fill(EEM[1]);

  if(mc<=0 and goodFSC){ // TBD check for MC and noise
    short unsigned int fsc_novthr = 0;
    for(short unsigned int fscch=0; fscch<8; fscch++){
      TH1F ** h = FSCm_TS_ha[fscch];
      //std::cout << "=> " << h << std::endl;
      for(short unsigned int ts=0; ts<10; ts++){
	h[cut]->Fill(ts, FSCm_TS[fscch][ts]);
      };
      FSCm_No_ha[fscch][cut]->Fill(FSCm_No[fscch]);
      FSCm_Si_ha[fscch][cut]->Fill(FSCm_Si[fscch]); 
      if(FSCm_flag[fscch]) fsc_novthr++;
    };
    FSCm_N_h[cut]->Fill(fsc_novthr);
    FSCm_No8_h[cut]->Fill(FSCm_No8);
    FSCm_Si8_h[cut]->Fill(FSCm_Si8);
    ZDCm_FSCmSi8_h[cut]->Fill(FSCm_Si8-FSCm_No8, Etotal[0]);
    ZDCm_FSCmN_h[cut]->Fill(fsc_novthr, Etotal[0]);
  };
  return true;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uaforward::NormalizeFSCts()
{
  if (mc>0) return;
  for(short unsigned int fscch=0; fscch<8; fscch++){
    TH1F ** h = FSCm_TS_ha[fscch];
    for(short unsigned int hn=0; hn<n_each_h1D; hn++){
      h[hn]->Sumw2();
      h[hn]->Scale(10./h[hn]->GetEntries());
    };
    h = FSCm_No_ha[fscch];
    for(short unsigned int hn=0; hn<n_each_h1D; hn++){
      h[hn]->Sumw2();
      h[hn]->Scale(1./h[hn]->GetEntries());
    };
    h = FSCm_Si_ha[fscch];
    for(short unsigned int hn=0; hn<n_each_h1D; hn++){
      h[hn]->Sumw2();
      h[hn]->Scale(1./h[hn]->GetEntries());
    };
  };

  for(short unsigned int hn=0; hn<n_each_h1D; hn++){
    FSCm_No8_h[hn]->Sumw2();
    FSCm_No8_h[hn]->Scale(1./FSCm_No8_h[hn]->GetEntries());
    FSCm_Si8_h[hn]->Sumw2();
    FSCm_Si8_h[hn]->Scale(1./FSCm_Si8_h[hn]->GetEntries());
    FSCm_N_h[hn]->Sumw2();
    FSCm_N_h[hn]->Scale(1./FSCm_N_h[hn]->GetEntries());
  };  
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uaforward::FillZDCWithMCtruth(const short unsigned int cut, const double E[2], const double EM[2])
{
  ZDCd_Etot_minus_h[cut]->Fill(E[0]);
  ZDCd_Etot_plus_h[cut] ->Fill(E[1]);
  ZDCd_EM_minus_h[cut]->Fill(EM[0]);
  ZDCd_EM_plus_h[cut] ->Fill(EM[1]);
}





//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uaforward::create_histos()
{
  TString title1, title2;

  n_each_h1D = n_cuts;
  n_each_h2D = n_cuts;

  ZDCd_Etot_minus_h = new TH1F * [n_each_h1D];
  ZDCd_Etot_plus_h  = new TH1F * [n_each_h1D];
  ZDCd_EM_minus_h   = new TH1F * [n_each_h1D];
  ZDCd_EM_plus_h    = new TH1F * [n_each_h1D];

  for(short unsigned int i=0; i<n_each_h1D; i++){
    title1 = "ZDCd_Etot_minus_h["; title1+=i; title1+="]";
    title2 = title1; title2+="; E^{-}_{tot} [GeV]";
    ZDCd_Etot_minus_h[i] = new TH1F(title1.Data(), title2.Data(), 4200, -2000, 40000);
    ZDCd_Etot_minus_h[i]->SetDirectory(directory); 

    title1 = "ZDCd_Etot_plus_h["; title1+=i; title1+="]";
    title2 = title1; title2+="; E^{+}_{tot} [GeV]";
    ZDCd_Etot_plus_h[i] = new TH1F(title1.Data(), title2.Data(), 4200, -20000, 400000);
    ZDCd_Etot_plus_h[i]->SetDirectory(directory);

    title1 = "ZDCd_EM_minus_h["; title1+=i; title1+="]";
    title2 = title1; title2+="; E^{-}_{tot} [GeV]";
    ZDCd_EM_minus_h[i] = new TH1F(title1.Data(), title2.Data(), 4200, -2000, 40000);
    ZDCd_EM_minus_h[i]->SetDirectory(directory);

    title1 = "ZDCd_EM_plus_h["; title1+=i; title1+="]";
    title2 = title1; title2+="; E^{+}_{tot} [GeV]";
    ZDCd_EM_plus_h[i] = new TH1F(title1.Data(), title2.Data(), 4200, -20000, 400000);  
    ZDCd_EM_plus_h[i]->SetDirectory(directory);
  }; 
  h1D->push_back(ZDCd_Etot_minus_h);
  h1D->push_back(ZDCd_Etot_plus_h);
  h1D->push_back(ZDCd_EM_minus_h);
  h1D->push_back(ZDCd_EM_plus_h);  

  if (mc<=0){ // TBD check
    for(short unsigned int fscch=0; fscch<8; fscch++){
      FSCm_TS_ha[fscch] = new TH1F * [n_each_h1D];
      FSCm_No_ha[fscch] = new TH1F * [n_each_h1D];
      FSCm_Si_ha[fscch] = new TH1F * [n_each_h1D];
      TH1F ** h  = FSCm_TS_ha[fscch];
      TH1F ** h1 = FSCm_No_ha[fscch];
      TH1F ** h2 = FSCm_Si_ha[fscch];
      for(short unsigned int i=0; i<n_each_h1D; i++){
	title1 = "FSCm_TS_ha"; title1+=fscch; title1+="["; title1+=i; title1+="]";
	title2 = title1; title2+=";TS#";
	h[i] = new TH1F(title1.Data(), title2.Data(), 10, 0, 10);
	h[i]->SetDirectory(directory); 
	title1 = "FSCm_No_ha"; title1+=fscch; title1+="["; title1+=i; title1+="]";
	title2 = title1; title2+=";q3 [fC]";
	h1[i] = new TH1F(title1.Data(), title2.Data(), 400, 0, 40000);
	h1[i]->SetDirectory(directory); 
	title1 = "FSCm_Si_ha"; title1+=fscch; title1+="["; title1+=i; title1+="]";
	title2 = title1; title2+=";q3 [fC]";
	h2[i] = new TH1F(title1.Data(), title2.Data(), 400, 0, 40000);
	h2[i]->SetDirectory(directory); 
      };
      h1D->push_back(h);
      h1D->push_back(h1);
      h1D->push_back(h2);
    };
    FSCm_N_h = new TH1F * [n_each_h1D];  
    FSCm_No8_h = new TH1F * [n_each_h1D];  
    FSCm_Si8_h = new TH1F * [n_each_h1D];  

    ZDCm_FSCmSi8_h = new TH2F * [n_each_h2D];
    ZDCm_FSCmN_h   = new TH2F * [n_each_h2D];
    for(short unsigned int i=0; i<n_each_h1D; i++){
      title1 = "FSCm_N_h["; title1+=i; title1+="]";
      title2 = title1; title2+="; number of FSC channels";
      FSCm_N_h[i] = new TH1F(title1.Data(), title2.Data(), 10, -1,9);
      FSCm_N_h[i]->SetDirectory(directory); 

      title1 = "FSCm_No8_h["; title1+=i; title1+="]";
      title2 = title1; title2+=";q6 [fC]";
      FSCm_No8_h[i] = new TH1F(title1.Data(), title2.Data(), 400, 0, 240000);
      FSCm_No8_h[i]->SetDirectory(directory); 

      title1 = "FSCm_Si8_h["; title1+=i; title1+="]";
      title2 = title1; title2+=";q6 [fC]";
      FSCm_Si8_h[i] = new TH1F(title1.Data(), title2.Data(), 400, 0, 240000);
      FSCm_Si8_h[i]->SetDirectory(directory); 

      title1 = "ZDCm_FSCmSi8_h["; title1+=i; title1+="]";
      title2 = title1; title2+=";q6(FSC-) [fC];E^{-}_{tot}(ZDC) [GeV]";
      ZDCm_FSCmSi8_h[i] = new TH2F(title1.Data(), title2.Data(), 200, -3000, 197000, 420, -2000, 40000);
      ZDCm_FSCmSi8_h[i]->SetDirectory(directory); 

      title1 = "ZDCm_FSCmN_h["; title1+=i; title1+="]";
      title2 = title1; title2+=";number of FSC channels;E^{-}_{tot}(ZDC) [GeV]";
      ZDCm_FSCmN_h[i] = new TH2F(title1.Data(), title2.Data(), 10, -1, 9, 420, -2000, 40000);
      ZDCm_FSCmN_h[i]->SetDirectory(directory); 
    };

    h1D->push_back(FSCm_N_h);
    h1D->push_back(FSCm_No8_h);
    h1D->push_back(FSCm_Si8_h);
    h2D->push_back(ZDCm_FSCmSi8_h);
    h2D->push_back(ZDCm_FSCmN_h);
  };  
}
