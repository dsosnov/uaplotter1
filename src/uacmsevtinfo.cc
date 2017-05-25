/*
 *
 */

#include "uacmsevtinfo.h"

#include "iostream"
#include "stdlib.h"

ClassImp(uacmsevtinfo)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uacmsevtinfo::uacmsevtinfo(TChain * tree, 
		  TDirectory * dir,    
		  const bool cmstotem, 
		  const short          int MC,  
		  const short unsigned int Ncuts 
		  ):uabase(cmstotem, MC, Ncuts, dir)
{
  CMSevtInfo  = 0;
  CMStrigInfo = 0;
  CMSHLT      = 0;
  if(tree_combined_flag){
    tree->SetBranchAddress("cmsEvtUA",               &CMSevtInfo);  
    tree->SetBranchAddress("cmsTrigUA",              &CMStrigInfo);  
    tree->SetBranchAddress("cmsHLTTrigUA",           &CMSHLT);  
  }else{
    tree->SetBranchAddress("evtId",                  &CMSevtInfo);  
    tree->SetBranchAddress("L1TrigOld",              &CMStrigInfo);  
    tree->SetBranchAddress("HLTrig",                 &CMSHLT);  
  };
  create_histos();
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

uacmsevtinfo::~uacmsevtinfo()
{
  std::cout << "uacmsevtinfo::~uacmsevtinfo(): deleting " 
	    << h1D->size() << "+" << h2D->size() << " histos" << std::endl;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void uacmsevtinfo::PrintEventInfo(const bool detailed)
{ 
  if(detailed){
    std::cout << "uacmsevtinfo::PrintEventInfo: L1:\n\t";
    for(short unsigned int b=0; b<128; b++)
      if(CMStrigInfo->PhysTrigWord[b])
	std::cout << b << "   ";
    std::cout << "\nuacmsevtinfo::PrintEventInfo: tt:\n\t";
    for(short unsigned int b=0; b<64; b++)
      if(CMStrigInfo->TechTrigWord[b])
	std::cout << b << "   ";
    std::cout << std::endl;
  }else{
    std::cout << "\nuacmsevtinfo::PrintEventInfo: L1 and tt of interest:\n\t";
    for(short unsigned int b=0; b<3; b++)
      std::cout << ALGO_TRIGGER_ARRAY[b] << "\t";
    for(short unsigned int b=0; b<5; b++)
      std::cout << "\t" << TT_TRIGGER_ARRAY[b] ;
    std::cout << "\n\t";
    for(short unsigned int b=0; b<3; b++)
      std::cout << CMStrigInfo->PhysTrigWord[ALGO_TRIGGER_ARRAY[b]] << "\t";
    for(short unsigned int b=0; b<5; b++)
      std::cout << "\t" << CMStrigInfo->TechTrigWord[TT_TRIGGER_ARRAY[b]] ;    
    std::cout << std::endl;
  };
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool uacmsevtinfo::FillLastEvent(const short unsigned int cut)
{
  if( cut>=n_cuts ){
    std::cout << "uacmsevtinfo::FillLastEvent: required cut number is larger that possible, do nothing. Please define larger uaforward::n_cut!\n"; 
    return false;
  }; 
  run_vs_bx_h[cut]->Fill(CMSevtInfo->Bunch,    CMSevtInfo->Run);
  run_vs_ls_h[cut]->Fill(CMSevtInfo->LumiSect, CMSevtInfo->Run);
  
  for(unsigned int tb=0; tb<3; tb++)
    if(CMStrigInfo->PhysTrigWord[ALGO_TRIGGER_ARRAY[tb]])
      triggers_h[cut]->Fill(tb);

  for(unsigned int tb=0; tb<5; tb++)
    if(CMStrigInfo->TechTrigWord[TT_TRIGGER_ARRAY[tb]])
      triggers_h[cut]->Fill(tb+5);
  return true;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uacmsevtinfo::ProceedEvent(const short unsigned int cut, const bool fill, const bool info)
{
  if( (!fill) ){
    std::cout << "uacmsevtinfo::ProceedEvent: fill flag is false, that makes no sense, do nothing!\n"; 
    return false;
  }; 
  FillLastEvent(cut);
  if(info) PrintEventInfo(true);
  return true;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uacmsevtinfo::CheckHLT(const char* path)
{
  bool on = false;
  for(std::map<string,bool>::iterator hlt=CMSHLT->HLTmap.begin(); hlt!=CMSHLT->HLTmap.end(); ++hlt){
    //std::cout << (*hlt).first << "\t" << (*hlt).second << std::endl;
    std::size_t found = (*hlt).first.find(path);
    if (found!=std::string::npos){
      on =  !((*hlt).second);
    };
  };
  return on;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uacmsevtinfo::create_histos()
{
  TString title1, title2;
  
  n_each_h1D  = n_cuts;
  n_each_h2D  = n_cuts;
  triggers_h  = new TH1F* [n_each_h1D];
  run_vs_bx_h = new TH2F* [n_each_h2D];
  run_vs_ls_h = new TH2F* [n_each_h2D];

  for(short unsigned int i=0;i<n_each_h1D; i++){
    title1 = "triggers_h["; title1+=i; title1+="]";
    title2 = title1; title2+="; trigger xxx";
    triggers_h[i] = new TH1F(title1.Data(),title2.Data(), 12, 0,12); 
    for(short unsigned int b=1; b<11; b++){
      TString label="";
      if(b<4) {
	label = "L1_"; label+=ALGO_TRIGGER_ARRAY[b-1];
      }else if (b>5){
	label = "tt_"; label+=TT_TRIGGER_ARRAY[b-6];
      };
      triggers_h[i]->GetXaxis()->SetBinLabel(b,label.Data());
    };
    triggers_h[i]->SetDirectory(directory);     
  };
  h1D->push_back(triggers_h);
  
  
  for(short unsigned int i=0;i<n_each_h2D; i++){  
    title1 = "run_vs_bx_h["; title1+=(i); title1+="]";
    title2 = title1; title2+=";BCN";
    run_vs_bx_h[i] = new TH2F(title1.Data(),title2.Data(), 3100, -1,3099, 1100, 210490, 211590);
    run_vs_bx_h[i]->SetDirectory(directory); 
    
    title1 = "run_vs_ls_h["; title1+=(i); title1+="]";
    title2 = title1; title2+=";LS";
    run_vs_ls_h[i] = new TH2F(title1.Data(),title2.Data(), 1000, -2, 1998, 1100, 210490, 211590);
    run_vs_ls_h[i]->SetDirectory(directory); 
  };
  h2D->push_back(run_vs_bx_h);
  h2D->push_back(run_vs_ls_h);
}
