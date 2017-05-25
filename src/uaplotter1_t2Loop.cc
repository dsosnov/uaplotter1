#include "uaplotter1.h"
#include "uathresholds.h"

#include "TString.h"
#include "TMath.h"

#include "iostream"
#include "stdlib.h"


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //! t2 studies on tt53 trigger
/*!
\param evts number of events to be prodused
\return number of events of the good evts
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int uaplotter1::t2Loop(const int evts){ 
  if (mc<0){
    std::cout << "can't run t2Loop on noise files...\n";
    return 0;
  }else if ( (mc==0) && (!tree_combined_flag) ){
    std::cout << "can't run t2Loop on CMS-only files...\n";
    return 0;    
  };
  
  int stat = chainTree->GetEntries(); 
  unsigned int nevts = evts;
  if(evts==-1 || evts>stat) nevts=stat;
  std::cout << "Total stat = " << stat << std::endl;
  std::cout << "uaplotter1::t2Loop() for " << nevts << " events\n";

  unsigned int tt53         =0;
  unsigned int tt52         =0;
  unsigned int al99         =0;
  unsigned int ttBOTH       =0;
  unsigned int BPTX         =0;
  unsigned int noPU         =0; // number of vertecis<2
  unsigned int ND           =0;
  unsigned int EL           =0;
  unsigned int SDm          =0;
  unsigned int SDp          =0;
  
  unsigned int kevt = 0;
  
  //==============================================================================
  for(long unsigned int i = 0; i<nevts; i++){  
    unsigned int kevt_current = i/1000;
    if(kevt_current>kevt){
      kevt = kevt_current;
      std::cout << ">>> proceeding " << kevt << std::endl;
    };
    chainTree->GetEntry(i);
    current_event = i;
    
    
    memset(sd_flag_central, 0, sizeof(sd_flag_central));
    memset(sd_flag_total,   0, sizeof(sd_flag_total));
    bool MB = false;
    //std::cout << ProceedTrigger(52,true) << "  " << ProceedTrigger(53,true) << std::endl;
    if(ProceedTrigger(53, true)){
      MB = true;
      tt53++;
    };
    
    if(ProceedTrigger(53, true)){
      if( MB ) ttBOTH++;
      //al99++;
      if( !CMSevtinfo->GetTechBit(7)){ //<=== !bptx quiet 
	BPTX++;
	
	CMStracking->ProceedEvent(dummy_cut, false, false);
	if(CMStracking->NverticesGood()<2){
	  noPU++;
	 
	  ProceedEvent(0, false, false);
	  FillLastEvent(0); // all evts
	  
	  
	  bool t2tr[2] = {false,false}; // minus, plus
	  // <==================================== T2 check
	  if(mc>0){
	   for(short unsigned int side=0; side<2; side++)
	      t2tr[side] = CMSmc->GetT2trigger(bool(side));
	  }else{
	    t2tr[0] = T2->NtracksMinus();
	    t2tr[1] = T2->NtracksPlus();
	  };
	  // <=====================================
	  
	  
	  short unsigned int evtKind = 0;
	  /* -------------> usual loop (before 20th of May)
	  if(t2tr[0]){
	    if(t2tr[1]){
	      ND++; evtKind = 1;
	    }else{
	      SDp++; evtKind = 4;
	    };
	  }else{
	    if(t2tr[1]){
	      SDm++; evtKind = 3;
	    }else{
	      EL++; evtKind = 2;
	    };
	  };
	  */ //<------------------
	  //if(t2tr[0]){ // t2TestPM minus
	  //if(T2->NtracksMinus()>0 && T2->NtracksPlus()==0){ 
	  if(T2->NtracksMinus()==1 && T2->NtracksPlus()==0){ 
	    SDp++; evtKind = 1; // Minus
	    FillLastEvent(evtKind);
	  };
	  //if(t2tr[1]){ // t2TestPM plus
	  //if(T2->NtracksPlus()>0 && T2->NtracksMinus()==0){
	  if(T2->NtracksPlus()==1 && T2->NtracksMinus()==0){
	    SDm++; evtKind = 2; // Plus
	    FillLastEvent(evtKind);
	  };
	  if(t2tr[0]==0 && t2tr[1]==0){
	    EL++; evtKind = 3; // "EL"
	    FillLastEvent(evtKind);
	  } 
	}; //end vtx cut
      }; // end bptx cut
    }; //end trigger cut
  }; // end loop ===============================================================
  std::cout << "total number of events         " << stat   << std::endl;
  std::cout << "worked over                    " << nevts  << std::endl;
  std::cout << "number of tt53 events          " << tt53   << std::endl;
  std::cout << "number of tt52 events          " << tt52   << std::endl;  
  std::cout << "number of al99 events          " << al99   << std::endl;  
  std::cout << "number of ttBOTH events        " << ttBOTH << std::endl;    
  std::cout << "number of BPTX events          " << BPTX   << std::endl;
  std::cout << "no PU events                   " << noPU   << std::endl;
  std::cout << "ND          " << ND   << std::endl;
  std::cout << "SD-         " << SDm  << std::endl;
  std::cout << "SD+         " << SDp  << std::endl;
  std::cout << "EL          " << EL   << std::endl;
  
  return noPU;
}