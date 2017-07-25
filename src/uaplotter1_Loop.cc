
#include "uaplotter1.h"

#include "TString.h"
#include "TMath.h"

#include "iostream"
#include "stdlib.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int uaplotter1::Loop(const int evts, const int trigger)
{
  if (mc < 0)
    return noiseLoop(evts);

  bool tech_bit   = true;
  int trigger_bit = DefineTrigger(trigger, tech_bit); // just translate initial trigger value

  int stat = chainTree->GetEntries();
  unsigned int nevts = evts;
  if (evts == -1 || evts > stat) nevts = stat;
  std::cout << "Total stat = " << stat << std::endl;
  std::cout << "uaplotter1::Loop(" << trigger_bit << ") for " << nevts << " events\n";

  unsigned int selected_evts[n_cuts];
  memset(selected_evts, 0, sizeof(selected_evts));

  unsigned int trigger_evts = 0;
  unsigned int t2prim_evts  = 0;
  unsigned int bptx_active  = 0;
  unsigned int vtx_cut      = 0;
  unsigned int sd_minus     = 0;
  unsigned int sd_plus      = 0;
  unsigned int sd_minus_eta[11];
  unsigned int sd_plus_eta[11];
  unsigned int goodFSC      = 0;
  memset(sd_minus_eta, 0, sizeof(sd_minus_eta));
  memset(sd_plus_eta,  0, sizeof(sd_plus_eta));

  unsigned int kevt = 0;
  for (long unsigned int i = 0; i < nevts; i++) {
    unsigned int kevt_current = i / 1000;
    if (kevt_current > kevt) {
      kevt = kevt_current;
      std::cout << kevt << std::endl;
    };
    chainTree->GetEntry(i);
    current_event = i;

    memset(sd_flag_central, 0, sizeof(sd_flag_central));
    memset(sd_flag_total,   0, sizeof(sd_flag_total));

    if (mc > 0) { // <============================  do MC loop here
      CMSmc->ProceedEvent(dummy_cut, false, false);
    };


    if (ProceedTrigger(trigger_bit, tech_bit)) {
      trigger_evts++;
      
      if( !CMSevtinfo->GetL1Bit(9)){ //<=== !bptx quiet 
	bptx_active++;
	CMStracking->ProceedEvent(dummy_cut, false, false);
	if(CMStracking->NverticesGood()<2){
	  vtx_cut++;
	 
	  // <==================================== T2 check
	  bool t2prim = false;
	  if(mc>0){
	    if( tech_bit && (trigger_bit==53) ){
	      t2prim = true; // we already checked in the ProceedTrigger; for MC it is the same!
	    }else{
	      t2prim = ProceedTrigger(53, true); 
	    };
	  }else if(tree_combined_flag && (mc==0)){
	    T2->ProceedEvent(dummy_cut, false, false);
	    t2prim = ( (T2->NPrimtracksMinus()>0) || (T2->NPrimtracksPlus()>0) );
	  };
	  if(t2prim) t2prim_evts++;
	  // <=====================================
	  
	  
	  ProceedEvent(0, false, false);
	  if(tree_digi_flag && CMSforward->FSCvalid())
	    goodFSC++;
	  
	  FillLastEvent(0); // -------------------------------->  all triggered events
	    
	  if(ppb){
	    if(RP->Valid()){
	      PrintEventInfo(true);
	      T2->PrintAllTracks();
	      //RP->PrintEventInfo(true);
	      //CMSforward->PrintEventInfo(true);
	    };
	  };

	  
	  bool proper_proton = false;
	  if(mc>0){
	    proper_proton = ( (CMSmc->IntactProton()!=0) && ( (CMSmc->IntactProton()==-1) == ppb ) );
	    if(proper_proton){   
	      FillLastEvent(1); // ---------------------------> all events with proton in proper direction
	      if(CMSmc->IntactProtonE()>3950){
		FillLastEvent(2); //--------------------------> energetic protons
		PrintEventInfo();
	      };
	    }; // end proper_proton
	  }else if (mc==0){
	    // Near and Far are of the same sign
	    proper_proton = ( (RP->Valid()) && ( (RP->trackValidUp()) || (RP->trackValidDn()) ) );
	    if(proper_proton){
	      FillLastEvent(2);
	    };
	  };
	  
	  
	  if(sd_flag_total[0]==4){ // elastic candidate
	    CalculateSDdiffMass(false);
	    FillLastEvent(3); // -----------------------------> elastic
	  };
	  
	  if(sd_flag_total[0]==0){ // ND
	      FillLastEvent(17); //---------------------------> ND
	  };
	  
	  if(sd_flag_total[0]==-1){
	    CalculateSDdiffMass(false);
	    FillLastEvent(4); // -----------------------------> all SD- events "4"
	    
	    /*
	    if(proper_proton){
	      PrintEventInfo(true);
	      //CMSmc->ProceedEvent(dummy_cut, false, true);
	    };*/
	    sd_minus++;
	    for(short unsigned int ii=0; ii<11; ii++){ // RG: [0,1)="5", [1,2)="6", ..., [10,11)="15"
	      short unsigned int cut   = 5+ii;
	      short unsigned int rgbin = 2*ii;     
	      if( (n_sd_minus_bins[0]>=rgbin) && (n_sd_minus_bins[0]<(rgbin+2)) ){
		FillLastEvent(cut);
		sd_minus_eta[ii]++;
	      };
	    };
	  }else if (sd_flag_total[0]==1){//-------------------> all SD+ events "19"
	    CalculateSDdiffMass(false);
	    FillLastEvent(19);
	    
	    /*
	    if(proper_proton){
	      PrintEventInfo(true);
	      CalculateSDdiffMass(true);
	      //CMSmc->ProceedEvent(dummy_cut, false, true);
	    };*/
	    sd_plus++;
	    for(short unsigned int ii=0; ii<11; ii++){// RG: [0,1)="20", [1,2)="21", ..., [10,11)="30"
	      short unsigned int cut   = 20+ii;
	      short unsigned int rgbin = 2*ii;
	      if( (n_sd_plus_bins[0]>=rgbin) && (n_sd_plus_bins[0]<(rgbin+2)) ){
		FillLastEvent(cut);
		sd_plus_eta[ii]++;
	      };
	    };
	  }; // end SD cases
	  
	}; // end vertices<2
      }; // end !bptx quiet
    }; // end trigger

  };// end loop
  std::cout << "Acceptance: [" << ETA_BIN_L[first_central_bin] << "," << ETA_BIN_L[last_central_bin] + ETA_BIN_W << "]\n";
  std::cout << "Total evts in chain       : " << stat << std::endl;
  std::cout << "Proceeded evts            : " << (current_event + 1)     << std::endl;
  std::cout << "Triggered evts            : " << trigger_evts << std::endl;
  std::cout << "T2 active (not selecting) : " << t2prim_evts  << std::endl;
  std::cout << "Active (!bptx quiet)      : " << bptx_active << std::endl;
  std::cout << "No PU vertices            : " << vtx_cut << std::endl;
  std::cout << "Good FSC evt              : " << goodFSC << std::endl;
  std::cout << "\tsd- candidates          : " << sd_minus << "\t\tsd+ candidates      : " << sd_plus << std::endl;
  for (short unsigned int ii = 0; ii < 11; ii++) {
    std::cout << ii << "<= |deta| < " << ii + 1 << "\t\t" << sd_minus_eta[ii] << "\t\t" << sd_plus_eta[ii] << std::endl;
  };
  return nevts;
}
