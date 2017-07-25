/*
 *
 */

#include "uacastor.h"

#include "iostream"

ClassImp(uacastor)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uacastor::uacastor(TChain* tree, 
		     TDirectory * dir,
		     const bool   cmstotem, 
		     const bool   cmsdigis,
		     const short          int MC, 
		     const short unsigned int Ncuts):
      uabase(cmstotem, MC, Ncuts, dir), digi(cmsdigis)
{
  Castor = 0;
  if(tree_combined_flag){
      tree->SetBranchAddress("cmsCastorRecHitsUA",  &Castor);
      if(digi){
      };
  }else{
      tree->SetBranchAddress("castorRecHits",   &Castor);
      if(digi){
      };
  };

  // threshold values
  // igor : ch_thr - bare fC; tower_thr - abs intercalibrated val
  // only individual threshold are from Igor
  for(unsigned int s = 0; s<castor::CSectors; s++){
    for(unsigned int m = 0; m<castor::CModules; m++)
      // convert channel_threshold to GeV
      channel_threshold[s][m] = (castor::noise_channel_gaumean[s][m]+3*castor::noise_channel_gausigma[s][m])*castor::channelGainQE[s][m]*castor::absCasEscaleFactor;
  };

  create_histos();
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uacastor::~uacastor()
{
  std::cout << "uacastor::~uacastor(): deleting " 
	    << h1D->size() << "+" << h2D->size() << " histos" << std::endl;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uacastor::PrintEventInfo(const bool detailed)
{
  std::cout << "uacastor::PrintEventInfo:\n";
  std::cout << "\tcas: E " << total_energy << "\tEem "<< total_em_energy << "\tNch " << n_ch_fired << "\tNtow " << n_towers << "\tNtow5 " << n_towers_5 << "\tNemSec "<< n_towers_em << std::endl;
  std::cout << "\tenergy per sector: \n";
  for(unsigned int i=0; i<castor::CSectors; i++)
    std::cout << "\t"<< 360./16.*i + (360./16./2.) << ": " << tower_energy[i] << " (" << em_energy[i] << ")";
  std::cout << std::endl;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uacastor::FillLastEvent(const short unsigned int cut)
{
  if( cut>=n_cuts ){
    std::cout << "uacastor::FillLastEvent: required cut number is larger that possible, do nothing. Please define larger uaforward::n_cut!\n"; 
    return false;
  }; 
  castor_etotal_h[cut]->Fill(total_energy);
  if(total_energy>0)
    castor_emfract_h[cut]->Fill(total_em_energy/total_energy);
  castor_towers_h[cut]->Fill(n_towers);
  castor_towers5_h[cut]->Fill(n_towers_5);
  castor_towersEM_h[cut]->Fill(n_towers_em);
  return true;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uacastor::ProceedEvent(const short unsigned int cut, const bool fill, const bool info)
{
  if(fill && (cut>n_cuts) ){
    std::cout << "uacastor::ProceedEvent: required cut number is larger that possible, please define larger uaforward::n_cut!\n"; 
    return false;
  }; 

  total_energy = 0;
  total_em_energy    = 0;
  n_ch_fired   = 0;
  n_towers     = 0;
  n_towers_5   = 0;
  n_towers_em  = 0;

  memset(em_energy,     0, sizeof(em_energy));
  memset(tower_energy,  0, sizeof(tower_energy));
  memset(tower_energy5, 0, sizeof(tower_energy5));

  memset(channel_above_threshold,   false, sizeof(channel_above_threshold));
  memset(tower_above_threshold,     false, sizeof(tower_above_threshold));
  memset(tower_above_threshold_5,   false, sizeof(tower_above_threshold_5));  
  memset(tower_above_threshold_em,  false, sizeof(tower_above_threshold_em));  
  memset(towers_modules_above_threshold, 0, sizeof(towers_modules_above_threshold));

    unsigned int ch = 0;
    for(std::vector<MyCastorRecHit>::iterator it=Castor->begin(); it!=Castor->end(); ++it){
      unsigned int m = (*it).mod-1;
      unsigned int s = (*it).sec-1;
      if(m<0 || m>=castor::CModules || s<0 || s>=castor::CSectors)
	break;

      if(castor::channelQuality[s][m]){
	double energy = (*it).energy;
	if(mc<=0) //data 
	  energy*= castor::channelGainQE[s][m]*castor::absCasEscaleFactor; // GeV

	total_energy   +=energy;
	tower_energy[s]+=energy;
	if(m<2){
	  em_energy[s]+=energy;
	  total_em_energy+=energy;
	};
	if(m<5){
	  tower_energy5[s]+=energy;
	}
/*	
	// igor : ch_thr - bare fC; tower_thr - abs intercalibrated val
	// here channel_threshold are already recalculated to GeV (see constructor)
	if(energy>channel_threshold[s][m]){
	  channel_above_threshold[s][m] = 1;
	  n_ch_fired++;
	};
*/
      }; // end channelQuality check
      ch++;
    }; // end loop reco
  if(digi){      
  };

  if(ch!=castor::CChannels){
    std::cout << "uacastor::ProceedEvent: wrong number of good channels!\n"; 
    return 0;
  };

  for(unsigned int s = 0; s<castor::CSectors; s++){

    if(tower_energy[s]>castor::tower_threshold[s]){
      tower_above_threshold[s] = 1;
      n_towers++;
    };
    if(tower_energy5[s]>castor::tower_threshold_5[s]){
      tower_above_threshold_5[s] = 1;
      n_towers_5++;
    };

    // here we need to check HAD later
    if(em_energy[s]>castor::tower_threshold_em[s]){
      tower_above_threshold_em[s] = 1;
    };

    // loop over modules
    bool hadveto = true;
    for(unsigned int m = 0; m<castor::CModules; m++){
      if(m<5){
	if(channel_above_threshold[s][m]){
	  towers_modules_above_threshold[s]++;
	}
      }; // end m<5
      if(m>2 && m<12)
	hadveto = (hadveto && (!channel_above_threshold[s][m]));
    }; // end module loop

    // finish em
    tower_above_threshold_em[s] = (tower_above_threshold_em && hadveto);
    if(tower_above_threshold_em[s])
      n_towers_em++;      
  };// end sector loop

  if(info) PrintEventInfo(true);
  if(fill) FillLastEvent(cut);
  return true;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uacastor::create_histos()
{
  TString title1, title2;
  n_each_h1D = n_cuts;
  castor_etotal_h    = new TH1F*[n_each_h1D];
  castor_emfract_h   = new TH1F*[n_each_h1D];
  castor_towers_h    = new TH1F*[n_each_h1D];
  castor_towers5_h   = new TH1F*[n_each_h1D];
  castor_towersEM_h  = new TH1F*[n_each_h1D];

  for(unsigned int i=0; i<n_each_h1D; i++){
    // total energy [GeV]
    title1 = "castor_etotal_h["; title1+=i; title1+="]";
    title2 = "castor_etotal_h["; title2+=i; title2+="]; E_{tot} [GeV]";
    castor_etotal_h[i] = new TH1F(title1.Data(), title2.Data(), 6500, -500,6000);
    castor_etotal_h[i]->SetDirectory(directory);

    // em fraction
    title1 = "castor_emfract_h["; title1+=i; title1+="]";
    title2 = "castor_emfract_h["; title2+=i; title2+="]; E_{EM}/E_{tot}";
    castor_emfract_h[i] = new TH1F(title1.Data(), title2.Data(), 200, -1,2);
    castor_emfract_h[i]->SetDirectory(directory);

    // n towers above threshold - totel, 5 modules and EM only
    title1 = "castor_towers_h["; title1+=i; title1+="]";
    title2 = "castor_towers_h["; title2+=i; title2+="] total towers; N_{towers}";
    castor_towers_h[i] = new TH1F(title1.Data(), title2.Data(), 18, -1,17);
    castor_towers_h[i]->SetDirectory(directory);

    // n towers above threshold - totel, 5 modules and EM only
    title1 = "castor_towers5_h["; title1+=i; title1+="]";
    title2 = "castor_towers5_h["; title2+=i; title2+="] 5-mod towers; N_{towers}";
    castor_towers5_h[i] = new TH1F(title1.Data(), title2.Data(), 18, -1,17);
    castor_towers5_h[i]->SetDirectory(directory);

    // n towers above threshold - totel, 5 modules and EM only
    title1 = "castor_towersEM_h["; title1+=i; title1+="]";
    title2 = "castor_towersEM_h["; title2+=i; title2+="] EM towers; N_{towers}";
    castor_towersEM_h[i] = new TH1F(title1.Data(), title2.Data(), 18, -1,17);    
    castor_towersEM_h[i]->SetDirectory(directory); 
  };
  h1D->push_back(castor_etotal_h);
  h1D->push_back(castor_emfract_h);
  h1D->push_back(castor_towers_h);
  h1D->push_back(castor_towers5_h);
  h1D->push_back(castor_towersEM_h);
}
