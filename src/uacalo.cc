/*
 *
 */

#include "uacalo.h"

#include <map>
#include "iostream"

ClassImp(uacalo)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uacalo::uacalo(TChain* tree, 
	       TDirectory* dir, 
	       const bool cmstotem, 
	       const short int MC, 
	       const short unsigned int Ncuts): uabasecentral(cmstotem, MC, Ncuts, dir)
{
  HF     = 0;
  Towers = 0;
  if(tree_combined_flag){
    tree->SetBranchAddress("cmshfRecHitsUA",   &HF);
    tree->SetBranchAddress("cmsCaloTowersUA",  &Towers);
  }else{
    tree->SetBranchAddress("hfRecHits",   &HF);
    tree->SetBranchAddress("caloTowers",  &Towers);
  };

  create_histos();
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uacalo::~uacalo()
{
  std::cout << "uacalo::~uacalo(): deleting " 
	    << h1D->size() << "+" << h2D->size() << " histos" << std::endl;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uacalo::FillLastEvent(const short unsigned int cut)
{
  if( cut>=n_cuts ){
    std::cout << "uatracking::FillLastEvent: required cut number is larger that possible, do nothing. Please define larger uaforward::n_cut!\n"; 
    return false;
  }; 
  hf_etotal_minus_h[cut]->Fill(hf_total_energy_tower[0]);
  hf_etotal_plus_h [cut]->Fill(hf_total_energy_tower[1]);
  for(unsigned int bin = 0; bin<N_ETA_BINS; bin++)
    if( (bin>=BINMIN) && (bin<=BINMAX) )
      calotower_e_eta_h[cut]->Fill(ETA_BIN_L[bin], energy[bin]);
  
  double MaxHFtow[2] = {0,0};
  for(short unsigned int side=0; side<2; side++)
    for(short unsigned int indx=0; indx<4; indx++){
      if(MaxHFtow[side]<hf_max_energy_tower[side][indx])
	MaxHFtow[side]=hf_max_energy_tower[side][indx];
    };
  hf_max_towerE_minus_h[cut]->Fill(MaxHFtow[0]);
  hf_max_towerE_plus_h[cut]->Fill(MaxHFtow[1]);
  return true;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uacalo::PrintEventInfo(const bool detailed)
{  
  std::cout << "uacalo::PrintEventInfo:\n\t";
  std::cout << "tower trigger HF- " << hf_trigger_tower[0] 
	    << "\t; HF+ "             << hf_trigger_tower[1] << std::endl;
  std::cout << "\tenergy HF-         " << hf_total_energy_tower[0] 
	    << "\t; HF+ "            << hf_total_energy_tower[1] << std::endl;
  PrintActivity(false);    
  if(detailed){
    std::cout << "\tnumber of calotowers per bin:\n\t";
    for(unsigned short int bin=0; bin<N_ETA_BINS; bin++){
      std::cout << calotowers[bin] << "  ";
    };std::cout << std::endl;
    std::cout << "\tenergy per bin:\n\t";
    for(unsigned short int bin=0; bin<N_ETA_BINS; bin++){
      std::cout << energy[bin] << "\t";
    };std::cout << std::endl;
  };
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool uacalo::ProceedEvent(const short unsigned int cut, const bool fill, const bool info)
{
  
  memset(hf_total_energy_tower,  0, sizeof(hf_total_energy_tower));
  memset(hf_max_energy_tower, 0, sizeof(hf_max_energy_tower));
  memset(hf_trigger_tower,   false, sizeof(hf_trigger_tower));
  memset(calotowers,             0, sizeof(calotowers));
  PrepareArrays(); // cleans eta arrays
  
  for(std::vector<MyCaloTower>::iterator it=Towers->begin(); it!=Towers->end(); ++it){
      if( (*it).Pt()>0 ){
	int bin = find_eta_bin((*it).eta());
	if(bin>=0){
	  energy[bin] += (*it).energy();
	  pz[bin]     += (*it).Pz();
	  pt[bin]     += (*it).Pt();
	  if((*it).hasHF){
	    if((*it).energy()>HF_TOWER_THR)
		calotowers[bin]++;
	  }else{
	    if((*it).energy()>CENT_TOWER_THR)
		calotowers[bin]++;   
	  };
	}; // end if bin
	
	//___________________________________________________________________
	if(((*it).hasHF) && (fabs((*it).eta())<=CALO_ETA_ACC) ){
	  bool minus = (bin<10);
	  short unsigned int indx = 10;
	  if(minus){
	    indx = bin - 3;
	  }else{
	    indx = 22 - bin;
	  };
	  short unsigned int side = int( not(minus)); // 0 - minus, 1 - plus
	  //std::cout << (*it).PseudoRapidity() << "\t" << bin << "\t" << indx << std::endl;
	  if( (*it).energy()>hf_max_energy_tower[side][indx]){
	    hf_max_energy_tower[side][indx]=(*it).energy();
	    //std::cout << (*it).PseudoRapidity() << "\t" <<bin << "\t" << side << " " << indx << "\t" << hf_max_energy_tower[side][indx] << std::endl;
	  };
	  // here HF "trigger" only 
	  if(fabs((*it).eta())<HF_ETA_MAX){
	    short unsigned ind = int((*it).zside>0);
	    hf_total_energy_tower[ind]+=(*it).energy();
	    if(!hf_trigger_tower[ind] && ((*it).energy()>HF_TOWER_THR))
	      hf_trigger_tower[ind] = true;
	  };
	}; // end HF
	//___________________________________________________________________
      };    
  };// end tower loop
  
  // <===================================================compare to threshold
  for(unsigned int bin = 0; bin<N_ETA_BINS; bin++)
    if( (THR_CALO[bin]>0) && (energy[bin]>THR_CALO[bin]) ){
      activity_loose[bin]=true;
      activity_tight[bin]=true;
    };
  
  if(info)
    PrintEventInfo(true);
  if(fill)  
    FillLastEvent(cut);

  return true;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uacalo::create_histos()
{
  TString title1, title2;
  
  n_each_h1D = n_cuts;
  n_each_h2D = n_cuts;

  hf_etotal_minus_h = new TH1F * [n_each_h1D];
  hf_etotal_plus_h  = new TH1F * [n_each_h1D];
  hf_max_towerE_minus_h = new TH1F * [n_each_h1D];
  hf_max_towerE_plus_h  = new TH1F * [n_each_h1D];
  //hf_towers_vs_rechits_h = new TH2F * [n_each_h2D];
  calotower_e_eta_h        = new TH2F * [n_each_h2D];

  // rechits are not done yet
  for(unsigned int i=0; i<n_cuts; i++){
    title1 = "hf_etotal_minus_h["; title1+=i; title1+="]";
    title2 = title1; title2+="; E_{HFtowers} [GeV]";
    hf_etotal_minus_h[i] = new TH1F(title1.Data(), title2.Data(), 7100, -100, 7000);
    hf_etotal_minus_h[i]->SetDirectory(directory);
    
    title1 = "hf_max_towerE_minus_h["; title1+=i; title1+="]";
    title2 = title1; title2+="; E^{max}_{HFtower} [GeV]";
    hf_max_towerE_minus_h[i] = new TH1F(title1.Data(), title2.Data(), 7100, -100, 7000);
    hf_max_towerE_minus_h[i]->SetDirectory(directory);
    
    title1 = "hf_etotal_plus_h["; title1+=i; title1+="]";
    title2 = title1; title2+="; E_{HFtowers} [GeV]";
    hf_etotal_plus_h[i] = new TH1F(title1.Data(), title2.Data(), 7100, -100, 7000);
    hf_etotal_plus_h[i]->SetDirectory(directory);
    
    title1 = "hf_max_towerE_plus_h["; title1+=i; title1+="]";
    title2 = title1; title2+="; E^{max}_{HFtower} [GeV]";
    hf_max_towerE_plus_h[i] = new TH1F(title1.Data(), title2.Data(), 7100, -100, 7000);
    hf_max_towerE_plus_h[i]->SetDirectory(directory);
    
//     title1 = "hf_towers_vs_rechits_h["; title1+=i; title1+="]";
//     title2 = title1; title2+="; E_{HFrechits} [GeV]; E_{HFtowers}";
//     hf_towers_vs_rechits_h[i] = new TH2F(title1.Data(), title2.Data(), 2100, -100, 2000, 2100, -100, 2000);    
//     hf_towers_vs_rechits_h[i]->SetDirectory(directory);
    
    title1 = "calotower_e_eta_h["; title1+=i; title1+="]";
    title2 = title1; title2+="; #eta; E_{towers}";
    calotower_e_eta_h[i] = new TH2F(title1.Data(), title2.Data(), 28, -7,7, 21000, -100, 2000);    
    calotower_e_eta_h[i]->SetDirectory(directory);
  };
  h1D->push_back(hf_etotal_minus_h);
  h1D->push_back(hf_etotal_plus_h);
  h1D->push_back(hf_max_towerE_minus_h);
  h1D->push_back(hf_max_towerE_plus_h);
  //h2D->push_back(hf_towers_vs_rechits_h);
  h2D->push_back(calotower_e_eta_h);
}
