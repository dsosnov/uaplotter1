#ifndef UACALO_H
#define UACALO_H

#include "uabasecentral.h"

#include "MyHFRecHit.h"
#include "MyCaloTower.h"

#include "vector"

class uacalo :  public uabasecentral
{
public:
    uacalo(TChain     * tree, 	 //!<tree of ua format
	   TDirectory * dir,    //!<directory in the output root file
	   const bool cmstotem, //!<true for merged, false for CMS only
	   const short          int MC,  //!< -1, 0 - data; >0 MC
	   const short unsigned int Ncuts //!< number of cuts
	  );
   ~uacalo();
  bool ProceedEvent(const short unsigned int cut=0, const bool fill=false, const bool info=false);
  bool FillLastEvent(const short unsigned int cut);
  void PrintEventInfo(const bool detailed = false);
  bool GetHFtowerTrig(bool plus){return hf_trigger_tower[int(plus)];};
  double GetHFmaxTower(short unsigned int etaBin){
    bool plus = (etaBin >10);
    short unsigned int indx = 10;
    if(plus){
      indx = 22 - etaBin;
    }else{
      indx = etaBin - 3;
    };
    return hf_max_energy_tower[int(plus)][indx];
  };
private:
  std::vector<MyHFRecHit>*     HF;     //!< not implemented yet
  std::vector<MyCaloTower>*    Towers;

  double hf_total_energy_tower[2];  
  double hf_max_energy_tower[2][4];  
  bool   hf_trigger_tower[2];

  unsigned int calotowers[N_ETA_BINS];

  void create_histos();

  TH1F** hf_etotal_minus_h;
  TH1F** hf_etotal_plus_h;

  TH1F** hf_max_towerE_minus_h;
  TH1F** hf_max_towerE_plus_h;

  TH2F** hf_towers_vs_rechits_h;
  TH2F** calotower_e_eta_h;

  double EHF[2];

  ClassDef(uacalo,2);
};

#endif // UACALO_H
