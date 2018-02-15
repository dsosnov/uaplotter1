#ifndef UACALO_H
#define UACALO_H

#include "uabasecentral.h"

#include "MyHFRecHit.h"
#include "MyCaloTower.h"

#include "vector"

class uacalo :  public uabasecentral {
public:
  uacalo(TChain*              tree,     //!<tree of ua format
         TDirectory*          dir,      //!<directory in the output root file
         const bool           cmstotem, //!<true for merged, false for CMS only
         const short          MC,       //!< -1, 0 - data; >0 MC
         const unsigned short Ncuts     //!< number of cuts
  );
  ~uacalo();
  bool ProceedEvent(const short unsigned int cut = 0, const bool fill = false, const bool info = false);
  bool FillLastEvent(const short unsigned int cut);
  void PrintEventInfo(const bool detailed = false);
  bool GetHFtowerTrig(bool plus) {return hf_trigger_tower[int(plus)];};
  double GetHFmaxTower(unsigned short etaBin)
  {
    double eta = find_eta(etaBin);
    short unsigned int indx = find_eta_bin(CALO_ETA_ACC) - find_eta_bin(fabs(eta));
    short unsigned int side = int(eta>0);
    return hf_max_energy_tower[side][indx];
  }
  double GetHFmax(unsigned short side){
    double MaxHFtow = 0;
    for (short unsigned int indx = 0; indx < nBinsHF; indx++) {
      if (MaxHFtow < hf_max_energy_tower[side][indx])
        MaxHFtow = hf_max_energy_tower[side][indx];
    }
    return MaxHFtow;
  }
private:
  std::vector<MyHFRecHit>     *HF;     //!< not implemented yet
  std::vector<MyCaloTower>    *Towers;

  double hf_total_energy_tower[2];
  uint nBinsHF = find_eta_bin(CALO_ETA_ACC) - find_eta_bin(HF_ETA_MIN) + 1;
  double* hf_max_energy_tower[2];
  double hf_max_energy_tower_ill[2];
  bool   hf_trigger_tower[2];

  unsigned int calotowers[N_ETA_BINS];
  double energyMax[N_ETA_BINS];

  void create_histos();

  TH1F **hf_etotal_minus_h;
  TH1F **hf_etotal_plus_h;

  TH1F **hf_max_towerE_minus_h;
  TH1F **hf_max_towerE_plus_h;
  TH1F **hf_max_towerE_minus_ill_h;                       // Ill at tower \eta=-37 [-4.363, -4.191], \phi=9 [-1.7453, -1.5708] [-M_PI+M_PI/18*8, -M_PI+M_PI/18*9]
  TH1F **hf_max_towerE_plus_ill_h;                        // Ill at tower \eta=35  [3.839, 4.013],   \phi=6 [-2.2689, -2.0944] [-M_PI+M_PI/18*5, -M_PI+M_PI/18*6]

  TH2F **hf_towers_vs_rechits_h;
  TH2F **calotower_e_eta_h;
  TH2F **calotower_eMaxTower_eta_h;

  TH2F **calotower_e_eta_phi_h;
  TH2F **calotower_e_eta_phi_minus_h;
  TH2F **calotower_e_eta_phi_plus_h;
  TH2F **hf_maxTowerE_minus_plus_h;

  double EHF[2];

  ClassDef(uacalo, 2);
};

#endif // UACALO_H
