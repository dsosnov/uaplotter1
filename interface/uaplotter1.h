#ifndef UAPLOTTER1_H
#define UAPLOTTER1_H
#include "uabase.h"

#include "uamc.h"
#include "uacmsevtinfo.h"
#include "uatracking.h"
#include "uapf.h"
#include "uacalo.h"
#include "uacastor.h"
#include "uaforward.h"
#include "uat2.h"
#include "uarp.h"

#include "TChain.h"
#include "TFile.h"
#include "TString.h"

class uaplotter1: public uabase {
public:
  uaplotter1(const bool               cmstotem,             //!< — see uabase::uabase()
             const bool               cmsdigis,             //!< → #tree_digi_flag
             const bool               CASTORp,              //!< → #ppb
             const short int          MC,                   //!< — see uabase::uabase()
             const short unsigned int Ncuts = 2             //!< — see uabase::uabase()
            );
  ~uaplotter1();
  int Loop(const int evts, const int trigger);
  int noiseLoop(const int evts, bool bptxQuiet = true);
  int t2Loop(const int evts);

private:
  const bool tree_digi_flag;                                //!< true - with cms digis
  const bool ppb;                                           //!< true for pPb (p->CASTOR), false for Pbp
  bool zdc56;                                               //!< flag for ZDCdigi proper TS, initialized in initializeChain
  const short unsigned int dummy_cut;                       //!< service variable for processing without histo filling
  unsigned int current_event;

  TChain     *chainTree;
  TFile      *outputFile;
  TString     initializeChain();

  uamc          *CMSmc;
  uacmsevtinfo  *CMSevtinfo;
  uatracking    *CMStracking;
  uapf          *CMSpf;
  uacalo        *CMScalo;
  uacastor      *CMScastor;
  uaforward     *CMSforward;
  uat2          *T2;
  uarp          *RP;

  int DefineTrigger(int trigger_bit_common, bool &tech_bit);
  bool ProceedTrigger(int trigger_bit, bool tech_bit);

  bool ProceedEvent(const short unsigned int cut = 0, const bool fill = false, const bool info = false);
  bool FillLastEvent(const short unsigned int cut);
  void PrintEventInfo(const bool detailed = false);

  void IniRapGapRange(); //!< to be called during construction. Sets value to first_central_bin and last_central_bin according to uathresholds
  bool FindRapGap(bool RECO = true); //!< if true - works on RECO, false - on MCtruth
  bool TotalRapGap(short unsigned int ind, bool active_minus, bool active_plus);
  void PrintRapGap();
  short unsigned int first_central_bin;
  short unsigned int last_central_bin;
  bool combined_central_activity[N_ETA_BINS];

  // [0]reco, [1]mctruth eta_binning, [2] mctruch total
  short int sd_flag_central[2]; //!<central_sd_flag: 0 - ND, -1 - SD-, +1 - SD+, 2 DD~central gap(s), 3  cd cand, 4 "elastic"
  short int sd_flag_total[3];   //!< all MCtruth for MC, T2 events for RECO
  short unsigned int n_sd_minus_bins[2];
  short unsigned int n_sd_plus_bins[2];
  short unsigned int n_dd_rg_bins[2];

  void CalculateSDdiffMass(bool info = false); //!< for SD only!!
  void PrintSDdiffMass(bool detailed);

  // [0]reco, [1]mctruth
  double xi_pf[2];
  double xi_calo[2] ;
  double xi_cas[2]  ;
  double xi_zdc[2]  ;
  double xi_full[2] ;
  double xi_mc_out;   // this includes ZDC and CASTOR also!
  double xi_mc_total; // everything from first active bin in reco

  void create_histos();
  TH2F **diff_flag_mc_full_reco_central_h;  //!< X:mc_total; Y:reco_central
  TH2F **diff_flag_mc_full_reco_full_h;    //!< X:mc_total; Y:reco_total
  TH2F **diff_flag_mc_full_mc_central_h;    //!< X:mc_full; Y:mc_central
  TH2F **diff_flag_mc_total_mc_central_h;    //!< X:mc_total; Y:mc_central
  TH2F **n_sd_minus_bins_mc_reco_h;          //!< X:mc; Y:reco
  TH2F **n_sd_plus_bins_mc_reco_h;          //!< X:mc; Y:reco

  TH2F **xi_mc_p_mc_total_h;       //!< X xi(p), Y:xi_mc_total
  TH2F **xi_mc_p_reco_full_h;       //!< X xi(p), Y:xi_full reco (but does not care about diffraction); MC only
  TH2F **xi_p_reco_full_h;          //!< X:RP xi or MC proton; Y:reco_total xi (combines the above and data case with RP)

  TH2F **xi_mc_total_mc_full_h;
  TH2F **xi_mc_total_reco_full_h;

  TH2F **xi_calo_mc_reco_h;
  TH2F **xi_pf_mc_reco_h;
  TH2F **xi_cas_mc_reco_h;
  TH2F **xi_zdc_mc_reco_h;

  TH1F **n_sd_minus_bins_h;
  TH1F **n_sd_plus_bins_h;
  TH1F **central_activity_h;
  TH1F **central_activity_mc_h;

  TH1F **xi_reco_full_h;

  TH2F **zdcM_vs_castor_h;
  TH2F **zdcM_vs_T2primM_h;
  TH2F **ZDCm_vs_xiRP_h;
  TH2F **ZDCp_vs_xiRP_h;
  TH2F **FSCmSi8_vs_xiRP_h;
  TH2F **FSCmN_vs_xiRP_h;
  TH2F **FSCmN_vs_castor_h;

  ClassDef(uaplotter1, 2);
};

#endif // UAPLOTTER1_H
