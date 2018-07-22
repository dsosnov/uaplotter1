#ifndef __UATHRESHOLDS_H__
#define __UATHRESHOLDS_H__

namespace uathresholds {
  constexpr short int RNDT      = -100;

  //this is wrong for MC
  //const double ENERGYN      = 3999.98095727174496; // energy in GeV from sqrt(S_NN) = 5.023
  constexpr double MOMBEAM        = 3997.7190062; // From EPOS, for Run 1, 5 TeV

  constexpr float HF_ETA_MAX    = 5.20; // for HF "trigger" only
  constexpr float HF_ETA_MIN    = 3.15; // for HF "trigger" only

  constexpr float  CAS_ETA_MIN  = -6.6;
  constexpr float  CAS_ETA_MAX  = -5.2;
  constexpr float  CAS_ETA_AVE  = -5.9;
  constexpr double CAS_TANH_AVE = -9.99984990996805823e-01;

  constexpr float T2_ABSETA_MIN = 5.3;
  constexpr float T2_ABSETA_MAX = 6.5;
  constexpr float T2_PT_THR     = 0.04;

  constexpr float CENT_TOWER_THR  = 2.0;        // energy in GeV
  constexpr float HF_TOWER_THR[2] = {8.0, 8.0}; // energy in GeV
  constexpr float ZDC_PXY_THR     = 0.495;      // hernja!!!
  constexpr float MIN_ZDC_ETA     = 8.;         //8.2?

  constexpr float MC_E_THR       = 0.;          // energy in GeV

  constexpr float CENT_ETA_ACC   = 2.0;         // eta
  //const float CENT_ETA_ACC   = 3.0; // 16Feb15 - test

  constexpr float CALO_ETA_ACC   = HF_ETA_MAX; // eta *was 5

  constexpr float RG_ETA_ACC   = 3.0; // Maximum eta for getting RG in CMS. Default: CALO_ETA_ACC. Change when you want to find rapidity gaps for \eta = (-RG_ETA_ACC, RG_ETA_ACC)

  constexpr float TRCK_ETA_MAX = 2.5; //eta = 2.5
  constexpr float TRCK_ETA_MIN = -2.5;  //eta = -2.5
  constexpr int   TRCK_ETA_MAX_BIN = 17; //eta = 2.5
  constexpr int   TRCK_ETA_MIN_BIN = 8;  //eta = -2.5

  constexpr float TRCK_PT_THR = 0.2; // pT threshold in GeV

  constexpr float ETA_BIN_W = 0.5;
  constexpr float ETA_MAX   = 6.5;//CALO_ETA_ACC;//ETA_BIN_L[N_ETA_BINS-1]+ETA_BIN_W;

  constexpr unsigned int N_ETA_BINS = ETA_MAX * 2 / ETA_BIN_W /*26*/;
  constexpr float ETA_BIN_L[N_ETA_BINS] = {-6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6};

  constexpr int BINMIN = (ETA_MAX - CALO_ETA_ACC) / ETA_BIN_W    ;
  constexpr int BINMAX = (ETA_MAX + CALO_ETA_ACC) / ETA_BIN_W - 1;

  static int _energy = 8;
  static bool _p_to_castor   = true;
  static int _run    = 2;

  // TODO remove this perversions and replace to normal vay to define (by classes, for example)

  // Thresholds for Calo Towers (E_{summary})
  static std::vector<float> _thr_calo_sume;
  inline std::vector<float> const THR_CALO_SUME(){
    std::vector<float> t0        = { 0, 0, 0,    0,    0,    0,    0,    0,   0,    0,    0,   0,   0,   0,   0,    0,    0,   0,    0,    0,    0,    0,    0, 0, 0, 0 };
    std::vector<float> c_1       = { 0, 0, 0,  2.9,  3.9,  3.4,  3.1, 12.8, 9.3, 11.2, 14.8, 8.4, 8.9, 8.2, 8.0, 13.7, 11.5, 8.8, 11.9,  3.4,  4.0,  4.7,  3.7, 0, 0, 0 }; // Run1, 5 TeV
    std::vector<float> c_2_5_pbp = { 0, 0, 0, 12.8, 23.1, 22.4, 30.0, 22.6, 5.1,  4.8,  4.0, 3.3, 5.0, 4.6, 3.1,  4.1,  5.4, 5.4, 16.6, 43.1, 39.2, 33.3, 19.6, 0, 0, 0 }; // Run2, 5TeV Pbp (Pb-> <-p)
    std::vector<float> c_2_5_ppb = { 0, 0, 0, 12.7, 22.8, 22.2, 29.7, 22.3, 4.9,  4.6,  3.8, 3.2, 5.0, 4.6, 3.1,  4.0,  5.3, 5.2, 16.2, 40.0, 38.8, 32.0, 18.9, 0, 0, 0 }; // Run2, 5TeV pPb (p-> <-Pb)
    std::vector<float> c_2_8_pbp = { 0, 0, 0, 14.5, 24.7, 23.9, 31.1, 22.4, 5.1,  4.8,  4.0, 3.3, 5.0, 4.6, 3.1,  4.1,  5.4, 5.4, 16.7, 43.7, 41.0, 34.9, 21.3, 0, 0, 0 }; // Run2, 8TeV Pbp (Pb-> <-p)
    std::vector<float> c_2_8_ppb = { 0, 0, 0, 16.1, 26.5, 25.3, 32.3, 22.1, 5.0,  4.7,  3.8, 3.3, 5.0, 4.6, 3.0,  4.0,  5.2, 5.3, 16.2, 40.6, 40.3, 33.7, 20.6, 0, 0, 0 }; // Run2, 8TeV pPb (p-> <-Pb)
    auto thr_p = &t0;
    switch (_run){
      case 1: {
        thr_p = &c_1;
        break;
      }
      case 2:{
        switch (_energy){
          case 5: {
            thr_p = (_p_to_castor) ?  &c_2_5_pbp:  &c_2_5_ppb;
            break;
          }
          case 8: {
            thr_p = (_p_to_castor) ?  &c_2_8_pbp:  &c_2_8_ppb;
            break;
          }
        }
      }
    }
    if(_thr_calo_sume != *thr_p) _thr_calo_sume = *thr_p;
    return _thr_calo_sume;
  }

  // Thresholds for Calo Towers (E_{MaxTower})
  static std::vector<float> _thr_calo_twre;
  inline std::vector<float> const THR_CALO_TWRE(){
    std::vector<float> t0        = { 0, 0, 0,   0,   0,   0,   0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,    0,   0,    0,   0,   0, 0, 0, 0 };
    std::vector<float> c_2_5_pbp = { 0, 0, 0, 3.2, 6.2, 3.6, 4.2, 13.4, 2.4, 1.9, 1.6, 1.5, 1.7, 1.5, 1.3, 1.9, 2.1, 2.7, 12.1, 6.6, 12.7, 4.8, 4.5, 0, 0, 0 }; // Run2, 5TeV Pbp (Pb-> <-p)
    std::vector<float> c_2_5_ppb = { 0, 0, 0, 3.1, 6.2, 3.3, 4.1, 13.3, 2.3, 1.9, 1.6, 1.4, 1.7, 1.4, 1.3, 1.9, 2.1, 2.6, 11.8, 6.6, 14.2, 4.8, 4.4, 0, 0, 0 }; // Run2, 5TeV pPb (p-> <-Pb)
    std::vector<float> c_2_8_pbp = { 0, 0, 0, 6.3, 7.1, 6.7, 7.3, 13.4, 2.4, 2.0, 1.6, 1.5, 1.7, 1.5, 1.3, 1.9, 2.1, 2.8, 12.1, 8.0, 14.1, 8.1, 7.5, 0, 0, 0 }; // Run2, 8TeV Pbp (Pb-> <-p)
    std::vector<float> c_2_8_ppb = { 0, 0, 0, 7.2, 7.7, 7.8, 8.9, 13.4, 2.4, 2.0, 1.6, 1.4, 1.7, 1.5, 1.3, 1.8, 2.1, 2.7, 11.9, 8.1, 14.3, 8.2, 7.3, 0, 0, 0 }; // Run2, 8TeV pPb (p-> <-Pb)
    auto thr_p = &t0;
    switch (_run){
//       case 1: { // TODO find thresholds for maximum energy in tower from first run
//         thr_p = &c_1;
//         break;
//       }
      case 2:{
        switch (_energy){
          case 5: {
            thr_p = (_p_to_castor) ?  &c_2_5_pbp:  &c_2_5_ppb;
            break;
          }
          case 8: {
            thr_p = (_p_to_castor) ?  &c_2_8_pbp:  &c_2_8_ppb;
            break;
          }
        }
      }
    }
    if(_thr_calo_twre != *thr_p) _thr_calo_twre = *thr_p;
    return _thr_calo_twre;

  }

  // Thresholds for Particle Flow (E_{summary}): neutral particles
  static std::vector<float> _thr_pfen_sume;
  inline std::vector<float> const THR_PFEN_SUME(){
    std::vector<float> t0         = { 0, 0,    0,    0,    0,    0,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,    0,    0,    0,    0,    0, 0, 0, 0 };
    std::vector<float> pf_1 =       { 0, 0,    0,  2.2,  2.5,  2.7,  2.4, 16.4, 9.5, 8.9, 7.0, 2.9, 4.7, 3.8, 2.6, 5.7, 9.2, 8.9, 15.5,  2.8,  3.0,  3.1,  2.9, 0, 0, 0 }; // Run1, 5 TeV
    std::vector<float> pf_2_5_pbp = { 0, 0,    0,  5.7,  7.1,  5.7,  7.9, 27.4, 5.6, 3.5, 2.4, 1.9, 2.7, 2.3, 1.7, 2.4, 3.8, 5.9, 21.8, 11.2, 27.7, 11.6,  9.7, 0, 0, 0 }; // Run2, 5TeV Pbp (Pb-> <-p)
    std::vector<float> pf_2_5_ppb = { 0, 0,    0,  5.7,  7.2,  5.4,  7.8, 27.0, 5.3, 3.4, 2.4, 1.9, 2.7, 2.3, 1.6, 2.3, 3.7, 5.6, 21.2, 10.8, 30.3, 11.2,  9.4, 0, 0, 0 }; // Run2, 5TeV pPb (p-> <-Pb)
    std::vector<float> pf_2_8_pbp = { 0, 0, 10.3, 11.2, 13.4, 12.3, 13.7, 27.8, 5.6, 3.5, 2.5, 1.9, 2.7, 2.3, 1.7, 2.4, 3.9, 5.9, 22.5, 16.3, 31.2, 17.6, 14.2, 0, 0, 0 };
    std::vector<float> pf_2_8_ppb = { 0, 0, 13.5, 13.7, 16.8, 15.7, 17.5, 27.4, 5.5, 3.6, 2.5, 1.9, 2.7, 2.3, 1.7, 2.4, 3.8, 5.8, 21.9, 16.3, 31.7, 17.5, 14.0, 0, 0, 0};
    std::vector<float> pf_2_8_pbp_unified = { 0, 0, 0, 0, 0, 0, 0, 13.4, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 13.4, 0, 0, 0, 0, 0, 0, 0 }; // Run2, 8TeV Pbp (Pb-> <-p) [new thresholds, unified for central region]
    std::vector<float> pf_2_8_ppb_unified = { 0, 0, 0, 0, 0, 0, 0, 13.4, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 13.4, 0, 0, 0, 0, 0, 0, 0 }; // Run2, 8TeV pPb (p-> <-Pb) [new thresholds, unified for central region]
    auto thr_p = &t0;
    switch (_run){
      case 1: {
        thr_p = &pf_1;
        break;
      }
      case 2:{
        switch (_energy){
          case 5: {
            thr_p = (_p_to_castor) ?  &pf_2_5_pbp:  &pf_2_5_ppb;
            break;
          }
          case 8: {
            thr_p = (_p_to_castor) ?  &pf_2_8_pbp_unified:  &pf_2_8_ppb_unified;
            break;
          }
        }
      }
    }
    if(_thr_pfen_sume != *thr_p) _thr_pfen_sume = *thr_p;
    return _thr_pfen_sume;
  }

  // Thresholds for Particle Flow (E_{MaxTower}), neutral particles
  static std::vector<float> _thr_pfen_twre;
  inline std::vector<float> const THR_PFEN_TWRE(){
    std::vector<float> t0         = { 0, 0, 0,    0,    0,    0,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,    0,    0,    0,    0,    0, 0, 0, 0 };
    std::vector<float> pf_2_5_pbp = { 0, 0, 0,  4.5,  6.3,  5.1,  6.2, 17.2, 3.1, 2.1, 1.7, 1.6, 2.0, 1.7, 1.5, 1.9, 2.2, 3.5, 16.0,  6.6, 25.1,  6.4,  6.0, 0, 0, 0 }; // Run2, 5TeV Pbp (Pb-> <-p)
    std::vector<float> pf_2_5_ppb = { 0, 0, 0,  4.3,  6.6,  4.6,  6.2, 17.1, 3.0, 2.1, 1.7, 1.6, 2.1, 1.7, 1.5, 1.9, 2.2, 3.4, 15.6,  6.4, 28.3,  6.4,  6.0, 0, 0, 0 }; // Run2, 5TeV pPb (p-> <-Pb)
    std::vector<float> pf_2_8_pbp = { 0, 0, 0,  9.9, 10.9, 11.0, 11.9, 17.3, 3.2, 2.1, 1.8, 1.7, 2.1, 1.7, 1.5, 1.9, 2.3, 3.6, 16.3, 12.4, 27.9, 12.8, 11.1, 0, 0, 0 }; // Run2, 8TeV Pbp (Pb-> <-p)
    std::vector<float> pf_2_8_ppb = { 0, 0, 0, 11.4, 12.6, 12.9, 14.8, 17.5, 3.2, 2.1, 1.9, 1.6, 2.1, 1.7, 1.5, 1.9, 2.3, 3.5, 16.1, 12.5, 28.6, 12.9, 10.8, 0, 0, 0 }; // Run2, 8TeV pPb (p-> <-Pb)
    auto thr_p = &t0;
    switch (_run){
//       case 1: { // TODO find thresholds for maximum energy in tower from first run
//         thr_p = &pf_1;
//         break;
//       }
      case 2:{
        switch (_energy){
          case 5: {
            thr_p = (_p_to_castor) ?  &pf_2_5_pbp:  &pf_2_5_ppb;
            break;
          }
          case 8: {
            thr_p = (_p_to_castor) ?  &pf_2_8_pbp:  &pf_2_8_ppb;
            break;
          }
        }
      }
    }
    if(_thr_pfen_twre != *thr_p) _thr_pfen_twre = *thr_p;
    return _thr_pfen_twre;
  }

  inline void updateThresholds(const int energy = 8, const bool p_to_castor = true, const int run = 2){
    _energy = energy;
    _p_to_castor = p_to_castor;
    _run = run;
  }
  // END TODO

  // version from Reco!!
  constexpr float THR_ZDC_RECO_TOT_M = 100.;
  constexpr float THR_ZDC_REC0_TOT_P = 320.;
  constexpr float THR_ZDC_RECO_EM_M  = 10.;
  constexpr float THR_ZDC_REC0_EM_P  = 10.;

  // dummy value for Digis!!
  constexpr float THR_ZDC_DIGI_TOT_M = 0.;
  constexpr float THR_ZDC_DIGI_TOT_P = 0.;
  constexpr float THR_ZDC_DIGI_EM_M  = 0.;
  constexpr float THR_ZDC_DIGI_EM_P  = 0.;

  //if(zdc56){sig = q.at(5) + q.at(6) - q.at(3) - q.at(4);}
  //else     {sig = q.at(4) + q.at(5) - q.at(2) - q.at(3);};
  constexpr short unsigned int FSC_NNOTS = 3; // number of first TSs
  constexpr short unsigned int FSC_NSITS = 3; // number of first TSs
  constexpr float FSC_SITHR[8] = {2586.69, 2511.54, 2487.4, 2591.24, 2620.62, 2606.23,10000, 10000};

  inline float find_eta(const unsigned int bin){
      float eta=ETA_BIN_L[bin]+ETA_BIN_W/2;
      return eta;
  }

  inline int find_eta_bin(const double eta){
    if(fabs(eta)>ETA_MAX)
      return -1;
    unsigned int n;
    bool search = true;
    if(eta>=0){
      n=N_ETA_BINS;
      while(search && n>0){
        n--;
        search = (eta<ETA_BIN_L[n]);
      };
    }else{
      n = 0;
      while(search && n<N_ETA_BINS){
        search = (eta>ETA_BIN_L[n]);
        n++;
      };
      n = n-2;
    };
    return n;
  }
};

#endif
