#ifndef __uaetabin_H__
#define __uaetabin_H__
#include "math.h"

const short int RNDT      = -100;

//this is wrong for MC
//const double ENERGYN      = 3999.98095727174496; // energy in GeV from sqrt(S_NN) = 5.023
const double MOMBEAM        = 3997.7190062; // From EPOS

const float HF_ETA_MAX    = 5.1; // for HF "trigger" only

const float CAS_ETA_MIN   = -6.6;
const float CAS_ETA_MAX   = -5.2;

const float CAS_ETA_AVE   = -5.9;
const double CAS_TANH_AVE = -9.99984990996805823e-01;

const float T2_ABSETA_MIN = 5.3;
const float T2_ABSETA_MAX = 6.5;
const float T2_PT_THR     = 0.04;

const float CENT_TOWER_THR = 2.0;   // energy in GeV
const float HF_TOWER_THR   = 5.0;   // energy in GeV
const float ZDC_PXY_THR    = 0.495;   // hernja!!!
const float MIN_ZDC_ETA    = 8.; //8.2?

const float MC_E_THR       = 0.;    // energy in GeV

const float CENT_ETA_ACC   = 2.0;   // eta
//const float CENT_ETA_ACC   = 3.0; // 16Feb15 - test 

const float CALO_ETA_ACC   = 5;   // eta *was 5

const float ETA_BIN_W = 0.5;
const float ETA_MAX   = 6.5;//CALO_ETA_ACC;//ETA_BIN_L[N_ETA_BINS-1]+ETA_BIN_W;

const unsigned int N_ETA_BINS = ETA_MAX * 2 / ETA_BIN_W /*26*/;
const float ETA_BIN_L[N_ETA_BINS] = {-6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6};

const int BINMIN = 3;  /* (ETA_MAX - CALO_ETA_ACC) / ETA_BIN_W / */
const int BINMAX = 22; /* (ETA_MAX + CALO_ETA_ACC) / ETA_BIN_W - 1 */

// for range [-5,5]
//const float THR_CALO[N_ETA_BINS]  = {   0., 0.,   0., 3.,   4.,  4.,  4.,  13., 10., 12., 15., 9., 9., 9.,  9., 14.,12., 9.,12.,4.,5., 5., 4., 0., 0., 0.};    
//new
const float THR_CALO[N_ETA_BINS]  = {   0., 0.,   0.,   2.9,  3.9, 3.4,  3.1, 12.8, 9.3, 11.2, 14.8, 8.4, 8.9, 8.2,  8.0, 13.7, 11.5, 8.8, 11.9, 3.4, 4.0, 4.7, 3.7, 
  0., 0., 0.};    
  // neutral particles
//const float THR_PFEN[N_ETA_BINS]  = {   0., 0.,   0., 2.5,  3.,  3., 2.5, 16.5, 10., 9.,  7.5, 3., 5., 4.,  3., 6., 9.5, 9.,16.,3.,3.5,3.5,3., 0., 0., 0.};
// new
  const float THR_PFEN[N_ETA_BINS]  = {   0., 0.,   0., 
  2.2,  2.5,  2.7, 2.4, 16.4, 9.5, 8.9,  7.0, 2.9, 4.7, 3.8,  2.6, 5.7, 9.2, 8.9, 15.5, 2.8, 3.0, 3.1,2.9, 
  0., 0., 0.};

// version from Reco!!
const float THR_ZDC_RECO_TOT_M = 100.;
const float THR_ZDC_REC0_TOT_P = 320.;
const float THR_ZDC_RECO_EM_M  = 10.;
const float THR_ZDC_REC0_EM_P  = 10.;

// dummy value for Digis!!
const float THR_ZDC_DIGI_TOT_M = 0.;
const float THR_ZDC_DIGI_TOT_P = 0.;
const float THR_ZDC_DIGI_EM_M  = 0.;
const float THR_ZDC_DIGI_EM_P  = 0.;

//if(zdc56){sig = q.at(5) + q.at(6) - q.at(3) - q.at(4);}
//else     {sig = q.at(4) + q.at(5) - q.at(2) - q.at(3);};
const short unsigned int FSC_NNOTS = 3; // number of first TSs
const short unsigned int FSC_NSITS = 3; // number of first TSs
const float FSC_SITHR[8] = {2586.69, 2511.54, 2487.4, 2591.24, 2620.62, 2606.23,10000, 10000};

inline float find_eta(unsigned int bin){
    float eta=ETA_BIN_L[bin]+ETA_BIN_W/2;
    return eta;
};

inline int find_eta_bin(double eta){
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
      //std::cout << "\t\t\t" << n << "   " << search << std::endl;
      n++;
    }; n = n-2; 
  };
  //std::cout << "\t" << eta << "\t\t" << n << std::endl;
  return n;
};



#endif
