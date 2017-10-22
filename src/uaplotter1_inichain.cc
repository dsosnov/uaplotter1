#include "uaplotter1.h"

#include "iostream"
#include "stdlib.h"

// also sets zdc56
TString uaplotter1::initializeChain()
{

  TString str = "";

  switch (mc) {
    // NOISE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    case -1: {
      zdc56 = false;
      // for new data
      if (ppb) {
        str = "noise_PARun2016_pPb";
        chainTree->Add("uatree_noise_pPb.root");
      } else {
        str = "noise_PARun2016_Pbp";
        chainTree->Add("uatree_noise_Pbp.root");
      }
      break;
     };

    case 0: { // DATA ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zdc56 = false;
      // for new data
      if (ppb) {
        str = "data_PARun2016_pPb";
        chainTree->Add("uatree_data_pPb.root");
      } else {
         str = "data_PARun2016_Pbp";
         chainTree->Add("uatree_data_Pbp.root");
      }
      break;
    };

    case 1: { // Hijing ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zdc56 = false;
      if (ppb) {
        str = "Hijing_pPb";
        chainTree->Add("uatree_hijing_pPb.root");
      } else {
        str = "Hijing_Pbp";
        chainTree->Add("uatree_hijing_Pbp.root");
      };
      break;
    };

    case 2: { // EPOS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zdc56 = false;
      if (ppb) {
        str = "EPOS_pPb";
        chainTree->Add("uatree_epos_pPb.root");
      } else {
        str = "EPOS_Pbp";
        chainTree->Add("uatree_epos_Pbp.root");
      };
      break;
    };

    case 3: { // QGSJET ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zdc56 = false;
      if (ppb) {
        str = "QGSJET_pPb";
        chainTree->Add("uatree_qgsjet_pPb.root");
      } else {
        str = "QGSJET_Pbp";
        chainTree->Add("uatree_qgsjet_Pbp.root");
      };
      break;
    };

    default: {
      zdc56 = false;
      std::cout << "uaplotter1::initializeChain(): unknown value for MC type = " << mc << std::endl;
      std::cout << "uaplotter1::initializeChain(): the TChain will not be initialized; exiting.\n";
      exit(10);
    };
  };
  return str;
};
