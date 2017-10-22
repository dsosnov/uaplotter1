#ifndef UACMSEVTINFO_H
#define UACMSEVTINFO_H

#include "uabase.h"

#include "MyEvtId.h"
#include "MyL1TrigRun2.h"
#include "MyL1MenuRun2.h"
#include "MyHLTrig.h"

const std::vector<unsigned int>L1_TRIGGER_ARRAY = {
  0,   // L1_ZeroBias
  1,   // L1_BptxPlus_NotBptxMinus
  2,   // L1_BptxMinus_NotBptxPlus
  3,   // L1_MinimumBiasHF0_OR_BptxAND
  4,   // L1_MinimumBiasHF0_AND_BptxAND
  5,   // L1_MinimumBiasHF0_OR
  6,   // L1_MinimumBiasHF0_AND
  9,   // L1_NotBptxOR
  10,  // L1_BptxPlus
  11,  // L1_BptxMinus
  209, // L1_BptxOR
  220  // L1_BptxXOR
};

// const unsigned int ALGO_TRIGGER_ARRAY[3] = {99,  // L1_CastorEm_TotemLowMultiplicity
//              105, // L1_DoubleJet20_TotemDiffractive
//              122};// L1_HcalHfSingleChannel_BptxAND_Instance1 (passive)
// const unsigned int TT_TRIGGER_ARRAY[5]   ={ 0, //  L1Tech_BPTX_plus_AND_minus.v0 (passive) //Run2: 10 && 11 or Or anything else
//              7, // L1Tech_BPTX_quiet.v0 (passive) // Run2: 9
//              9, //   L1Tech_HCAL_HF_coincidence_PM.v2 (passive)
//              52, // L1Tech_TOTEM_Diffractive.v0 (passive)
//              53};  // L1Tech_TOTEM_MinBias.v0
class uacmsevtinfo :  public uabase {
public:
  uacmsevtinfo(TChain      *tree,     //!<tree of ua format
               TDirectory *dir,     //!<directory in the output root file
               const bool cmstotem, //!<true for merged, false for CMS only
               const short          int MC,  //!< -1, 0 - data; >0 MC
               const short unsigned int Ncuts //!< number of cuts);
              );
  ~uacmsevtinfo();
  bool GetL1Bit(short unsigned int b) {return CMStrigInfo->triggerResults.at(b);};
  bool GetL1Bit(string name) {
    for(auto i: CMSmenuInfo->menu)
      if(!name.compare(i.second.name))
        return CMStrigInfo->triggerResults.at(i.first);
    return false;
  };
  vector<bool> GetL1AllBits() {return CMStrigInfo->triggerResults;};
  /*!
   * \param path
   * \param compareType type of comparing: \li 0 - find any occurrence \li 1 - paths should be equal \li 2 - regex compare
   * \return if this path was in dataset
   */
  bool CheckHLT(const char *path, int compareType = 0);
  void PrintEventInfo(const bool detailed = false);
  bool FillLastEvent(const short unsigned int cut);
  bool ProceedEvent(const short unsigned int cut, const bool fill, const bool info);

private:
  void create_histos();

  MyEvtId      *CMSevtInfo;
  MyL1TrigRun2 *CMStrigInfo;
  MyL1MenuRun2 *CMSmenuInfo;
  MyHLTrig     *CMSHLT;

  TH1F **triggers_h;
  TH2F **run_vs_bx_h;
  TH2F **run_vs_ls_h;

  ClassDef(uacmsevtinfo, 3);
};

#endif // UACMSEVTINFO_H
