#ifndef UACMSEVTINFO_H
#define UACMSEVTINFO_H

#include "uabase.h"

#include "MyEvtId.h"
#include "MyL1TrigOld.h"
#include "MyHLTrig.h"

const unsigned int ALGO_TRIGGER_ARRAY[3] = {99,  // L1_CastorEm_TotemLowMultiplicity
					    105, // L1_DoubleJet20_TotemDiffractive
					    122};// L1_HcalHfSingleChannel_BptxAND_Instance1 (passive)
const unsigned int TT_TRIGGER_ARRAY[5]   ={ 0, //	L1Tech_BPTX_plus_AND_minus.v0 (passive)
					    7, //	L1Tech_BPTX_quiet.v0 (passive)
					    9, // 	L1Tech_HCAL_HF_coincidence_PM.v2 (passive)
					    52, // L1Tech_TOTEM_Diffractive.v0 (passive)
					    53};
class uacmsevtinfo :  public uabase
{
public:
    uacmsevtinfo( TChain     * tree, 	 //!<tree of ua format
		  TDirectory * dir,    //!<directory in the output root file
		  const bool cmstotem, //!<true for merged, false for CMS only
		  const short          int MC,  //!< -1, 0 - data; >0 MC
		  const short unsigned int Ncuts //!< number of cuts);
		  );
    ~uacmsevtinfo();
    bool GetTechBit(short unsigned int b){return CMStrigInfo->TechTrigWord[b];};
    bool GetAlgoBit(short unsigned int b){return CMStrigInfo->PhysTrigWord[b];};
    bool CheckHLT(const char * path);
    void PrintEventInfo(const bool detailed);
    bool FillLastEvent(const short unsigned int cut);
    bool ProceedEvent(const short unsigned int cut, const bool fill, const bool info);

private:
    void create_histos();
    
    MyEvtId*      CMSevtInfo;
    MyL1TrigOld*  CMStrigInfo;
    MyHLTrig*     CMSHLT;
    
    TH1F ** triggers_h;
    TH2F ** run_vs_bx_h;
    TH2F ** run_vs_ls_h;

    ClassDef(uacmsevtinfo,3);
};

#endif // UACMSEVTINFO_H
