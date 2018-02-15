#ifndef UAFORWARD_H
#define UAFORWARD_H

#include "uabase.h"

#include "MyZDCDigi.h"
#include "MyFSCDigi.h"

#include "TObject.h"
#include "TChain.h"

/* ZDC and FSC treatment
 * for digis only (at the moment there is no reliable reco)
 * IMPORTANT for ZDC: side = -1 for +Z!
 * ZDC is configured to have firstSample be TS5 in runs 210498-210676, and TS4 in runs 210737-211831
 */
class uaforward: public uabase{

public:
  uaforward(TChain     * tree, 	 //!<tree of ua format
	    TDirectory * dir,    //!<directory in the output root file
	    const bool cmstotem, //!<true for merged, false for CMS only
	    const bool ZDC56,    //!<flag for different TS
	    const short          int MC,  //!< -1, 0 - data; >0 MC
	    const short unsigned int Ncuts //!< number of cuts
	   );
  ~uaforward();
  bool ProceedEvent(const short unsigned int cut=0, const bool fill=false, const bool info=false);
  bool FillLastEvent(const short unsigned int cut);
  void PrintEventInfo(const bool detailed = false);
  double GetZDCEtotal(bool plus){return Etotal[int(plus)];};
  void FillZDCWithMCtruth(const short unsigned int cut, const double E[2], const double EM[2]); //!< bypass for MC
  void NormalizeFSCts();

  bool FSCvalid(){return goodFSC;};
  float GetFSCmSignal8(){return FSCm_Si8-FSCm_No8;};
  short unsigned int GetFSCmN(){short unsigned int n=0; for(short unsigned int i=0;i<8;i++){if(FSCm_flag[i]>0)n++;};return n;}
private:
  const bool zdc56;              //!< flag for ZDCdigi proper TS, 
  short unsigned int startFSCts;
  std::vector<MyZDCDigi>*    ZDCdigis;
  std::vector<MyFSCDigi>*    FSCdigis;

  void create_histos();

  TH1F ** ZDCd_Etot_minus_h;
  TH1F ** ZDCd_Etot_plus_h;
  TH1F ** ZDCd_EM_minus_h;
  TH1F ** ZDCd_EM_plus_h;

  TH1F ** FSCm_TS_ha[8];
  TH1F ** FSCm_No_ha[8];
  TH1F ** FSCm_Si_ha[8];
  TH1F ** FSCm_N_h;
  TH1F ** FSCm_No8_h;
  TH1F ** FSCm_Si8_h;

  TH2F ** ZDCm_FSCmSi8_h;
  TH2F ** ZDCm_FSCmN_h;

  float   FSCm_TS[8][10];  // 8 stations 10 TSs
  float   FSCm_No[8];
  float   FSCm_Si[8];
  bool    FSCm_flag[8];
  float   FSCm_Si8;
  float   FSCm_No8;
  bool    goodFSC;

  double Etotal[2];
  double EEM[2];
  ClassDef(uaforward,3);
};

#endif // UAFORWARD_H
