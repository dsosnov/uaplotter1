#ifndef UAPF_H
#define UAPF_H

#include "uabasecentral.h"

#include "MyPFCand.h"

class uapf :  public uabasecentral
{
public:
    uapf(TChain     * tree, 	 //!<tree of ua format
	   TDirectory * dir,    //!<directory in the output root file
	   const bool cmstotem, //!<true for merged, false for CMS only
	   const short          int MC,  //!< -1, 0 - data; >0 MC
	   const short unsigned int Ncuts //!< number of cuts
	  );
   ~uapf();
    void PrintEventInfo(const bool detailed);
    bool FillLastEvent(const short unsigned int cut);
    bool ProceedEvent(const short unsigned int cut, const bool fill, const bool info);
private:
  std::vector<MyPFCand>*     PFCand;
  unsigned int nPfCand;

  void create_histos();
  TH1F ** pfcand_h;
  TH2F ** pf_e_eta_h;
  //TH2F ** pf_pt_eta_h;
  TH2F ** pf_e_hch_h;
  TH2F ** pf_e_h0_h;
  TH2F ** pf_e_em0_h;
  double energyH0[N_ETA_BINS];
  double energyHch[N_ETA_BINS];
  double energyEM0[N_ETA_BINS];

  ClassDef(uapf,2);
};

#endif // UAPF_H
