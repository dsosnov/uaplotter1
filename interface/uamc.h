#ifndef UAMC_H
#define UAMC_H

#include "uabasecentral.h"

#include "MyGenPart.h"

class uamc :  public uabasecentral
{
public:
  uamc(TChain*              tree,     //!<tree of ua format
       TDirectory*          dir,      //!<directory in the output root file
       const bool           cmstotem, //!<true for merged, false for CMS only
       const bool           pPb,      //!< true for pPb (p->Castor)
       const short          MC,       //!< -1, 0 - data; >0 MC
       const unsigned short Ncuts     //!< number of cuts
    );
  ~uamc();
  void PrintEventInfo(const bool detailed);
  void PrintProtonInfo();
  bool FillLastEvent(const short unsigned int cut);
  bool ProceedEvent(const short unsigned int cut, const bool fill, const bool info);

  bool         GetT2trigger () const {return ( (trksT2[0]>0) || (trksT2[1]>0) ); };
  bool         GetT2trigger(bool plus) const {return trksT2[int(plus)]>0;};
  unsigned int GetNT2trk(bool plus) const {return trksT2[int(plus)];};

  double GetOuterE(bool plus) const {return outRangeE[int(plus)];};
  double GetOuterPz(bool plus) const {return outRangePz[int(plus)];};

  short  IntactProton()   const {return protonSign;};
  double IntactProtonE()  const {return protonE;};
  double IntactProtonPz() const {return protonPz;};
  double IntactProtonXi() const {return protonXi;};
  double IntactProtonPt() const {return protonPt;};

  double GetZDCEn(bool plus)  const {return eZDCn[int(plus)];};
  double GetZDCEg(bool plus)  const {return eZDCg[int(plus)];};
  double GetZDCPzn(bool plus) const {return pzZDCn[int(plus)];};
  double GetZDCPzg(bool plus) const {return pzZDCg[int(plus)];};

  double GetCastorE()  const {return castorE;};
  double GetCastorPz() const {return castorPz;};

  double GetHFE(bool side)   const {return hfE[int(side)]    ;};
  double GetHFPz(bool side)  const {return hfPz[int(side)]   ;};
  double GetHFmax(bool side) const {return hfE_max[int(side)];};
  double GetBlindSpotMax(bool side) const {return blindSpot_max[int(side)];};

  int GetProcessID() const { for(auto p: *MCthuth) if(p.processID) return p.processID; return 0; }

private:
  const bool ppb;

  std::vector<MyGenPart>*  MCthuth;

  double energyMax[N_ETA_BINS];
  double ptChargedMax[N_ETA_BINS];

  unsigned int trksT2[2];

  double castorE;
  double castorPz;

  double eZDCn[2];
  double eZDCg[2];
  double pzZDCn[2];
  double pzZDCg[2];

  double hfE[2];
  double hfPz[2];
  double hfE_max[2];
  double blindSpot_max[2]; // [3.0, 3.15)

  short int protonSign;
  double    protonXi;
  double    protonEta;
  double    protonPt;
  double    protonPz;
  double    protonE;

  double totalE;
  double totalPz;
  double totalXi;

  double outRangeE[2];
  double outRangePz[2];

  void create_histos();
  TH1F ** proton_pt2_h;    //!< pt^2 = -t2 (pt of p`)
  TH1F ** proton_e_h;      //!< E(p`)
  TH2F ** proton_xi_mc_h;  //!< X: xi(p`);  Y: xi(X)
  TH2F ** proton_pt2_xi_h; //!< X: pt^2;    Y: xi(p`)
//  TH2F ** proton_e_eta_h;  //!<  X: |eta|;  Y: E(p`)
//  TH2F ** proton_e_pt2_h;  //!<  X: pt^2;   Y: E(p`)

  ClassDef(uamc,2);
};

#endif // UAMC_H
