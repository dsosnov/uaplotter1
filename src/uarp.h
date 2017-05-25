#ifndef __UARP_H__
#define __UARP_H__

#include "uabase.h"

#include "RPRootDumpReconstructedProton.h"
#include "RPRootDumpTrackInfo.h"

#include <TObject.h>
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"

#include "vector"
#include "string"

/* RP treatment
 * CASTORp = true; // RP right
 * CASTORp = false; // RP - left
 */

class uarp: public uabase{
public :
  uarp(TChain     * tree, 	 //!<tree of ua format
       TDirectory * dir,    //!<directory in the output root file
       const bool cmstotem, //!<true for merged, false for CMS only
       const bool CASTORp,  //!<true proton to -Z (pPb)
       const short          int MC,  //!< -1, 0 - data; >0 MC
       const short unsigned int Ncuts //!< number of cuts
      );
  ~uarp();
  
  double  Xi()   {return proton_xi;};
  double  t()    {return proton_t;};
  bool    Valid(){return proton_valid;};
  bool trackValidUp(){return track_valid[0];};
  bool trackValidDn(){return track_valid[1];};
  float FriciVar();
  bool ProceedEvent(const short unsigned int cut=0, const bool fill=false, const bool info=false);
  bool FillLastEvent(const short unsigned int cut);
  void PrintEventInfo(const bool detailed = false);
  void FillRPWithMCtruth(const short unsigned int cut, const double xi, const double t); //!< bypass for MC

private:
  double proton_xi;
  double proton_t;
  double proton_y;
  bool   proton_valid;
  double y[2][2]; // [F,N][U,D]
  bool   track_valid[2]; // [U,D]
  RPRootDumpReconstructedProton*     RPproton;
  RPRootDumpTrackInfo*               RPtrackNU;
  RPRootDumpTrackInfo*               RPtrackND;
  RPRootDumpTrackInfo*               RPtrackFU;
  RPRootDumpTrackInfo*               RPtrackFD;
  
  void create_histos();
  TH1F ** xi_valid_h; 
  TH1F ** t_valid_h; 
  TH2F ** xi_vs_t_valid_h; 
  TH2F ** friciU_h;
  TH2F ** friciD_h;
  TH2F ** yp_vs_yN;
  
  ClassDef(uarp,2);
  
};

#endif