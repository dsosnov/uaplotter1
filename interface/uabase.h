#ifndef UABASE_H
#define UABASE_H

#include "TObject.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"

class uabase :  public TObject {
public:
  uabase(const bool               cmstotem = true,          //!< → #tree_combined_flag
         const short int          MC       = 0,             //!< → #mc
         const short unsigned int Ncuts    = 0,             //!< → #n_cuts
         TDirectory              *dir      = NULL           //!< → #directory
        );
  ~uabase();

protected:
  const bool               tree_combined_flag;              //!< true - with totem, false cms only
  const short int          mc;                              //!< \li -3 — Noise for Run2 \li -2 — Noise \li -1 — ? \li 0 — data, \li 1 — Hijing, \li 2 — EPOS
  const short unsigned int n_cuts;                          //!< total numbers of cuts supposed
  short unsigned int       n_each_h1D;                      //!< number of 1D histos per one cut
  short unsigned int       n_each_h1Da;                     //!< number of 1D histos per one cut
  short unsigned int       n_each_h2D;                      //!< number of 2D histos per one cut
  TDirectory              *directory;                       //!< root directory for histos
  std::vector<TH1F **>    *h1D;                             //!< vector of 1D histos[]
  std::vector<TH1F **>    *h1Da;                            //!< vector of 1D histos[]
  std::vector<TH2F **>    *h2D;                             //!< vector of 2D histos[]

  virtual bool ProceedEvent(const short unsigned int cut  = 0,
                            const bool               fill = false,
                            const bool               info = false
                           ) = 0;
  virtual bool FillLastEvent(const short unsigned int cut) = 0;
  virtual void PrintEventInfo(const bool detailed = false) = 0;
private:

  void delete_histos(TH1F **h);                             //!< deletes 1D histos (for ~uabase)
  void delete_histos(TH2F **h);                             //!< deletes 2D histos (for ~uabase)
  void resete_histos(const unsigned int cut);               //!<resets histos (TBD)

  virtual void create_histos() = 0;

  ClassDef(uabase, 2);
};

#endif // UABASE_H
