#ifndef UATRACKING_H
#define UATRACKING_H

#include "uabasecentral.h"

#include "MyVertex.h"
#include "MyTracks.h"

#include "TH1F.h"

#include "vector"

class uatracking : public uabasecentral {
public:
  uatracking(TChain      *tree,    //!<tree of ua format
             TDirectory *dir,     //!<directory in the output root file
             const bool cmstotem, //!<true for merged, false for CMS only
             const short          int MC,  //!< -1, 0 - data; >0 MC
             const short unsigned int Ncuts //!< number of cuts
            );
  ~uatracking();
  bool ProceedEvent(const short unsigned int cut = 0, const bool fill = false, const bool info = false);
  bool FillLastEvent(const short unsigned int cut);
  void PrintEventInfo(const bool detailed = false);

  unsigned int NverticesGood() {return nVtxGoog;};
  unsigned int Ntracks() {return nTracks;}
  unsigned int NtracksGood() {return nTracksGood;};
private:
  std::vector<MyVertex>     *Vertices;
  std::vector<MyTracks>     *Tracks;

  unsigned int nVtx;
  unsigned int nVtxGoog;
  unsigned int nTracks;
  unsigned int nTracksGood;

  unsigned int n_tracks_bin_all [N_ETA_BINS];
  unsigned int n_tracks_bin_good[N_ETA_BINS];
  float ptMax_bin_all [N_ETA_BINS];
  float ptMax_bin_good[N_ETA_BINS];

  void create_histos();

  TH1F **vertices_h;
  //TH1F ** vertices_X_h;
  //TH1F ** vertices_Y_h;
  //TH1F ** vertices_Z_h;
  TH1F **tracks_h;
  TH1F **tracks_pt_h;
  TH2F **tracks_eta_n_h;
  TH2F **tracks_eta_pt_h;

  ClassDef(uatracking, 2);

};

#endif // UATRACKING_H
