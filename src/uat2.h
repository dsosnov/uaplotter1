#ifndef __uat2_H__
#define __uat2_H__

#include "T2Event.h"
#include "TriggerData.h"// TBD to be moved outside

#include <TObject.h>
#include "TDirectory.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"

#include "vector"
#include "string"

class uat2:public TObject{
  public :
    uat2(TChain*tree, 
 	 TDirectory * dir,    //!<directory in the output root file
	 const bool combined=true, const int Ncuts=2);
    ~uat2();
  
    int ProceedEvent(unsigned int cut=0, const bool fill=false, const bool info=false);
    int Plot_histots(bool D1=true);
    void PrintPrimTracks(bool minus=true);
    void PrintAllTracks(bool minus=true);
    
    unsigned int NtracksPlus(){return n_tracks_plus;};
    unsigned int NtracksMinus(){return n_tracks_minus;};
    unsigned int NPrimtracksPlus(){return n_prim_tracks_plus;};
    unsigned int NPrimtracksMinus(){return n_prim_tracks_minus;};
    
    TH1F ** n_tracks_plus_h; 
    TH1F ** n_tracks_minus_h;
    TH1F ** n_tracks_prim_plus_h; 
    TH1F ** n_tracks_prim_minus_h;
    TH1F ** z_impact_plus_h;
    TH1F ** z_impact_minus_h;
    TH1F ** z_impact_prim_plus_h;
    TH1F ** z_impact_prim_minus_h;
    // trigger problem studies
    TH1F *** n_clusters_plus_h; 
    TH1F *** n_clusters_minus_h;
    TH1F *** n_clusters_plusF_h; 
    TH1F *** n_clusters_minusF_h;
    TH1F *** n_clusters_plusN_h; 
    TH1F *** n_clusters_minusN_h;
    TH2F *** n_cluster_2h;
    TH2F *** n_qtracks_2h;
    TH1F *** trigger_h; // TBD to be moved outside
    
  void delete_histo(TH1F **h,  int N);
  void delete_histo(TH2F ***h, int N, int M);
  void delete_histo(TH1F ***h, int N, int M);
  void resethistos(unsigned int cut);
  
  T2Event*    T2tracks;
  T2Event*    T2primaryPlus;
  T2Event*    T2primaryMinus;
  TriggerData* TOTtrig;// TBD to be moved outside
private:
    TDirectory * directory;
    unsigned int n_cuts;
    unsigned int n_t2_h;
    unsigned int n_tracks_plus;
    unsigned int n_tracks_minus;
    unsigned int n_qtracks_plus[5]; // number of tracks per a quarter
    unsigned int n_qtracks_minus[5];// number of tracks per a quarter 
    unsigned int n_prim_tracks_minus;
    unsigned int n_prim_tracks_plus;
    
    static const unsigned short int NTCASES = 5; 
    unsigned int n_clust_minus[5];
    unsigned int n_clust_plus[5];
    
    bool trigger[15]; // TBD to be moved outside
    bool T2oneArm;    // TBD to be moved outside
    bool T2trig;      // TBD to be moved outside
    void addT2primaryTrack(unsigned int i, unsigned int cut=0, const bool fill=false, const bool info = false);
    void resetT2primary(bool plus);
    ClassDef(uat2,3);
};

#endif