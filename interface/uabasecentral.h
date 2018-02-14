#ifndef UABASECENTRAL_H
#define UABASECENTRAL_H

#include "uabase.h"
#include "uathresholds.h"
class uabasecentral :  public uabase
{
public:
    uabasecentral(const bool cmstotem, const short int MC, const short unsigned int Ncuts, TDirectory* dir);
    ~uabasecentral();
    bool GetActivityLoose(short unsigned int bin){return activity_loose[bin];};
    bool GetActivityTight(short unsigned int bin){return activity_tight[bin];};
    double GetE(short unsigned int bin){return energy[bin];};
    double GetPz(short unsigned int bin){return pz[bin];};
protected:
    virtual void PrintEventInfo(const bool detailed) =0;
    bool FillLastEvent(const short unsigned int cut);
    virtual bool ProceedEvent(const short unsigned int cut, const bool fill, const bool info)=0;
    void PrepareArrays();
    void PrintActivity(bool tight);
    void create_histos();
    bool activity_loose[N_ETA_BINS];
    bool activity_tight[N_ETA_BINS];
    double energy[N_ETA_BINS];
    double pz[N_ETA_BINS];
    double pt[N_ETA_BINS];

    TH1F **activity_loose_h;
    TH1F **activity_tight_h;

    ClassDef(uabasecentral,2);
};

#endif // UABASECENTRAL_H
