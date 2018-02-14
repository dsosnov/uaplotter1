/*
 *
 */

#include "uabasecentral.h"
#include "iostream"

ClassImp(uabasecentral)

uabasecentral::uabasecentral( const bool cmstotem,
                              const short int MC,
                              const short unsigned int Ncuts,
                              TDirectory* dir):uabase(cmstotem, MC, Ncuts, dir)
{
}

uabasecentral::~uabasecentral()
{

}

void uabasecentral::PrepareArrays()
{
  memset(activity_loose, false, sizeof(activity_loose));
  memset(activity_tight, false, sizeof(activity_tight));
  memset(energy, 0, sizeof(energy));
  memset(pt,     0, sizeof(pt));
  memset(pz,     0, sizeof(pz));
}

bool uabasecentral::FillLastEvent(const short unsigned int cut)
{
  if (cut >= n_cuts) {
    std::cout << "uabasecentral::FillLastEvent: required cut number is larger that possible, do nothing. Please define larger uabasecentral::n_cut!\n";
    return false;
  };
  for (unsigned int bin = 0; bin < N_ETA_BINS; bin++){
    if(activity_loose[bin])
      activity_loose_h[cut]->Fill(find_eta(bin));
    if(activity_tight[bin]) activity_tight_h[cut]->Fill(find_eta(bin));
  }
  return true;
}


void uabasecentral::PrintActivity(bool tight)
{
  if(tight){
    std::cout << "\ttight activity per eta bin:\n\t";
    for(unsigned short int bin=0; bin<N_ETA_BINS; bin++)
      std::cout << activity_loose[bin] << "  ";
  }else{
    std::cout << "\tloose activity per eta bin:\n\t";
    for(unsigned short int bin=0; bin<N_ETA_BINS; bin++)
      std::cout << activity_loose[bin] << "  ";
  };
  std::cout << std::endl;
}

void uabasecentral::create_histos()
{
  TString title1, title2;

  if( !(n_each_h1D&&n_each_h2D) ){
    n_each_h1D = n_cuts;
    n_each_h2D = n_cuts;
  }

  activity_loose_h      = new TH1F * [n_each_h1D];
  activity_tight_h      = new TH1F * [n_each_h1D];

  for (unsigned int i = 0; i < n_each_h1D; i++) {
    title1 = "activity_loose_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ; #eta";
    activity_loose_h[i] = new TH1F(title1.Data(), title2.Data(), 28, -7, 7);
    activity_loose_h[i]->SetDirectory(directory);

    title1 = "activity_tight_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ; #eta";
    activity_tight_h[i] = new TH1F(title1.Data(), title2.Data(), 28, -7, 7);
    activity_tight_h[i]->SetDirectory(directory);
  }

  h1D->push_back(activity_loose_h);
  h1D->push_back(activity_tight_h);
}
