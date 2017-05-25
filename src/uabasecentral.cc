/*
 *
 */

#include "uabasecentral.h"
#include "iostream"

ClassImp(uabasecentral)

uabasecentral::uabasecentral(const bool cmstotem, 
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
