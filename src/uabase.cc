/*
 *
 */

#include "uabase.h"
#include "iostream"

ClassImp(uabase)



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uabase::uabase(const bool cmstotem, const short int MC, const short unsigned int Ncuts, TDirectory *dir):
  tree_combined_flag(cmstotem),
  mc(MC),
  n_cuts(Ncuts),
  n_each_h1D(1),
  n_each_h1Da(0),
  n_each_h2D(1),
  directory(dir)
{
  h1D = new std::vector<TH1F **>;
  h2D = new std::vector<TH2F **>;
  //std::cout << "uabase::uabase(Ncuts=" << n_cuts << ")" << std::endl;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uabase::~uabase()
{
  std::cout << "uabase::~uabase()\n";
  if (h1D != 0) {
    std::cout << "uabase::~uabase() : vector 1D size: " << h1D->size() << std::endl;
    for (std::vector<TH1F **>::iterator i = h1D->begin(); i != h1D->end(); ++i)
      delete_histos((*i));
    delete h1D;
  };
  if (h2D != 0) {
    std::cout << "uabase::~uabase() : vector 2D size: " << h2D->size() << std::endl;
    for (std::vector<TH2F **>::iterator i = h2D->begin(); i != h2D->end(); ++i)
      delete_histos((*i));
    delete h2D;
  };
  if (h1Da != 0) {

  }
  if (directory != 0)
    delete directory;
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uabase::delete_histos(TH1F **h)
{
  if (h) {
    for (unsigned int i = 0; i < n_each_h1D; i++) delete h[i];
    delete [] h;
//     std::cout << "ok\n";
  }
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uabase::delete_histos(TH2F **h)
{
  if (h) {
    for (unsigned int i = 0; i < n_each_h2D; i++) delete h[i];
    delete [] h;
//     std::cout << "ok\n";
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uabase::resete_histos(const unsigned int cut)
{
  if (cut > n_cuts) {
    std::cout << "uabase::resethistos(" << cut << ") for n_cuts=" << n_cuts << "\n";
    return;
  };
}
