#include "uaplotter1.h"
#include "uathresholds.h"

#include "TString.h"
#include "TMath.h"

#include "iostream"
#include "stdlib.h"


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! noise studies on rnd trigger
/*!
\param evts number of events to be prodused
\return number of events of the quiet BPTX
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int uaplotter1::noiseLoop(const int evts)
{
  if (mc >= 0) {
    std::cout << "can't run noiseLoop on MC/data files...\n";
    return 0;
  };

  int stat = chainTree->GetEntries();
  unsigned int nevts = evts;
  if (evts == -1 || evts > stat) nevts = stat;
  std::cout << "Total stat = " << stat << std::endl;
  std::cout << "uaplotter1::noiseLoop() for " << nevts << " events\n";

  unsigned int rndHLT       = 0;
  unsigned int not53        = 0;
  unsigned int not53notBPTX = 0;

  unsigned int kevt = 0;

  for (long unsigned int i = 0; i < nevts; i++) {
    unsigned int kevt_current = i / 1000;
    if (kevt_current > kevt) {
      kevt = kevt_current;
      std::cout << ">>> proceeding " << kevt << std::endl;
    };
    chainTree->GetEntry(i);
    current_event = i;


    memset(sd_flag_central, 0, sizeof(sd_flag_central));
    memset(sd_flag_total,   0, sizeof(sd_flag_total));

    // THIS WILL BE VALID FOR NEW NOISE DATA ONLY! (mc=-2)
    //if(CMSevtinfo->CheckHLT("HLT_PARandom_v1")){ they are not set :(
    rndHLT++;
    if(true /*!CMSevtinfo->GetTechBit(53)*/){
      not53++;
      if(CMSevtinfo->GetL1Bit(9)){
	not53notBPTX++;
	ProceedEvent(0, false, false);
	FillLastEvent(0);
      }
    };
    //};

  }; // end loop
  std::cout << "total number of events         " << stat   << std::endl;
  std::cout << "worked over                    " << nevts  << std::endl;
  std::cout << "number of rndHLT events        " << rndHLT << std::endl;
  std::cout << "rndHLT without tt53            " << not53  << std::endl;
  std::cout << "rndHLT without tt53 and noBPTX " << not53notBPTX << std::endl;
  return not53notBPTX;
}
