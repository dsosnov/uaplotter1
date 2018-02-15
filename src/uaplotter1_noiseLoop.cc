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

  unsigned int rndHLT    = 0;
  unsigned int not53     = 0;
  unsigned int bptxQuiet = 0;
  unsigned int bptxXor   = 0;

  unsigned int kevt = 0;

  for (long unsigned int i = 0; i < nevts; i++) {
    unsigned int kevt_current = i / 1000;
    if (kevt_current > kevt) {
      kevt = kevt_current;
      std::cout << ">>> proceeding " << kevt << std::endl;
    };
    chainTree->GetEntry(i);
    current_event = i;

    // THIS WILL BE VALID FOR NEW NOISE DATA ONLY! (mc=-2)
    rndHLT++;
    if (true /*!CMSevtinfo->GetTechBit(53)*/) {
      not53++;
      bool isQuiet = CMSevtinfo->GetL1Bit(9);
      bool isXor = CMSevtinfo->GetL1Bit(1) ||  CMSevtinfo->GetL1Bit(2);
      if (isQuiet || isXor){
        ProceedEvent(dummy_cut, false, false);
//         CMStracking->ProceedEvent(dummy_cut,false,true);
        if(isXor)                     FillLastEvent(0);
        if(CMSevtinfo->GetL1Bit(1))   FillLastEvent(1); // L1_BptxPlus_NotBptxMinus
        if(CMSevtinfo->GetL1Bit(2))   FillLastEvent(2); // L1_BptxMinus_NotBptxPlus
        if(isQuiet)                   FillLastEvent(3);
        if (CMStracking->NtracksGood() == 0){ // Without good tracks
          if(isQuiet) bptxQuiet++; else bptxXor++;
          if(isXor)                   FillLastEvent(4);
          if(CMSevtinfo->GetL1Bit(1)) FillLastEvent(5); // L1_BptxPlus_NotBptxMinus
          if(CMSevtinfo->GetL1Bit(2)) FillLastEvent(6); // L1_BptxMinus_NotBptxPlus
          if(isQuiet)                 FillLastEvent(7);
        }
      }
    }
  }
  std::cout << "total number of events         " << stat      << std::endl;
  std::cout << "worked over                    " << nevts     << std::endl;
  std::cout << "number of rndHLT events        " << rndHLT    << std::endl;
  std::cout << "rndHLT without tt53            " << not53     << std::endl;
  std::cout << "noBPTX without good tracks     " << bptxQuiet << std::endl;
  std::cout << "BPTX_XOr without good tracks   " << bptxXor   << std::endl;
  PrintEventInfo(true);
  return bptxQuiet+bptxXor;
}
