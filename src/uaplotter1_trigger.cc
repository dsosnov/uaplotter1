
#include "uaplotter1.h"
#include "uathresholds.h"

#include "TString.h"
#include "TMath.h"

#include "iostream"
#include "stdlib.h"

#define L1_Run2
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //! translates the user input (XX - tech bits; 1YY - 100+L1 trigger)
/*!
\param trigger_bit_common given trigger word: XX - tech bits; 2YY - 100+L1 trigger
\param tech_bit           set here, flag if the trigger bit was technical trigger, othervise algo trigger
\return trigger bit value (regardless tt or l1) or RNDT=-100 if mc<0 (random trigger), or -1 if not valid.
*/

int uaplotter1::DefineTrigger(int trigger_bit_common, bool &tech_bit){
  // ramdom trigger
  if(mc<0)
    return RNDT;

  int trigger_bit = trigger_bit_common;
#ifndef L1_Run2
  if(trigger_bit>64){
    tech_bit = false;
    trigger_bit-=100;
  };
  if(trigger_bit<0 || trigger_bit>127) trigger_bit = -1;

  std::cout << "uaplotter1::DefineTrigger: trigger selection on ";
  if(trigger_bit!=-1){
    if(tech_bit) {
      std::cout << "tech bit # ";
    }else{
      std::cout << "algo bit # ";
    }; std::cout << trigger_bit << std::endl;
  }else{
    std::cout << "...NOT FOUND, running on everything!\n";
  };
#else
  tech_bit = false;
#endif
  return trigger_bit;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //! Check trigger for DATA or MC (at the moment only tt53 for MC)
/*!
\param trigger_bit trigger bit; if < 0 => returns true always  // TODO check. Always or if not random trigger?
\param tech_bit    flag if the trigger bit was technical trigger, othervise algo trigger
\return trigger bit status (DATA) or emulation (MC, tt53) or noise studies conditions
*/
bool uaplotter1::ProceedTrigger(int trigger_bit, bool tech_bit) // TODO find invocations and disable work with bit 53
{
  bool proceed = true;

  if(mc<0){/////////////////// random trigger                                          
    proceed = ( (!CMSevtinfo->CheckHLT("HLT_L1Tech53_MB")) // TODO check
                && (!CMSevtinfo->GetL1Bit(10)) // ! L1Tech_BPTX_plus.v0
                && (!CMSevtinfo->GetL1Bit(1)) // ! L1_BptxPlus_NotBptxMinus
                // && (!CMSevtinfo->GetAlgoBit(2)) // ! L1_BeamGas_Hf_BptxPlusPostQuiet
                && (CMSevtinfo->GetL1Bit(9)) //
                // && (!CMSevtinfo->GetTechBit(53))   // ! L1Tech_TOTEM_MinBias.v0  // paranojac
      );
  }else{
#ifndef NEW_L1
    if(trigger_bit>=0) {
      if(tech_bit){///////////// tech trigger 
        if( (trigger_bit==53) && (mc>0) ){ // T2 trigger with MC //TODO disable
          proceed = CMSmc->GetT2trigger(); //TODO to view on this function
        }else{
          proceed   = CMSevtinfo->GetL1Bit(trigger_bit);
        }
      }else{ ////////////////// phys trigger
        proceed   = CMSevtinfo->GetL1Bit(trigger_bit);
      }
    }
#else
    proceed   = CMSevtinfo->GetL1Bit(trigger_bit);
#endif
  }
  return proceed;
}
