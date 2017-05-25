
#include "uaplotter1.h"
#include "uathresholds.h"

#include "TString.h"
#include "TMath.h"

#include "iostream"
#include "stdlib.h"


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
  return trigger_bit;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //! Check trigger for DATA or MC (at the moment only tt53 for MC)
/*!
\param trigger_bit trigger bit; if < 0 => returns true always
\param tech_bit    flag if the trigger bit was technical trigger, othervise algo trigger
\return trigger bit status (DATA) or emulation (MC, tt53) or noise studies conditions
*/
bool uaplotter1::ProceedTrigger(int trigger_bit, bool tech_bit)
{
  bool proceed = true;
  
  if(mc<0){/////////////////// random trigger                                          
    proceed = ( (!CMSevtinfo->CheckHLT("HLT_L1Tech53_MB")) && 
		(!CMSevtinfo->GetTechBit(1)) && (!CMSevtinfo->GetAlgoBit(2)) && 
	        (CMSevtinfo->GetTechBit(7)) && (!CMSevtinfo->GetTechBit(53)) // paranojac
	      );
  }else if(trigger_bit>=0) {
    if(tech_bit){///////////// tech trigger
      if( (trigger_bit==53) && (mc>0) ){ // T2 trigger with MC
	proceed = CMSmc->GetT2trigger();
      }else{
	proceed = CMSevtinfo->GetTechBit(trigger_bit);
      };
    }else{ ////////////////// phys trigger
      proceed   = CMSevtinfo->GetAlgoBit(trigger_bit);
    };
  };
  return proceed;
}
