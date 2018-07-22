
#include "uaplotter1.h"

#include "TString.h"
#include "TMath.h"

#include "iostream"
#include "stdlib.h"

#include <algorithm>

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*!
 * \param hlt_path_regexps vector of regexps for HLT paths, that should be enabled at the same time.
 */
int uaplotter1::Loop(const int evts, const int trigger, vector<string> hlt_path_regexps) {
  if (mc < 0)
    return noiseLoop(evts);

  bool tech_bit   = true;
  int trigger_bit = DefineTrigger(trigger, tech_bit); // just translate initial trigger value

  int stat = chainTree->GetEntries();
  unsigned int nevts = evts;
  if (evts == -1 || evts > stat) nevts = stat;
  std::cout << "Total stat = " << stat << std::endl;
  std::cout << "uaplotter1::Loop(" << trigger_bit << ") for " << nevts << " events\n";

  unsigned int selected_evts[n_cuts];
  memset(selected_evts, 0, sizeof(selected_evts));

  unsigned int trigger_evts = 0;
  unsigned int t2prim_evts  = 0;
  unsigned int bptx_active  = 0;
  unsigned int hf_inelastic_count = 0;
  unsigned int sd_minus     = 0;
  unsigned int sd_plus      = 0;
  unsigned int sd_minus_eta[11];
  unsigned int sd_plus_eta[11];
  unsigned int rg_minus[28];     // in all topology id's
  unsigned int rg_plus[28];      // in all topology id's
  unsigned int rg_minus_hfp[28]; // in all topology id's
  unsigned int rg_plus_hfm[28];  // in all topology id's
  unsigned int goodFSC      = 0;
  unsigned int elastic_cand = 0;
  unsigned int cd_cand      = 0;
  unsigned int nd_cand      = 0;
  unsigned int dd_cand      = 0;
  memset(sd_minus_eta, 0, sizeof(sd_minus_eta));
  memset(sd_plus_eta,  0, sizeof(sd_plus_eta));
  memset(rg_minus, 0, sizeof(rg_minus));
  memset(rg_plus,  0, sizeof(rg_plus));
  memset(rg_minus_hfp, 0, sizeof(rg_minus_hfp));
  memset(rg_plus_hfm,  0, sizeof(rg_plus_hfm));

  short unsigned int bins_merged = 1; // Number of bins will be merged for filling histograms. "1" is histograms for every rg bin, "2" is standard "per one \eta" //TODO
  short unsigned int cuts_count = int(ceil((last_central_bin - first_central_bin + 1 + 1) / double(bins_merged)));

  double hf_inelastic_cut = 7.0;   // in GeV
  double hf_veto_cut      = 0.005; // in GeV -- seems like empty on particle level
  double hf_veto_cut_data = 3.0; // in GeV -- seems like empty on detector level

  unsigned int kevt = 0;
  std::vector<bool>hlts;
  for (long unsigned int i = 0; i < nevts; i++) {
    unsigned int kevt_current = i / 1000;
    if (kevt_current > kevt) {
      kevt = kevt_current;
      std::cout << "Event: " << i << std::endl;
    };
    chainTree->GetEntry(i);
    current_event = i;

    memset(sd_flag_central, processID::pid_undefined, sizeof(sd_flag_central));
    memset(sd_flag_total,   processID::pid_undefined, sizeof(sd_flag_total));
    memset(hf_inelastic,    false,                    sizeof(hf_inelastic));
    memset(hf_emptyHF,      false,                    sizeof(hf_emptyHF));
    memset(hf_blindFilled,      false,                    sizeof(hf_blindFilled));

    if (mc > 0) { // <============================  do MC loop here
      CMSmc->ProceedEvent(dummy_cut, false, false);
    };

    if (CMSevtinfo->GetL1Bit(9)) continue; //<=== if bptx quiet
    bptx_active++;

    hlts.resize(hlt_path_regexps.size(), false);
    std::transform(hlt_path_regexps.begin(), hlt_path_regexps.end(), hlts.begin(), [this](string s) {return this->CMSevtinfo->CheckHLT(s.c_str(), 1);});
    if ( !( ProceedTrigger(trigger_bit, tech_bit) &&
            std::all_of(hlts.begin(), hlts.end(), [](bool b) {return b == true;})) ) continue; //Check for all hlt paths and l1 bits
    trigger_evts++;

    CMScalo->ProceedEvent(dummy_cut, false, false);
    hf_inelastic[0] = CMScalo->GetHFmax(0) > hf_inelastic_cut;
    hf_inelastic[1] = CMScalo->GetHFmax(1) > hf_inelastic_cut;
    if(mc > 0){
      hf_emptyHF[0] = CMSmc->GetHFmax(0) < hf_veto_cut;
      hf_emptyHF[1] = CMSmc->GetHFmax(1) < hf_veto_cut;
      hf_blindFilled[0] = CMSmc->GetBlindSpotMax(0) > hf_inelastic_cut;
      hf_blindFilled[1] = CMSmc->GetBlindSpotMax(1) > hf_inelastic_cut;
    } else {
      hf_emptyHF[0] = CMScalo->GetHFmax(0) < hf_veto_cut_data;
      hf_emptyHF[1] = CMScalo->GetHFmax(1) < hf_veto_cut_data;
    }
    if ( !hf_inelastic[0] && !hf_inelastic[1] ) continue; //<=== If maximumTower in both HF lesser then 4GeV (inelastic cut)
    hf_inelastic_count++;

//     // <==================================== T2 check
//     bool t2prim = false;
//     if (mc > 0) {
//       if (tech_bit && (trigger_bit == 53)) {
//         t2prim = true; // we already checked in the ProceedTrigger; for MC it is the same!
//       } else {
//         t2prim = ProceedTrigger(53, true);
//       };
//     } else if (tree_combined_flag && (mc == 0)) {
//       T2->ProceedEvent(dummy_cut, false, false);
//       t2prim = ((T2->NPrimtracksMinus() > 0) || (T2->NPrimtracksPlus() > 0));
//     };
//     if (t2prim) t2prim_evts++;
//     // <=====================================
//
    ProceedEvent(dummy_cut, false, false);
//     if (tree_digi_flag && CMSforward->FSCvalid()) goodFSC++;

    FillLastEvent(0); // -------------------------------->  all triggered events

    if(hf_inelastic[0]) FillLastEvent(9);
    if(hf_inelastic[1]) FillLastEvent(10);
    if(hf_emptyHF[1] && hf_inelastic[0]) FillLastEvent(11);
    if(hf_emptyHF[0] && hf_inelastic[1]) FillLastEvent(12);

    switch(sd_flag_total[0]){
      case processID::pid_elastic:
//         CalculateSDdiffMass(false);
        FillLastEvent(3);
        elastic_cand++;
        break;
      case processID::pid_sdm:
//         CalculateSDdiffMass(false);
        FillLastEvent(4);
        sd_minus++;
        break;
      case processID::pid_sdp:
//         CalculateSDdiffMass(false);
        FillLastEvent(5);
        sd_plus++;
        break;
      case processID::pid_cd:
        FillLastEvent(6);
        cd_cand++;
        break;
      case processID::pid_nd:
        FillLastEvent(7);
        nd_cand++;
        break;
      case processID::pid_dd:
        FillLastEvent(8);
        dd_cand++;
        break;
    }

    int cuts_start = 13;
    auto sd_minus_cut = static_cast<unsigned int>(floor(n_sd_minus_bins[0] / static_cast<double>(bins_merged)));
    auto sd_plus_cut  = static_cast<unsigned int>(floor(n_sd_plus_bins[0] / static_cast<double>(bins_merged)));
    if(hf_inelastic[1]){ // HF>(inelastic cut) at plus side (RG at minus side).
      FillLastEvent(cuts_start + sd_minus_cut); //
      rg_minus_hfp[sd_minus_cut]++;
    }
    if(hf_inelastic[0]){  // HF>(inelastic cut) at minus side (for RG at plus side).
      FillLastEvent(cuts_start + sd_plus_cut + cuts_count);
      rg_plus_hfm[sd_plus_cut]++;
    }

  };// end loop
  std::cout << "Acceptance: [" << ETA_BIN_L[first_central_bin] << "," << ETA_BIN_L[last_central_bin] + ETA_BIN_W << "]\n";
  std::cout << "Total evts in chain       : " << stat << std::endl;
  std::cout << "Proceeded evts            : " << (current_event + 1)     << std::endl;
  std::cout << "Active (!bptx quiet)      : " << bptx_active << std::endl;
  std::cout << "Triggered evts            : " << trigger_evts << std::endl;
  std::cout << "After inelastic cut       : " << hf_inelastic_count << std::endl;
  std::cout << "\telastic candidates      : " << elastic_cand << "\t\tnd candidates       : " << nd_cand << std::endl;
  std::cout << "\tcd candidates           : " << cd_cand      << "\t\tdd candidates       : " << dd_cand << std::endl;
  std::cout << "\tsd- candidates          : " << sd_minus     << "\t\tsd+ candidates      : " << sd_plus << std::endl;

  std::cout   << "\n****                In all events                    ****" << std::endl;
  std::cout   << "\t\t" << "    Minus side" << "\t" << "     Plus side" << std::endl;
  for (short unsigned int ii = 0; ii < cuts_count; ii++) {
    std::cout << ii << " <= |deta| < " << ii + 1 << "\t\t" << rg_minus_hfp[ii] << "\t\t" << rg_plus_hfm[ii] << std::endl;
  };
  return nevts;
}
