#include "uaplotter1.h"

#include "uamc.h"

#include "iostream"
#include "stdlib.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//sets first_central_bin and last_central_bin according to uathresholds
void uaplotter1::IniRapGapRange()
{
  short unsigned int bin = 0;
  while (ETA_BIN_L[bin] < -RG_ETA_ACC)bin++;
  first_central_bin = bin;
  bin = N_ETA_BINS - 1;
  while (ETA_BIN_L[bin] + ETA_BIN_W > RG_ETA_ACC)bin--;
  last_central_bin = bin;
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uaplotter1::FindRapGap(bool RECO)
{

  bool rg = false;

  int ind = int(!RECO);
  if (RECO)
    memset(combined_central_activity, false, sizeof(combined_central_activity));

  n_sd_minus_bins[ind] = 0; // [0]-RECO, [1]-MCtruth
  n_sd_plus_bins[ind]  = 0;
  n_dd_rg_bins[ind]    = 0;

  bool diff_cand[3] = {true, false, false}; // -, central, +
  short unsigned int rg_bins[3]   = {0,   0,   0};     //
  short unsigned int max_central_rgbins = 0;
  short unsigned int n_active_bins   = 0;

  for (unsigned int bin = first_central_bin; bin <= last_central_bin; bin++) {

    bool binactivity;
    if (RECO) { //<====================== RG on RECO
      if(7 < bin && bin < 18) // tracking works good only in [-2.5,2.5]
        combined_central_activity[bin] = (CMSpf->GetActivityLoose(bin) || CMStracking->GetActivityTight(bin));   // GetActivityLoose or GetActivityTight
      else
        combined_central_activity[bin] = CMSpf->GetActivityLoose(bin);
      binactivity = combined_central_activity[bin];
    } else { //=========================== RG on MCtruth
      binactivity = CMSmc->GetActivityLoose(bin);
    };

    if (!binactivity) {
      if (diff_cand[0]) {
        rg_bins[0]++;
      } else {
        diff_cand[2] = true;
        rg_bins[2]++;
      };
    } else {
      n_active_bins++;
      if (diff_cand[0]) {
        diff_cand[0] = false;
      } else {
        if (diff_cand[2]) {
          diff_cand[1] = true;
          rg_bins[1]   = rg_bins[2];
          if (rg_bins[1] > max_central_rgbins) max_central_rgbins = rg_bins[1];
          diff_cand[2] = false;
          rg_bins[2]   = 0;
        };
      };
    }

  };//end loop

  if (n_active_bins == 0) {                           // "~elastic"
    sd_flag_central[ind] = processID::pid_elastic;
    rg_bins[2] = rg_bins[0];
  } else if ((rg_bins[0] == 0) && (rg_bins[2] == 0)){
    if (max_central_rgbins > 1){
      sd_flag_central[ind] = processID::pid_dd;       // DD candidate (it was >1 earlies)
    } else if(max_central_rgbins == 0){               // ND
      sd_flag_central[ind] = processID::pid_nd;
    } else if(max_central_rgbins == 1){               // It was an special case, TODO check how it will be better
      sd_flag_central[ind] = processID::pid_dd;
    }
  } else if(rg_bins[2] == 0){
    sd_flag_central[ind] = processID::pid_sdm;        // SD-
  } else if(rg_bins[0] == 0){
    sd_flag_central[ind] = processID::pid_sdp;        // SD+
  } else if(rg_bins[0] > 0 && rg_bins[2] > 0){
    sd_flag_central[ind] = processID::pid_cd;         // CD candidate
  }

  n_sd_minus_bins[ind] = rg_bins[0];
  n_sd_plus_bins[ind]  = rg_bins[2];
  n_dd_rg_bins[ind]    = max_central_rgbins;

  sd_flag_total[ind] = sd_flag_central[ind];
// WARNING temporary disabled finding of total sd flag
//   ///////////////////////////////////////////////
//   // take into account also sides
//   bool outer_activity_minus;
//   bool outer_activity_plus;
//   if (RECO) {
//     if (mc > 0) {
//       outer_activity_minus = CMSmc->GetT2trigger(false); // minus
//       outer_activity_plus  = CMSmc->GetT2trigger(true);  // plus
//     } else if (mc == 0 && tree_digi_flag) {
//       outer_activity_minus = (T2->NPrimtracksMinus() > 0);
//       outer_activity_plus  = (T2->NPrimtracksPlus() > 0);
//     } else { // fake case noise
//       outer_activity_minus = false;
//       outer_activity_plus  = false;
//     };
//     TotalRapGap(0, outer_activity_minus, outer_activity_plus);
//   } else {
//     outer_activity_minus = false;
//     outer_activity_plus  = false;
//     for (unsigned int bin = 0; bin < first_central_bin; bin++)
//       outer_activity_minus = (outer_activity_minus || CMSmc->GetActivityLoose(bin));
//     for (unsigned int bin = N_ETA_BINS - 1; bin > last_central_bin; bin--)
//       outer_activity_plus = (outer_activity_plus || CMSmc->GetActivityLoose(bin));
//     TotalRapGap(1, outer_activity_minus, outer_activity_plus);
//
//     outer_activity_minus = (outer_activity_minus || CMSmc->GetOuterE(false));
//     outer_activity_plus = (outer_activity_plus || CMSmc->GetOuterE(true));
//     TotalRapGap(2, outer_activity_minus, outer_activity_plus);
//   };


  // below is just to prepare values for diff mass calculation and make sure the values are ok.
  memset(xi_pf,   0, sizeof(xi_pf));
  memset(xi_calo, 0, sizeof(xi_calo));
  memset(xi_cas,  0, sizeof(xi_cas));
  memset(xi_zdc,  0, sizeof(xi_zdc));
  memset(xi_full, 0, sizeof(xi_full));

  xi_mc_out   = 0; // this includes ZDC and CASTOR also!
  xi_mc_total = 0; // everything from first active bin in reco
//   cout << ":2 RapidityGaps type: " << sd_flag_central[ind] << " -> " <<  sd_flag_total[ind] << std::endl;
  return rg;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uaplotter1::TotalRapGap(short unsigned int ind, bool active_minus, bool active_plus)
{
  bool rg = false;
  //std::cout << "=====> " << ind << "  " << active_minus << " " << active_plus << std::endl;
  switch (sd_flag_central[ind]) {
    case processID::pid_nd: // ND
      if (!active_minus) {
        if (!active_plus) {
          sd_flag_total[ind] = processID::pid_cd;
        } else {
          sd_flag_total[ind] = processID::pid_sdm;
        };
      } else {
        if (!active_plus) {
          sd_flag_total[ind] = processID::pid_sdp;
        } else {
          sd_flag_total[ind] = processID::pid_nd;
        }
      };
      break;
    case processID::pid_dd:
      if (!active_minus) {
        if (!active_plus) {
          sd_flag_total[ind] = processID::pid_cd;
        } else {
          sd_flag_total[ind] = processID::pid_sdm;
        };
      } else {
        if (!active_plus) {
          sd_flag_total[ind] = processID::pid_sdp;
        } else {
          sd_flag_total[ind] = processID::pid_dd;
        }
      };
    case processID::pid_cd:
      if (!active_minus) {
        if (!active_plus) {
          sd_flag_total[ind] = processID::pid_cd;
        } else {
          sd_flag_total[ind] = processID::pid_sdm;
        };
      } else {
        if (!active_plus) {
          sd_flag_total[ind] = processID::pid_sdp;
        } else {
          sd_flag_total[ind] = processID::pid_dd;
        }
      };
      break;
    case processID::pid_sdm:
      if (!active_minus) {
        sd_flag_total[ind] = processID::pid_sdm;
      } else {
        sd_flag_total[ind] = processID::pid_dd;
      };
      break;
    case processID::pid_sdp:
      if (!active_plus) {
        sd_flag_total[ind] = processID::pid_sdp;
      } else {
        sd_flag_total[ind] = processID::pid_dd;
      };
      break;
    case processID::pid_elastic:
      if (!active_minus) {
        if (!active_plus) {
          sd_flag_total[ind] = processID::pid_elastic;
        } else {
          sd_flag_total[ind] = processID::pid_sdm;
        };
      } else {
        if (!active_plus) {
          sd_flag_total[ind] = processID::pid_sdp;
        } else {
          sd_flag_total[ind] = processID::pid_dd;
        }
      };
      break;
    case processID::pid_undefined:
      cout << "uaplotter1::TotalRapGap: undefined PID\n";
      break;
  }; // end switch
  rg = (sd_flag_total[ind] != processID::pid_undefined);
  return rg;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uaplotter1::PrintRapGap()
{
  std::cout << "combined central activity:\n\t";
  for (unsigned int bin = 0; bin < N_ETA_BINS; bin++) {
    std::cout << combined_central_activity[bin] << "   ";
  }; std::cout << std::endl;
  std::cout << "diffractive central flag: " << sd_flag_central[0] << "; MC: " << sd_flag_central[1] << std::endl;
  std::cout << "diffractive total   flag: " << sd_flag_total[0] << "; MC: " << sd_flag_total[1] << std::endl;
  std::cout << "max sd bins - : "    << n_sd_minus_bins[0];
  if (mc > 0) std::cout << " (" << n_sd_minus_bins[1] << ")";
  std::cout << "; max sd bins + : " << n_sd_plus_bins[0];
  if (mc > 0) std::cout << " (" << n_sd_plus_bins[1]  << ")";
  std::cout << "; max dd bins : " << n_dd_rg_bins[0];
  if (mc > 0) std::cout << " (" << n_dd_rg_bins[1] << ")";
  std::cout << std::endl;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uaplotter1::CalculateSDdiffMass(bool info)
{

  memset(xi_pf,   0, sizeof(xi_pf));
  memset(xi_calo, 0, sizeof(xi_calo));
  memset(xi_cas,  0, sizeof(xi_cas));
  memset(xi_zdc,  0, sizeof(xi_zdc));
  memset(xi_full, 0, sizeof(xi_full));

  xi_mc_out   = 0; // this includes ZDC and CASTOR also!
  xi_mc_total = 0; // everything from first active bin in reco

  //==================================
  // ZDC contribution

  // -1=> rg at -side, X @ +Z => sign -; considered ZDC+: plus = true
  // +1=> rg at +side, X @ -Z => sign +; considered ZDC-: plus = false
  int sign_ = 0;
  switch (sd_flag_total[0]){
    case processID::pid_sdm:       sign_ = -1; break;
    case processID::pid_nd :       sign_ =  0; break;
    case processID::pid_sdp:       sign_ =  1; break;
    case processID::pid_dd :       sign_ =  2; break;
    case processID::pid_cd :       sign_ =  3; break;
    case processID::pid_elastic:   sign_ =  4; break;
    case processID::pid_undefined: 
      cout << "uaplotter1::CalculateSDdiffMass: undefined PID\n";
      break; 
  }
  bool Xplus = !(sign_ > 0);
  // zdc for X side xi ~ (E-Pz)/2Ep ~ 0 (if we don't know pz)
  // zdc for opposite side we can't use
  if (mc > 0) {
    xi_zdc[0] = 0;//2*(CMSmc->GetZDCEn(Xplus)+CMSmc->GetZDCEg(Xplus));
    xi_zdc[1] = CMSmc->GetZDCEn(Xplus) + CMSmc->GetZDCEg(Xplus) + sign_ * (CMSmc->GetZDCPzn(Xplus) + CMSmc->GetZDCPzg(Xplus));
//     std::cout << "ZDC:" << CMSmc->GetZDCEn(Xplus)+CMSmc->GetZDCEg(Xplus) << "\t"
//      << (CMSmc->GetZDCPzn(Xplus)+CMSmc->GetZDCPzg(Xplus)) << std::endl;
  } else {
    xi_zdc[0] = 0;//2*CMSforward->GetZDCEtotal(Xplus);
  };

  //==================================
  // castor
  // if SD = -1 => makes no sense!
  // if SD = +1
  if (sign_ > 0)
    xi_cas[0] = CMScastor->GetE() * (1 + sign_ * CAS_TANH_AVE);

  if (mc > 0) xi_cas[1] = CMSmc->GetCastorE() + sign_ * CMSmc->GetCastorPz();

  //==================================
  // central region
  short unsigned int start = first_central_bin;
  short unsigned int stop  = last_central_bin;
  if (Xplus) { // SD-
    start = first_central_bin + n_sd_minus_bins[0]; // based on RECO RG!
  } else {     // SD+
    stop  = last_central_bin  - n_sd_plus_bins[0];  // based on RECO RG!
  };

  for (unsigned int bin = start; bin <= stop; bin++) {
    if ((ETA_BIN_L[bin] >= -1 * CENT_ETA_ACC) && ((ETA_BIN_L[bin] + ETA_BIN_W) <= CENT_ETA_ACC)) {
      xi_pf[0] += (CMSpf->GetE(bin) + sign_ * CMSpf->GetPz(bin));
      if (mc > 0)
        xi_pf[1] += (CMSmc->GetE(bin) + sign_ * CMSmc->GetPz(bin));
    } else {
      xi_calo[0] += (CMScalo->GetE(bin) + sign_ * CMScalo->GetPz(bin));
      if (mc > 0)
        xi_calo[1] += (CMSmc->GetE(bin) + sign_ * CMSmc->GetPz(bin));
    };
  };

  if (mc > 0) {
    //std::cout << "1)" << xi_mc_out/(4000.*2) << std::endl;
    float eee = 0; float pzzz = 0;
    if (Xplus) { // RG @ -; X @ +
      for (unsigned int bin = N_ETA_BINS - 1; bin > last_central_bin; bin--) {
        xi_mc_out += (CMSmc->GetE(bin) + sign_ * CMSmc->GetPz(bin));
        eee += CMSmc->GetE(bin); pzzz += CMSmc->GetPz(bin);
      }
    } else {
      for (unsigned int bin = 0; bin < first_central_bin; bin++) {
        xi_mc_out += (CMSmc->GetE(bin) + sign_ * CMSmc->GetPz(bin));
        eee += CMSmc->GetE(bin); pzzz += CMSmc->GetPz(bin);
      };
    };
    //std::cout << "2)" << xi_mc_out/(4000.*2) << "\te " << eee << "\tpz " << pzzz << std::endl;
    xi_mc_out += (CMSmc->GetOuterE(Xplus) + sign_ * CMSmc->GetOuterPz(Xplus));
    //std::cout << "3)" << xi_mc_out/(4000.*2) << "\te " << eee+CMSmc->GetOuterE(Xplus) << "\tpz " << pzzz+CMSmc->GetOuterPz(Xplus) << std::endl;
    short int ip = CMSmc->IntactProton();
    if ((ip != 0) && (sign_ != ip)) { // proton at X side => add it to X
      xi_mc_out += (CMSmc->IntactProtonE() + sign_ * CMSmc->IntactProtonPz());
    };
  };


  if (info) {
    std::cout << "uaplotter1::CalculateSDdiffMass: sign_, Xplus : " << sign_ << " " << Xplus;
    if (mc > 0) {
      std::cout << "; IntactProton sign: " << CMSmc->IntactProton();
    };
    std::cout << std::endl;
    if (mc > 0)
      std::cout << "\tZDC MCtruth: Etot, Pztot:" << CMSmc->GetZDCEn(Xplus) + CMSmc->GetZDCEg(Xplus)
                << "; " << CMSmc->GetZDCPzn(Xplus) + CMSmc->GetZDCPzg(Xplus) << std::endl;
    std::cout << "start/stop bin for central activity: [" << start << " , " << stop << "]" << std::endl;
  };

  //==================================
  // normalization
  double denom = 2 * 4000.;
//   if( CASTORp ){ // CM to the positive side
//     denom = 2*4000.; // 2*Ep
//   }else{
//     denom = 2*4000.; // 2*Z(Pb)*Ep
//   };

  for (short unsigned int ii = 0; ii < 2; ii++) {
    xi_pf[ii] /= denom;
    xi_calo[ii] /= denom;
    xi_cas[ii] /= denom;
    xi_zdc[ii] /= denom;
    xi_full[ii]     = xi_pf[ii] + xi_calo[ii] + xi_cas[ii] + xi_zdc[ii];
  };
  xi_mc_out /= denom;
  xi_mc_total = xi_pf[1] + xi_calo[1] + xi_mc_out;

  /*
  if(fill){
    diffmass_pf_h[current_cut]->Fill(xsi_pf, xsi_mc_pf);
    diffmass_calo_h[current_cut]->Fill(xsi_calo, xsi_mc_calo);
    diffmass_castor_h[current_cut]->Fill(xsi_cas, xsi_mc_cas);
    diffmass_total_h[current_cut]->Fill(xsi_total, xsi_mc_total);
    diffmass_proton_h[current_cut]->Fill(xsi_mc_total+xsi_mc_out, CMSmc->GetXsiP());
  };
  */
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uaplotter1::PrintSDdiffMass(bool detailed)
{
  std::cout << "==> XI: \n";
  std::cout << "full visible RECO : " << xi_full[0];
  if (mc > 0) std::cout << ";  full visible MCtruth: " << xi_full[1] << "; total MCthruth: " << xi_mc_total;
  std::cout << std::endl;
  if (detailed) {
    std::cout << "\t pf   RECO: " << xi_pf[0];
    if (mc > 0)
      std::cout << "\t pf   MCtruth: " << xi_pf[1];
    std::cout << "\n\t calo RECO: " << xi_calo[0];
    if (mc > 0)
      std::cout << "\t calo MCtruth: " << xi_calo[1];
    std::cout << "\n\t cas  RECO: " << xi_cas[0];
    if (mc > 0)
      std::cout << "\t cas  MCtruth: " << xi_cas[1];
    std::cout << "\n\t zdc  RECO: " << xi_zdc[0];
    if (mc > 0)
      std::cout << "\t zdc  MCtruth: " << xi_zdc[1];
    std::cout << std::endl;
    if (mc > 0)
      std::cout << "\tMCtruth outer:" << xi_mc_out << std::endl;
  };
}
