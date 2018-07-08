/*
 *
 */

#include "uaplotter1.h"

#include "TString.h"
#include "TMath.h"

#include "iostream"
#include "stdlib.h"

ClassImp(uaplotter1)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uaplotter1::uaplotter1(const bool               cmstotem,
                       const bool               cmsdigis,
                       const bool               CASTORp,
                       const short int          MC,
                       const short unsigned int Ncuts,
                       const int                energy
                      ):
  uabase(cmstotem, MC, Ncuts),
  tree_digi_flag(cmsdigis),
  ppb(CASTORp),
  dummy_cut(Ncuts + 1),
  current_event(0),
  CMSmc(NULL),
  CMSevtinfo(NULL),
  CMStracking(NULL),
  CMSpf(NULL),
  CMScalo(NULL),
  CMScastor(NULL),
  CMSforward(NULL),
  T2(NULL),
  RP(NULL),
  first_central_bin(0),
  last_central_bin(0),
  hf_by_processID_t(NULL)
  {
  updateThresholds(energy, ppb);
  IniRapGapRange();
  if (tree_combined_flag) {
    chainTree = new TChain("cms_totem");
  } else {
    chainTree = new TChain("evt");
  };
  TString str = initializeChain();
  str.Prepend("uaplot_output_histos_"); str += ".root";

  outputFile = TFile::Open(str.Data(), "RECREATE");

  if (mc > 0) {
    TDirectory *dirMC = outputFile->mkdir("MC");
    CMSmc              = new uamc(chainTree, dirMC, tree_combined_flag, ppb, mc, n_cuts);
  };

  TDirectory *dirEVT = outputFile->mkdir("CMSinfo");
  CMSevtinfo          = new uacmsevtinfo(chainTree, dirEVT, tree_combined_flag, mc, n_cuts);

  TDirectory *dirTRK = outputFile->mkdir("CMStracking");
  CMStracking         = new uatracking(chainTree, dirTRK, tree_combined_flag, mc, n_cuts);

  TDirectory *dirCAL = outputFile->mkdir("CMScalo");
  CMScalo             = new uacalo(chainTree, dirCAL, tree_combined_flag, mc, n_cuts);

  TDirectory *dirPF  = outputFile->mkdir("CMSpf");
  CMSpf               = new uapf(chainTree, dirPF, tree_combined_flag, mc, n_cuts);

  TDirectory *dirCAS = outputFile->mkdir("CMScastor");
  CMScastor           = new uacastor(chainTree, dirCAS, tree_combined_flag, tree_digi_flag, mc, n_cuts);

  if (tree_digi_flag || (mc > 0)) {
    TDirectory *dirZDC = outputFile->mkdir("ZDC");
    CMSforward          = new uaforward(chainTree, dirZDC, tree_combined_flag, zdc56, mc, n_cuts);
  };

  if (tree_combined_flag) {
    TDirectory *dirT2   = outputFile->mkdir("T2");
    T2                  = new uat2(chainTree, dirT2, tree_combined_flag, n_cuts);
  };

  if (tree_combined_flag || (mc > 0)) {
    TDirectory *dirRP  = outputFile->mkdir("RP");
    RP                  = new uarp(chainTree, dirRP, tree_combined_flag, ppb, mc, n_cuts);
  };

  //directory = outputFile->mkdir("main");
  create_histos();

//   hf_by_processID_t = new TTree("hf_by_processID","hf_by_processID");
  if (hf_by_processID_t != NULL){
    hf_by_processID_t->Branch("hfMinus",    &(hf_by_processID.hfMinus));
    hf_by_processID_t->Branch("hfPlus",     &(hf_by_processID.hfPlus));
    hf_by_processID_t->Branch("processID",  &(hf_by_processID.processID));
    hf_by_processID_t->Branch("l1Triggers", &(hf_by_processID.l1Triggers));
  }
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uaplotter1::~uaplotter1()
{
  std::cout << "uaplotter1::~uaplotter1()\n";
  if (tree_digi_flag || (mc > 0)) {
    CMSforward->NormalizeFSCts();
  };
  outputFile->Write();
  chainTree->Delete();
  if (mc > 0)
    delete CMSmc;
  delete CMSevtinfo;
  delete CMStracking;
  delete CMSpf;
  delete CMScalo;
  delete CMScastor;
  if (tree_digi_flag || (mc > 0)) {
    delete CMSforward;
  };
  if (tree_combined_flag)
    delete T2;
  if (tree_combined_flag || (mc > 0))
    delete RP;
  if (hf_by_processID_t != NULL)
    hf_by_processID_t->Delete();
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uaplotter1::ProceedEvent(const short unsigned int cut, const bool fill, const bool info)
{
  bool proceed_mc = (mc > 0);

  if (proceed_mc)
    CMSmc->ProceedEvent(cut, fill, false);

  if (fill)
    CMSevtinfo->ProceedEvent(cut, fill, false); // does nothing except for fill

  CMStracking->ProceedEvent(cut, fill, false);
  CMSpf->ProceedEvent(cut, fill, false);
  CMScalo->ProceedEvent(cut, fill, false);
  CMScastor->ProceedEvent(cut, fill, false);

  // ********** ZDC/FSC ******************
  if (tree_digi_flag) {
    CMSforward->ProceedEvent(cut, fill, false);
  } else if (proceed_mc && fill) {
    double E[2]  = {0, 0};
    double EM[2] = {0, 0};
    EM[0] = CMSmc->GetZDCEg(false);         EM[1] = CMSmc->GetZDCEg(true);
    E[0]  = EM[0] + CMSmc->GetZDCEn(false); E[1]  = EM[1] + CMSmc->GetZDCEn(true);
    CMSforward->FillZDCWithMCtruth(cut, E, EM);
  };

  // ************** T2 ********************
  if (tree_combined_flag && fill) // here to be optimized
    T2->ProceedEvent(cut, fill, false);

  // ************* forward dependencies ****************
  if (fill) {
    if (tree_digi_flag) {
      zdcM_vs_castor_h[cut]->Fill(CMScastor->GetE(),       CMSforward->GetZDCEtotal(false));
      zdcM_vs_T2primM_h[cut]->Fill(T2->NPrimtracksMinus(), CMSforward->GetZDCEtotal(false));
    } else if (proceed_mc) {
      zdcM_vs_castor_h[cut]->Fill(CMScastor->GetE(), CMSmc->GetZDCEg(false) + CMSmc->GetZDCEn(false));
      zdcM_vs_T2primM_h[cut]->Fill(CMSmc->GetNT2trk(false), CMSmc->GetZDCEg(false) + CMSmc->GetZDCEn(false));
    };
  };
  // ************** RP ********************
  if (tree_combined_flag || proceed_mc) {
    RP->ProceedEvent(cut, fill, false); // if MC does nothing
    if (fill && proceed_mc) {
      if (CMSmc->IntactProton() != 0) {
        double protonpt = CMSmc->IntactProtonPt();
        RP->FillRPWithMCtruth(cut, CMSmc->IntactProtonXi(), -protonpt * protonpt);
      };
    };
  };

  // ********** RapGap ******************
  FindRapGap(0); // RG in RECO
  if (proceed_mc){
    FindRapGap(1); // RG in MCtruth
    FindRapGap(2); // RG in MCtruth
  }

  // ************* SD ******************************************************
  if (sd_flag_central[0] != processID::pid_nd && (sd_flag_total[0] == processID::pid_sdm || sd_flag_total[0] == processID::pid_sdp)) {
    CalculateSDdiffMass();

    if (fill) {
      if (tree_combined_flag) {
        if (RP->Valid())
          xi_p_reco_full_h[cut]->Fill(RP->Xi(), xi_full[0]);              // DATA -> X: RP xi, Y: xi reco full
      } else if (proceed_mc && (CMSmc->IntactProton() != 0)) {
        xi_p_reco_full_h[cut]->Fill(CMSmc->IntactProtonXi(), xi_full[0]); // MC   -> X: proton xi, Y: xi reco full
      };
    };
  };

  // ********** Short information about MC events ******************
  if(cut==0){;
    hf_by_processID.hfMinus = CMScalo->GetHFmax(0);
    hf_by_processID.hfPlus = CMScalo->GetHFmax(1);
    hf_by_processID.processID = (mc>0) ? CMSmc->GetProcessID() : 0;
    hf_by_processID.l1Triggers = CMSevtinfo->GetL1AllBits();
  }

  // ***********************************************************************
  if (info)
    PrintEventInfo(true);
  return true;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uaplotter1::FillLastEvent(const short unsigned int cut)
{
  if (cut >= n_cuts) {
    std::cout << "uaplotter1::FillLastEvent: cut number [" << cut << "] larger then n_cuts [" << n_cuts << "]\n";
    return false;
  };
  bool ok = true;
  if (mc > 0)
    ok *= CMSmc->FillLastEvent(cut);
  ok *= CMSevtinfo->FillLastEvent(cut);
  ok *= CMStracking->FillLastEvent(cut);
  ok *= CMSpf->FillLastEvent(cut);
  ok *= CMScalo->FillLastEvent(cut);
  ok *= CMScastor->FillLastEvent(cut);
  if (tree_digi_flag) {
    ok *= CMSforward->FillLastEvent(cut);
  } else if (mc > 0) {
    double E[2]  = {0, 0};
    double EM[2] = {0, 0};
    EM[0] = CMSmc->GetZDCEg(false);       EM[1] = CMSmc->GetZDCEg(true);
    E[0]  = EM[0] + CMSmc->GetZDCEn(false); E[1]  = EM[1] + CMSmc->GetZDCEn(true);
    CMSforward->FillZDCWithMCtruth(cut, E, EM);
  };
  if (mc > 0) {
    diff_flag_mc_full_reco_central_h[cut]->Fill(sd_flag_total[1], sd_flag_central[0]);
    diff_flag_mc_full_reco_full_h[cut]->Fill(sd_flag_total[1], sd_flag_total[0]);
    diff_flag_mc_full_mc_central_h[cut]->Fill(sd_flag_total[1], sd_flag_central[1]);
    diff_flag_mc_total_mc_central_h[cut]->Fill(sd_flag_total[2], sd_flag_central[1]);
    n_sd_minus_bins_mcreco_mctruth_h[cut]->Fill(n_sd_minus_bins[0], n_sd_minus_bins[1]);
    if(hf_emptyHF[0])
      n_sd_minus_bins_mcreco_mctruth_emptyHF_h[cut]->Fill(n_sd_minus_bins[0], n_sd_minus_bins[1]);
    n_sd_plus_bins_mcreco_mctruth_h[cut]->Fill(n_sd_plus_bins[0], n_sd_plus_bins[1]);
    if(hf_emptyHF[1])
      n_sd_plus_bins_mcreco_mctruth_emptyHF_h[cut]->Fill(n_sd_plus_bins[0], n_sd_plus_bins[1]);
    n_sd_minus_bins_mcreco_mctruthLoose_h[cut]->Fill(n_sd_minus_bins[0], n_sd_minus_bins[2]);
    if(hf_emptyHF[0])
      n_sd_minus_bins_mcreco_mctruthLoose_emptyHF_h[cut]->Fill(n_sd_minus_bins[0], n_sd_minus_bins[2]);
    n_sd_plus_bins_mcreco_mctruthLoose_h[cut]->Fill(n_sd_plus_bins[0], n_sd_plus_bins[2]);
    if(hf_emptyHF[1])
      n_sd_plus_bins_mcreco_mctruthLoose_emptyHF_h[cut]->Fill(n_sd_plus_bins[0], n_sd_plus_bins[2]);
    n_sd_minus_bins_mctruth_mctruthLoose_h[cut]->Fill(n_sd_minus_bins[1], n_sd_minus_bins[2]);
    if(hf_emptyHF[0])
      n_sd_minus_bins_mctruth_mctruthLoose_emptyHF_h[cut]->Fill(n_sd_minus_bins[1], n_sd_minus_bins[2]);
    n_sd_plus_bins_mctruth_mctruthLoose_h[cut]->Fill(n_sd_plus_bins[1], n_sd_plus_bins[2]);
    if(hf_emptyHF[1])
      n_sd_plus_bins_mctruth_mctruthLoose_emptyHF_h[cut]->Fill(n_sd_plus_bins[1], n_sd_plus_bins[2]);
    for (short unsigned int bin = 0; bin < N_ETA_BINS; bin++)
      if (CMSmc->GetActivityLoose(bin))
        central_activity_mc_h[cut]->Fill(find_eta(bin));
  };
  n_sd_minus_bins_h[cut]->Fill(n_sd_minus_bins[0]);
  n_sd_plus_bins_h[cut]->Fill(n_sd_plus_bins[0]);
  n_sd_minus_bins_plus_bins_h[cut]->Fill(n_sd_minus_bins[0], n_sd_plus_bins[0]);
  for (short unsigned int bin = 0; bin < N_ETA_BINS; bin++)
    if (combined_central_activity[bin])
      central_activity_h[cut]->Fill(find_eta(bin));

  if (mc > 0) {
    if (CMSmc->IntactProton() != 0) {
      xi_mc_p_mc_total_h[cut]->Fill(TMath::Log10(CMSmc->IntactProtonXi()), TMath::Log10(xi_mc_total));
      xi_mc_p_reco_full_h[cut]->Fill(TMath::Log10(CMSmc->IntactProtonXi()), TMath::Log10(xi_full[0]));
    };
    xi_mc_total_mc_full_h[cut]->Fill(TMath::Log10(xi_mc_total), TMath::Log10(xi_full[1]));
    xi_mc_total_reco_full_h[cut]->Fill(TMath::Log10(xi_mc_total), TMath::Log10(xi_full[0]));
    xi_calo_mc_reco_h[cut]->Fill(TMath::Log10(xi_calo[1]), TMath::Log10(xi_calo[0]));
    xi_pf_mc_reco_h[cut]->Fill(TMath::Log10(xi_pf[1]), TMath::Log10(xi_pf[0]));
    xi_cas_mc_reco_h[cut]->Fill(TMath::Log10(xi_cas[1]), TMath::Log10(xi_cas[0]));
    xi_zdc_mc_reco_h[cut]->Fill(TMath::Log10(xi_zdc[1]), TMath::Log10(xi_zdc[0]));
  };
  xi_reco_full_h[cut]->Fill(TMath::Log10(xi_full[0]));

  if (tree_combined_flag) {
    T2->ProceedEvent(cut, true, false); // to be optimized later!
    RP->FillLastEvent(cut);
  } else if (mc > 0) {
    if (CMSmc->IntactProton() != 0) {
      double ptp = CMSmc->IntactProtonPt();
      RP->FillRPWithMCtruth(cut, CMSmc->IntactProtonXi(), -ptp * ptp);
    };
  };
  if (tree_digi_flag) {
    zdcM_vs_castor_h[cut]->Fill(CMScastor->GetE(),       CMSforward->GetZDCEtotal(false));
    zdcM_vs_T2primM_h[cut]->Fill(T2->NPrimtracksMinus(), CMSforward->GetZDCEtotal(false));
    FSCmN_vs_castor_h[cut]->Fill(CMScastor->GetE(), CMSforward->GetFSCmN());
  } else if (mc > 0) {
    zdcM_vs_castor_h[cut]->Fill(CMScastor->GetE(), CMSmc->GetZDCEg(false) + CMSmc->GetZDCEn(false));
    zdcM_vs_T2primM_h[cut]->Fill(CMSmc->GetNT2trk(false), CMSmc->GetZDCEg(false) + CMSmc->GetZDCEn(false));
  };

  if (sd_flag_central[0] != processID::pid_nd && (sd_flag_total[0] == processID::pid_sdm || sd_flag_total[0] == processID::pid_sdp)) {
    if (tree_combined_flag) {
      if (RP->Valid() && (RP->trackValidUp() || RP->trackValidDn()))
        xi_p_reco_full_h[cut]->Fill(RP->Xi(), xi_full[0]);              // DATA -> X: RP xi, Y: xi reco full
    } else if ((mc > 0) && (CMSmc->IntactProton() != 0)) {
      xi_p_reco_full_h[cut]->Fill(CMSmc->IntactProtonXi(), xi_full[0]); // MC   -> X: proton xi, Y: xi reco full
    };
  };

  if (tree_digi_flag && RP->Valid() && (RP->trackValidUp() || RP->trackValidDn())) {
    ZDCm_vs_xiRP_h[cut]->Fill(RP->Xi(), CMSforward->GetZDCEtotal(false));
    ZDCp_vs_xiRP_h[cut]->Fill(RP->Xi(), CMSforward->GetZDCEtotal(true));
    FSCmSi8_vs_xiRP_h[cut]->Fill(RP->Xi(), CMSforward->GetFSCmSignal8());
    FSCmN_vs_xiRP_h[cut]->Fill(RP->Xi(), CMSforward->GetFSCmN());
  };

  sd_flag_central_reco_h[cut]->Fill(sd_flag_central[0]);
  sd_flag_total_reco_h[cut]->Fill(sd_flag_total[0]);

  int pid = 0;
  if(mc>0){
    if(CMSmc->GetProcessID() > 100) pid = CMSmc->GetProcessID() - 100;
    n_sd_minus_bins_mctruth_pid_inelastic_h[cut]->Fill(n_sd_minus_bins[1], pid);
    n_sd_plus_bins_mctruth_pid_inelastic_h[cut]->Fill(n_sd_plus_bins[1], pid);
    n_sd_minus_bins_mctruthLoose_pid_inelastic_h[cut]->Fill(n_sd_minus_bins[2], pid);
    n_sd_plus_bins_mctruthLoose_pid_inelastic_h[cut]->Fill(n_sd_plus_bins[2], pid);
    if(hf_emptyHF[0]){
      n_sd_minus_bins_reco_pid_inelastic_veto_h[cut]->Fill(n_sd_minus_bins[0], pid);
      n_sd_minus_bins_mctruth_pid_inelastic_veto_h[cut]->Fill(n_sd_minus_bins[1], pid);
      n_sd_minus_bins_mctruthLoose_pid_inelastic_veto_h[cut]->Fill(n_sd_minus_bins[2], pid);
      if(hf_blindFilled[0]){
        n_sd_minus_bins_reco_pid_inelastic_veto_blindFilled_h[cut]->Fill(n_sd_minus_bins[0], pid);
        n_sd_minus_bins_mctruth_pid_inelastic_veto_blindFilled_h[cut]->Fill(n_sd_minus_bins[1], pid);
      }
    }
    if(hf_emptyHF[1]){
      n_sd_plus_bins_reco_pid_inelastic_veto_h[cut]->Fill(n_sd_plus_bins[0], pid);
      n_sd_plus_bins_mctruth_pid_inelastic_veto_h[cut]->Fill(n_sd_plus_bins[1], pid);
      n_sd_plus_bins_mctruthLoose_pid_inelastic_veto_h[cut]->Fill(n_sd_plus_bins[2], pid);
      if(hf_blindFilled[1]){      
        n_sd_plus_bins_reco_pid_inelastic_veto_blindFilled_h[cut]->Fill(n_sd_plus_bins[0], pid);
        n_sd_plus_bins_mctruth_pid_inelastic_veto_blindFilled_h[cut]->Fill(n_sd_plus_bins[1], pid);
      }
    }
  }
  n_sd_minus_bins_reco_pid_inelastic_h[cut]->Fill(n_sd_minus_bins[0], pid);
  n_sd_plus_bins_reco_pid_inelastic_h[cut]->Fill(n_sd_plus_bins[0], pid);

  // ********** Short information about MC events ******************
  if((hf_by_processID_t != NULL) && (cut==0)){;
    hf_by_processID_t->Fill();
  }

  return ok;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uaplotter1::PrintEventInfo(const bool detailed)
{
  std::cout << "\n\n" << current_event
            << "\t=====================================================\n";
  if (mc > 0)
    CMSmc->PrintEventInfo(detailed);
  CMSevtinfo->PrintEventInfo(detailed);
  CMStracking->PrintEventInfo(detailed);
  CMSpf->PrintEventInfo(detailed);
  CMScalo->PrintEventInfo(detailed);
  if (tree_digi_flag)
    CMSforward->PrintEventInfo(detailed);
  CMScastor->PrintEventInfo(detailed);
  if (mc == 0 && tree_combined_flag)
    std::cout << "T2 info : T2-: " << T2->NPrimtracksMinus() <<  "(" << T2->NtracksMinus()
              << ");\tT2+: " << T2->NPrimtracksPlus() << "(" << T2->NtracksPlus() << ")\n";
  if (tree_combined_flag) {
    RP->PrintEventInfo(detailed);
  };
  PrintRapGap();
  if (sd_flag_total[0] != processID::pid_nd && (sd_flag_total[0] == processID::pid_sdm || sd_flag_total[0] == processID::pid_sdp)) {
    PrintSDdiffMass(detailed);
    if (mc > 0)
      CMSmc->PrintProtonInfo();
  };
}




void uaplotter1::create_histos()
{
  n_each_h2D = n_cuts;
  n_each_h1D = n_cuts;
  TString title1, title2;

  auto pids_count = processID::pid_max - processID::pid_min + 1;
  if (mc > 0) { //******************************************************************************************
    diff_flag_mc_full_reco_central_h  = new TH2F * [n_each_h2D];
    diff_flag_mc_full_reco_full_h     = new TH2F * [n_each_h2D];
    diff_flag_mc_total_mc_central_h   = new TH2F * [n_each_h2D];
    diff_flag_mc_full_mc_central_h    = new TH2F * [n_each_h2D];

    n_sd_minus_bins_mcreco_mctruth_h         = new TH2F * [n_each_h2D];
    n_sd_plus_bins_mcreco_mctruth_h          = new TH2F * [n_each_h2D];
    n_sd_minus_bins_mcreco_mctruth_emptyHF_h = new TH2F * [n_each_h2D];
    n_sd_plus_bins_mcreco_mctruth_emptyHF_h  = new TH2F * [n_each_h2D];
    n_sd_minus_bins_mcreco_mctruthLoose_h         = new TH2F * [n_each_h2D];
    n_sd_plus_bins_mcreco_mctruthLoose_h          = new TH2F * [n_each_h2D];
    n_sd_minus_bins_mcreco_mctruthLoose_emptyHF_h = new TH2F * [n_each_h2D];
    n_sd_plus_bins_mcreco_mctruthLoose_emptyHF_h  = new TH2F * [n_each_h2D];
    n_sd_minus_bins_mctruth_mctruthLoose_h         = new TH2F * [n_each_h2D];
    n_sd_plus_bins_mctruth_mctruthLoose_h          = new TH2F * [n_each_h2D];
    n_sd_minus_bins_mctruth_mctruthLoose_emptyHF_h = new TH2F * [n_each_h2D];
    n_sd_plus_bins_mctruth_mctruthLoose_emptyHF_h  = new TH2F * [n_each_h2D];

    xi_mc_p_mc_total_h     = new TH2F * [n_each_h2D];
    xi_mc_p_reco_full_h    = new TH2F * [n_each_h2D];
    xi_mc_total_mc_full_h  = new TH2F * [n_each_h2D];
    xi_mc_total_reco_full_h = new TH2F * [n_each_h2D];

    xi_calo_mc_reco_h      = new TH2F * [n_each_h2D];
    xi_pf_mc_reco_h        = new TH2F * [n_each_h2D];
    xi_cas_mc_reco_h       = new TH2F * [n_each_h2D];
    xi_zdc_mc_reco_h       = new TH2F * [n_each_h2D];

    n_sd_minus_bins_mctruth_pid_inelastic_h      = new TH2F * [n_each_h2D];
    n_sd_plus_bins_mctruth_pid_inelastic_h       = new TH2F * [n_each_h2D];
    n_sd_minus_bins_reco_pid_inelastic_veto_h    = new TH2F * [n_each_h2D];
    n_sd_plus_bins_reco_pid_inelastic_veto_h     = new TH2F * [n_each_h2D];
    n_sd_minus_bins_mctruth_pid_inelastic_veto_h = new TH2F * [n_each_h2D];
    n_sd_plus_bins_mctruth_pid_inelastic_veto_h  = new TH2F * [n_each_h2D];
    n_sd_minus_bins_mctruthLoose_pid_inelastic_h      = new TH2F * [n_each_h2D];
    n_sd_plus_bins_mctruthLoose_pid_inelastic_h       = new TH2F * [n_each_h2D];
    n_sd_minus_bins_mctruthLoose_pid_inelastic_veto_h = new TH2F * [n_each_h2D];
    n_sd_plus_bins_mctruthLoose_pid_inelastic_veto_h  = new TH2F * [n_each_h2D];
    n_sd_minus_bins_reco_pid_inelastic_veto_blindFilled_h    = new TH2F * [n_each_h2D];
    n_sd_plus_bins_reco_pid_inelastic_veto_blindFilled_h     = new TH2F * [n_each_h2D];
    n_sd_minus_bins_mctruth_pid_inelastic_veto_blindFilled_h = new TH2F * [n_each_h2D];
    n_sd_plus_bins_mctruth_pid_inelastic_veto_blindFilled_h  = new TH2F * [n_each_h2D];
    

    for (unsigned int i = 0; i < n_each_h2D; i++) {
      title1 = "diff_flag_mc_full_reco_central_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; diff_flag_full_MCtruth; diff_flag_central_RECO";
      diff_flag_mc_full_reco_central_h[i] = new TH2F(title1.Data(), title2.Data(),
                                                     pids_count, processID::pid_min - 0.5, processID::pid_max + 0.5,
                                                     pids_count, processID::pid_min - 0.5, processID::pid_max + 0.5);

      title1 = "diff_flag_mc_full_reco_full_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; diff_flag_full_MCtruth; diff_flag_full_RECO";
      diff_flag_mc_full_reco_full_h[i] = new TH2F(title1.Data(), title2.Data(),
                                                  pids_count, processID::pid_min - 0.5, processID::pid_max + 0.5,
                                                  pids_count, processID::pid_min - 0.5, processID::pid_max + 0.5);

      title1 = "diff_flag_mc_total_mc_central_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; diff_flag_total_MCtruth; diff_flag_central_MCtruth";
      diff_flag_mc_total_mc_central_h[i] = new TH2F(title1.Data(), title2.Data(),
                                                    pids_count, processID::pid_min - 0.5, processID::pid_max + 0.5,
                                                    pids_count, processID::pid_min - 0.5, processID::pid_max + 0.5);

      title1 = "diff_flag_mc_full_mc_central_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; diff_flag_full_MCtruth; diff_flag_central_MCtruth";
      diff_flag_mc_full_mc_central_h[i] = new TH2F(title1.Data(), title2.Data(),
                                                   pids_count, processID::pid_min - 0.5, processID::pid_max + 0.5,
                                                   pids_count, processID::pid_min - 0.5, processID::pid_max + 0.5);

      title1 = "n_sd_minus_bins_mcreco_mctruth_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; RECO sd- bins; MCTruth sd- bins";
      n_sd_minus_bins_mcreco_mctruth_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 30, -1, 29);

      title1 = "n_sd_plus_bins_mcreco_mctruth_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; RECO sd+ bins; MCTruth sd+ bins";
      n_sd_plus_bins_mcreco_mctruth_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 30, -1, 29);

      title1 = "n_sd_minus_bins_mcreco_mctruth_emptyHF_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; RECO sd- bins; MCTruth sd- bins";
      n_sd_minus_bins_mcreco_mctruth_emptyHF_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 30, -1, 29);

      title1 = "n_sd_plus_bins_mcreco_mctruth_emptyHF_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; RECO sd+ bins; MCTruth sd+ bins";
      n_sd_plus_bins_mcreco_mctruth_emptyHF_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 30, -1, 29);

      title1 = "n_sd_minus_bins_mcreco_mctruthLoose_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; RECO sd- bins; MCTruth (Loose) sd- bins";
      n_sd_minus_bins_mcreco_mctruthLoose_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 30, -1, 29);

      title1 = "n_sd_plus_bins_mcreco_mctruthLoose_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; RECO sd+ bins; MCTruth (Loose) sd+ bins";
      n_sd_plus_bins_mcreco_mctruthLoose_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 30, -1, 29);

      title1 = "n_sd_minus_bins_mcreco_mctruthLoose_emptyHF_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; RECO sd- bins; MCTruth (Loose) sd- bins";
      n_sd_minus_bins_mcreco_mctruthLoose_emptyHF_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 30, -1, 29);

      title1 = "n_sd_plus_bins_mcreco_mctruthLoose_emptyHF_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; RECO sd+ bins; MCTruth (Loose) sd+ bins";
      n_sd_plus_bins_mcreco_mctruthLoose_emptyHF_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 30, -1, 29);

      title1 = "n_sd_minus_bins_mctruth_mctruthLoose_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; MCTruth sd- bins; MCTruth (Loose) sd- bins";
      n_sd_minus_bins_mctruth_mctruthLoose_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 30, -1, 29);

      title1 = "n_sd_plus_bins_mctruth_mctruthLoose_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; MCTruth sd+ bins; MCTruth (Loose) sd+ bins";
      n_sd_plus_bins_mctruth_mctruthLoose_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 30, -1, 29);

      title1 = "n_sd_minus_bins_mctruth_mctruthLoose_emptyHF_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; MCTruth sd- bins; MCTruth (Loose) sd- bins";
      n_sd_minus_bins_mctruth_mctruthLoose_emptyHF_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 30, -1, 29);

      title1 = "n_sd_plus_bins_mctruth_mctruthLoose_emptyHF_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; MCTruth sd+ bins; MCTruth (Loose) sd+ bins";
      n_sd_plus_bins_mctruth_mctruthLoose_emptyHF_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 30, -1, 29);

      title1 = "xi_mc_p_mc_total_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; #xi_{p}; #xi_{MCtruth_total}";
      xi_mc_p_mc_total_h[i] = new TH2F(title1.Data(), title2.Data(), 200, -9., 1., 200, -9., 1.);

      title1 = "xi_mc_p_reco_full_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; #xi_{p}; #xi_{RECO_full}";
      xi_mc_p_reco_full_h[i] = new TH2F(title1.Data(), title2.Data(), 200, -9., 1., 200, -9., 1.);

      title1 = "xi_mc_total_mc_full_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; #xi_{MCtruth_total}; #xi_{MCtruth_full}";
      xi_mc_total_mc_full_h[i] = new TH2F(title1.Data(), title2.Data(), 200, -9., 1., 200, -9., 1.);

      title1 = "xi_mc_total_reco_full_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; #xi_{MCtruth_total}; #xi_{RECO_full}";
      xi_mc_total_reco_full_h[i] = new TH2F(title1.Data(), title2.Data(), 200, -9., 1., 200, -9., 1.);

      title1 = "xi_calo_mc_reco_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += "; MCtruth; RECO";
      xi_calo_mc_reco_h[i] = new TH2F(title1.Data(), title2.Data(), 200, -9., 1., 200, -9., 1.);

      title1 = "xi_pf_mc_reco_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += "; MCtruth; RECO";
      xi_pf_mc_reco_h[i] = new TH2F(title1.Data(), title2.Data(), 200, -9., 1., 200, -9., 1.);

      title1 = "xi_cas_mc_reco_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += "; MCtruth; RECO";
      xi_cas_mc_reco_h[i] = new TH2F(title1.Data(), title2.Data(), 200, -9., 1., 200, -9., 1.);

      title1 = "xi_zdc_mc_reco_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += "; MCtruth; RECO";
      xi_zdc_mc_reco_h[i] = new TH2F(title1.Data(), title2.Data(), 200, -9., 1., 200, -9., 1.);

      title1 = "n_sd_minus_bins_mctruth_pid_inelastic_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_minus_bins (MCtruth); ProcessID";
      n_sd_minus_bins_mctruth_pid_inelastic_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

      title1 = "n_sd_plus_bins_mctruth_pid_inelastic_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_plus_bins (MCtruth); ProcessID";
      n_sd_plus_bins_mctruth_pid_inelastic_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

      title1 = "n_sd_minus_bins_reco_pid_inelastic_veto_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_minus_bins (RECO) (with empty HF); ProcessID";
      n_sd_minus_bins_reco_pid_inelastic_veto_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

      title1 = "n_sd_plus_bins_reco_pid_inelastic_veto_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_plus_bins (RECO) (with empty HF); ProcessID";
      n_sd_plus_bins_reco_pid_inelastic_veto_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

      title1 = "n_sd_minus_bins_mctruth_pid_inelastic_veto_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_minus_bins (MCtruth) (with empty HF); ProcessID";
      n_sd_minus_bins_mctruth_pid_inelastic_veto_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

      title1 = "n_sd_plus_bins_mctruth_pid_inelastic_veto_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_plus_bins (MCtruth) (with empty HF); ProcessID";
      n_sd_plus_bins_mctruth_pid_inelastic_veto_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

      title1 = "n_sd_minus_bins_mctruthLoose_pid_inelastic_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_minus_bins (MCtruth[Loose]); ProcessID";
      n_sd_minus_bins_mctruthLoose_pid_inelastic_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

      title1 = "n_sd_plus_bins_mctruthLoose_pid_inelastic_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_plus_bins (MCtruth[Loose]); ProcessID";
      n_sd_plus_bins_mctruthLoose_pid_inelastic_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

      title1 = "n_sd_minus_bins_mctruthLoose_pid_inelastic_veto_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_minus_bins (MCtruth[Loose]) (with empty HF); ProcessID";
      n_sd_minus_bins_mctruthLoose_pid_inelastic_veto_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

      title1 = "n_sd_plus_bins_mctruthLoose_pid_inelastic_veto_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_plus_bins (MCtruth[Loose]) (with empty HF); ProcessID";
      n_sd_plus_bins_mctruthLoose_pid_inelastic_veto_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID
      /******************************************************************************/
      title1 = "n_sd_minus_bins_reco_pid_inelastic_veto_blindFilled_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_minus_bins (RECO) (with empty HF); ProcessID";
      n_sd_minus_bins_reco_pid_inelastic_veto_blindFilled_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

      title1 = "n_sd_plus_bins_reco_pid_inelastic_veto_blindFilled_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_plus_bins (RECO) (with empty HF); ProcessID";
      n_sd_plus_bins_reco_pid_inelastic_veto_blindFilled_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

      title1 = "n_sd_minus_bins_mctruth_pid_inelastic_veto_blindFilled_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_minus_bins (MCtruth) (with empty HF); ProcessID";
      n_sd_minus_bins_mctruth_pid_inelastic_veto_blindFilled_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

      title1 = "n_sd_plus_bins_mctruth_pid_inelastic_veto_blindFilled_h["; title1 += i; title1 += "]";
      title2 = title1; title2 += " ; n_sd_plus_bins (MCtruth) (with empty HF); ProcessID";
      n_sd_plus_bins_mctruth_pid_inelastic_veto_blindFilled_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID
    };

    h2D->push_back(diff_flag_mc_full_reco_central_h);
    h2D->push_back(diff_flag_mc_full_reco_full_h);
    h2D->push_back(diff_flag_mc_full_mc_central_h);
    h2D->push_back(diff_flag_mc_total_mc_central_h);
    h2D->push_back(n_sd_minus_bins_mcreco_mctruth_h);
    h2D->push_back(n_sd_plus_bins_mcreco_mctruth_h);
    h2D->push_back(n_sd_minus_bins_mcreco_mctruth_emptyHF_h);
    h2D->push_back(n_sd_plus_bins_mcreco_mctruth_emptyHF_h);
    h2D->push_back(n_sd_minus_bins_mcreco_mctruthLoose_h);
    h2D->push_back(n_sd_plus_bins_mcreco_mctruthLoose_h);
    h2D->push_back(n_sd_minus_bins_mcreco_mctruthLoose_emptyHF_h);
    h2D->push_back(n_sd_plus_bins_mcreco_mctruthLoose_emptyHF_h);
    h2D->push_back(n_sd_minus_bins_mctruth_mctruthLoose_h);
    h2D->push_back(n_sd_plus_bins_mctruth_mctruthLoose_h);
    h2D->push_back(n_sd_minus_bins_mctruth_mctruthLoose_emptyHF_h);
    h2D->push_back(n_sd_plus_bins_mctruth_mctruthLoose_emptyHF_h);
    h2D->push_back(xi_mc_p_mc_total_h);
    h2D->push_back(xi_mc_p_reco_full_h);
    h2D->push_back(xi_mc_total_mc_full_h);
    h2D->push_back(xi_mc_total_reco_full_h);
    h2D->push_back(xi_calo_mc_reco_h);
    h2D->push_back(xi_pf_mc_reco_h);
    h2D->push_back(xi_cas_mc_reco_h);
    h2D->push_back(xi_zdc_mc_reco_h);
    h2D->push_back(n_sd_minus_bins_mctruth_pid_inelastic_h);
    h2D->push_back(n_sd_plus_bins_mctruth_pid_inelastic_h);
    h2D->push_back(n_sd_minus_bins_reco_pid_inelastic_veto_h);
    h2D->push_back(n_sd_plus_bins_reco_pid_inelastic_veto_h);
    h2D->push_back(n_sd_minus_bins_mctruth_pid_inelastic_veto_h);
    h2D->push_back(n_sd_plus_bins_mctruth_pid_inelastic_veto_h);
    h2D->push_back(n_sd_minus_bins_mctruthLoose_pid_inelastic_h);
    h2D->push_back(n_sd_plus_bins_mctruthLoose_pid_inelastic_h);
    h2D->push_back(n_sd_minus_bins_mctruthLoose_pid_inelastic_veto_h);
    h2D->push_back(n_sd_plus_bins_mctruthLoose_pid_inelastic_veto_h);
    h2D->push_back(n_sd_minus_bins_reco_pid_inelastic_veto_blindFilled_h);
    h2D->push_back(n_sd_plus_bins_reco_pid_inelastic_veto_blindFilled_h);
    h2D->push_back(n_sd_minus_bins_mctruth_pid_inelastic_veto_blindFilled_h);
    h2D->push_back(n_sd_plus_bins_mctruth_pid_inelastic_veto_blindFilled_h);
  }; // end if mc *******************************************************************


  xi_p_reco_full_h       = new TH2F * [n_each_h2D];
  zdcM_vs_castor_h       = new TH2F * [n_each_h2D];
  zdcM_vs_T2primM_h      = new TH2F * [n_each_h2D];

  ZDCm_vs_xiRP_h      = new TH2F * [n_each_h2D];
  ZDCp_vs_xiRP_h      = new TH2F * [n_each_h2D];
  FSCmSi8_vs_xiRP_h   = new TH2F * [n_each_h2D];
  FSCmN_vs_xiRP_h     = new TH2F * [n_each_h2D];
  FSCmN_vs_castor_h     = new TH2F * [n_each_h2D];

  n_sd_minus_bins_plus_bins_h = new TH2F * [n_each_h2D];

  n_sd_minus_bins_reco_pid_inelastic_h = new TH2F * [n_each_h2D];
  n_sd_plus_bins_reco_pid_inelastic_h  = new TH2F * [n_each_h2D];


  for (unsigned int i = 0; i < n_each_h2D; i++) {
    title1 = "xi_p_reco_full_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ;#xi_{p}; #xi_{RECO_full}";
    xi_p_reco_full_h[i]    = new TH2F(title1.Data(), title2.Data(), 22, -1.1, 1.1, 22, -1.1, 1.1);

    title1 = "zdcM_vs_castor_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += ";CASTOR E [GeV]; ZDC- E [a.u.]";
    zdcM_vs_castor_h[i] = new TH2F(title1.Data(), title2.Data(), 65, -500, 6000, 42, -2000, 40000);

    title1 = "zdcM_vs_T2primM_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += ";N_{tr}(T2-); ZDC- E [a.u.]";
    zdcM_vs_T2primM_h[i] = new TH2F(title1.Data(), title2.Data(), 11, -10., 100., 42, -2000, 40000);

    title1 = "ZDCm_vs_xiRP_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += ";#xi_{p}; ZDC- E [a.u.]";
    ZDCm_vs_xiRP_h[i] = new TH2F(title1.Data(), title2.Data(), 22, -1.1, 1.1, 42, -2000, 40000);

    title1 = "ZDCp_vs_xiRP_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += ";#xi_{p}; ZDC+ E [a.u.]";
    ZDCp_vs_xiRP_h[i] = new TH2F(title1.Data(), title2.Data(), 22, -1.1, 1.1, 42, -20000, 400000);

    title1 = "FSCmSi8_vs_xiRP_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += ";#xi_{p};q6(FSC-) [fC]";
    FSCmSi8_vs_xiRP_h[i] = new TH2F(title1.Data(), title2.Data(), 22, -1.1, 1.1, 2, -3000, 197000);

    title1 = "FSCmN_vs_xiRP_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += ";#xi_{p};number of FSC channels";
    FSCmN_vs_xiRP_h[i] = new TH2F(title1.Data(), title2.Data(), 2200, -1.1, 1.1,  10, -1, 9);

    title1 = "FSCmN_vs_castor_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += ";CASTOR E [GeV];;number of FSC channels";
    FSCmN_vs_castor_h[i] = new TH2F(title1.Data(), title2.Data(), 6500, -500, 6000, 10, -1, 9);

    title1 = "n_sd_minus_bins_plus_bins_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ; RECO sd- bins; RECO sd+ bins";
    n_sd_minus_bins_plus_bins_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 30, -1, 29);

    title1 = "n_sd_minus_bins_reco_pid_inelastic_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ; n_sd_minus_bins (RECO); ProcessID";
    n_sd_minus_bins_reco_pid_inelastic_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

    title1 = "n_sd_plus_bins_reco_pid_inelastic_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ; n_sd_plus_bins (RECO); ProcessID";
    n_sd_plus_bins_reco_pid_inelastic_h[i] = new TH2F(title1.Data(), title2.Data(), 30, -1, 29, 6, -0.5, 5.5); //TODO change to ProcessID

  };
  h2D->push_back(xi_p_reco_full_h);
  h2D->push_back(zdcM_vs_castor_h);
  h2D->push_back(zdcM_vs_T2primM_h);
  h2D->push_back(ZDCm_vs_xiRP_h);
  h2D->push_back(ZDCp_vs_xiRP_h);
  h2D->push_back(FSCmSi8_vs_xiRP_h);
  h2D->push_back(FSCmN_vs_xiRP_h);
  h2D->push_back(FSCmN_vs_castor_h);
  h2D->push_back(n_sd_minus_bins_plus_bins_h);
  h2D->push_back(n_sd_minus_bins_reco_pid_inelastic_h);
  h2D->push_back(n_sd_plus_bins_reco_pid_inelastic_h);

  n_sd_minus_bins_h         = new TH1F * [n_each_h1D];
  n_sd_plus_bins_h          = new TH1F * [n_each_h1D];
  central_activity_h        = new TH1F * [n_each_h1D];
  central_activity_mc_h     = new TH1F * [n_each_h1D];
  xi_reco_full_h            = new TH1F * [n_each_h1D];
  sd_flag_central_reco_h    = new TH1F * [n_each_h1D];
  sd_flag_total_reco_h      = new TH1F * [n_each_h1D];
  for (unsigned int i = 0; i < n_each_h1D; i++) {
    title1 = "n_sd_minus_bins_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ; n_sd_minus_bins";
    n_sd_minus_bins_h[i] = new TH1F(title1.Data(), title2.Data(), 30, -1, 29);

    title1 = "n_sd_plus_bins_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ; n_sd_plus_bins";
    n_sd_plus_bins_h[i] = new TH1F(title1.Data(), title2.Data(), 30, -1, 29);

    title1 = "central_activity_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ; #eta";
    central_activity_h[i] = new TH1F(title1.Data(), title2.Data(), 28, -7, 7);

    title1 = "central_activity_mc_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ; #eta";
    central_activity_mc_h[i] = new TH1F(title1.Data(), title2.Data(), 28, -7, 7);

    title1 = "xi_reco_full_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "#xi_{RECO_full}";
    xi_reco_full_h[i] = new TH1F(title1.Data(), title2.Data(), 100, -7, 0);

    title1 = "sd_flag_central_reco_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ; diff_flag_central_RECO";
    sd_flag_central_reco_h[i] = new TH1F(title1.Data(), title2.Data(), pids_count, processID::pid_min - 0.5, processID::pid_max + 0.5);

    title1 = "sd_flag_total_reco_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ; diff_flag_central_RECO";
    sd_flag_total_reco_h[i] = new TH1F(title1.Data(), title2.Data(), pids_count, processID::pid_min - 0.5, processID::pid_max + 0.5);
  };
  h1D->push_back(n_sd_minus_bins_h);
  h1D->push_back(n_sd_plus_bins_h);
  h1D->push_back(central_activity_h);
  h1D->push_back(central_activity_mc_h);
  h1D->push_back(xi_reco_full_h);
  h1D->push_back(sd_flag_central_reco_h);
  h1D->push_back(sd_flag_total_reco_h);
}
