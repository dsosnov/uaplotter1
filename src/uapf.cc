/*
 *
 */

#include "uapf.h"

#include "iostream"

ClassImp(uapf)

/*!
 * \li 'true'  — for calculating activity maximum element at bin used;
 * \li 'false' — for calculating activity sum energy from bin used;
 */
const bool PF_ACTIVITY_USE_TOWER_ENERGY = false;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uapf::uapf(TChain *tree,
           TDirectory *dir,
           const bool cmstotem,
           const short int MC,
           const short unsigned int Ncuts): uabasecentral(cmstotem, MC, Ncuts, dir)
{
  PFCand = 0;
  if (tree_combined_flag) {
    tree->SetBranchAddress("cmsParticleFlowUA", &PFCand);
  } else {
    tree->SetBranchAddress("particleFlow",      &PFCand);
  };
  create_histos();
  uabasecentral::create_histos();
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uapf::~uapf()
{
  std::cout << "uapf::~uapf(): deleting "
            << h1D->size() << "+" << h2D->size() << " histos" << std::endl;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uapf::FillLastEvent(const short unsigned int cut)
{
  if (cut >= n_cuts) {
    std::cout << "uapf::FillLastEvent: required cut number is larger that possible, do nothing. Please define larger uaforward::n_cut!\n";
    return false;
  };
  uabasecentral::FillLastEvent(cut);
  pfcand_h[cut]->Fill(nPfCand);
  for (unsigned int bin = 0; bin < N_ETA_BINS; bin++) {
    if ((bin >= BINMIN) && (bin <= BINMAX)) {
      pf_e_eta_h[cut]->Fill(ETA_BIN_L[bin], energy[bin]);
      pf_eMaxTower_eta_h[cut]->Fill(ETA_BIN_L[bin], energyMax[bin]);
      //pf_pt_eta_h[cut]->Fill(ETA_BIN_L[bin], pt[bin]);
      pf_e_em0_h[cut]->Fill(ETA_BIN_L[bin], energyEM0[bin]);
      pf_eMaxTower_em0_h[cut]->Fill(ETA_BIN_L[bin], energyEM0Max[bin]);
      pf_e_h0_h[cut]->Fill(ETA_BIN_L[bin], energyH0[bin]);
      pf_eMaxTower_h0_h[cut]->Fill(ETA_BIN_L[bin], energyH0Max[bin]);
      pf_e_hch_h[cut]->Fill(ETA_BIN_L[bin], energyHch[bin]);
      pf_eMaxTower_hch_h[cut]->Fill(ETA_BIN_L[bin], energyHchMax[bin]);
    };
  };
  double energySum_pfneut[12]; memset(energySum_pfneut, 0, sizeof(energySum_pfneut));
  double energySum_pfcharged[12]; memset(energySum_pfcharged, 0, sizeof(energySum_pfcharged));
  double energySum_pftracks[12]; memset(energySum_pftracks, 0, sizeof(energySum_pftracks));
  int n_pfneut[12]; memset(n_pfneut, 0, sizeof(n_pfneut));
  int n_pfcharged[12]; memset(n_pfcharged, 0, sizeof(n_pfcharged));
  int n_pftracks[12]; memset(n_pftracks, 0, sizeof(n_pftracks));
  for (auto pf: *PFCand) {
    auto bin = static_cast<int>(floor((pf.eta() + 3.0) / 0.5));
    if(bin < 0 || bin >= 12) continue;
    if(pf.charge == 0){
      energySum_pfneut[bin] += pf.E();
      pfneut_eta_ePF_h[cut]->Fill(pf.eta(), pf.E());
      n_pfneut[bin]++;
    } else {
      energySum_pfcharged[bin] += pf.E();
      pfcharged_eta_smallpt_h[cut]->Fill(pf.eta(), pf.Pt());      
      pfcharged_eta_ePF_h[cut]->Fill(pf.eta(), pf.E());      
      n_pfcharged[bin]++;;
      if(pf.Pt() > TRCK_PT_THR){
        energySum_pftracks[bin] += pf.E();
        pftracks_eta_smallpt_h[cut]->Fill(pf.eta(), pf.Pt());      
        n_pftracks[bin]++;;
      }
    }
  }
  for (unsigned int i = 0; i < 12; i++) {
    auto etai = i * 0.5 - 3.0 + 0.25;    
    pfcharged_eta_eSum_h[cut]->Fill(etai, energySum_pfcharged[i]);
    pfneut_eta_eSum_h[cut]->Fill(etai, energySum_pfneut[i]);
    pfcharged_eta_n_h[cut]->Fill(etai, n_pfcharged[i]);
    pfneut_eta_n_h[cut]->Fill(etai, n_pfneut[i]);
    pftracks_eta_eSum_h[cut]->Fill(etai, energySum_pftracks[i]);
    pftracks_eta_n_h[cut]->Fill(etai, n_pftracks[i]);
  }
  return true;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uapf::PrintEventInfo(const bool detailed)
{
  std::cout << "uapf::PrintEventInfo: total " << nPfCand << " candidates\n\t";
  PrintActivity(false);
  if (detailed) {
    std::cout << "\tenergy per bin:\n\t";
    for (unsigned short int bin = 0; bin < N_ETA_BINS; bin++) {
      std::cout << energy[bin] << "\t";
    }; std::cout << std::endl;
  };
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uapf::ProceedEvent(const short unsigned int cut, const bool fill, const bool info)
{
  PrepareArrays();
  nPfCand = 0;
  memset(energyEM0,  0, sizeof(energyEM0));
  memset(energyH0, 0, sizeof(energyH0));
  memset(energyHch, 0, sizeof(energyHch));
  memset(energyMax,  0, sizeof(energyMax));
  memset(energyEM0Max,  0, sizeof(energyEM0Max));
  memset(energyH0Max, 0, sizeof(energyH0Max));
  memset(energyHchMax, 0, sizeof(energyHchMax));
  for (auto pf: *PFCand) {
    nPfCand++;
    if (info)
      std::cout << "\tpf#" << nPfCand << "  " << pf.particleId << " eta: " << pf.Eta() << "\tphi: " << pf.Phi() << "\tPt: " << pf.Pt() << "\tE:" << pf.energy() << std::endl;
    int bin = find_eta_bin(pf.eta());
    if (bin >= 0) {
      energy[bin] += pf.energy();
      if (pf.energy() > energyMax[bin]) energyMax[bin] = pf.energy();
      pt[bin] += pf.Pt();
      pz[bin] += pf.Pz();
      if ((pf.particleId == 5) || (pf.particleId == 6) || (pf.particleId == 0)) {// h0 and h_HF and unknown
        energyH0[bin] += pf.energy();
        if (pf.energy() > energyH0Max[bin]) energyH0Max[bin] = pf.energy();
      } else if (pf.particleId == 1) {                 // charged hadron
        energyHch[bin] += pf.energy();
        if (pf.energy() > energyHchMax[bin]) energyHchMax[bin] = pf.energy();
      } else if ((pf.particleId == 4) || (pf.particleId == 7)) {// gamma and egamma_HF
        energyEM0[bin] += pf.energy();
        if (pf.energy() > energyEM0Max[bin]) energyEM0Max[bin] = pf.energy();
      } else if ((pf.particleId == 2) || (pf.particleId == 3))  {
        if (info)
          std::cout << "lepton???" << std::endl
                    << "\tpf#" << nPfCand << "  " << pf.particleId << " eta: " << pf.Eta() << "\tphi: " << pf.Phi() << "\tPt: " << pf.Pt() << "\tE:" << pf.energy() << std::endl;
      }
    }; // end bin >=0
  }; // end loop
  for (unsigned int bin = 0; bin < N_ETA_BINS; bin++) {
    auto thr = (PF_ACTIVITY_USE_TOWER_ENERGY) ? THR_PFEN_TWRE().at(bin) : THR_PFEN_SUME().at(bin); // Current threshold
    double e_t, e_l;
    if(bin < TRCK_ETA_MIN_BIN || bin > TRCK_ETA_MAX_BIN){
      e_l = (PF_ACTIVITY_USE_TOWER_ENERGY) ? energyH0Max[bin]          : energyH0[bin];             // Energy for comparing with threshold for 'loose activity'
      e_t = (PF_ACTIVITY_USE_TOWER_ENERGY) ? energyH0Max[bin]          : energyH0[bin];             // Energy for comparing with threshold for 'tight activity'
    } else {
      e_l = (PF_ACTIVITY_USE_TOWER_ENERGY) ? energyMax[bin]          : energy[bin];             // Energy for comparing with threshold for 'loose activity'
      e_t = (PF_ACTIVITY_USE_TOWER_ENERGY) ? energyMax[bin]          : energy[bin];             // Energy for comparing with threshold for 'tight activity'
    }
    if(thr > 0){
      if(e_l > thr) activity_loose[bin] = true;
      if(e_t > thr) activity_tight[bin] = true;
    }
  }
  if (info)
    PrintEventInfo(true);
  if (fill)
    FillLastEvent(cut);
  return true;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uapf::create_histos()
{
  TString title1, title2;
  n_each_h1D  = n_cuts;
  n_each_h2D  = n_cuts;
  pfcand_h    = new TH1F*[n_each_h1D];
  pf_e_eta_h  = new TH2F*[n_each_h2D];
  pf_eMaxTower_eta_h  = new TH2F*[n_each_h2D];
  pfcharged_eta_n_h  = new TH2F*[n_each_h2D];
  pfcharged_eta_smallpt_h  = new TH2F*[n_each_h2D];
  pfcharged_eta_eSum_h  = new TH2F*[n_each_h2D];
  pfcharged_eta_ePF_h  = new TH2F*[n_each_h2D];
  pfneut_eta_n_h  = new TH2F*[n_each_h2D];
  pfneut_eta_eSum_h  = new TH2F*[n_each_h2D];
  pfneut_eta_ePF_h  = new TH2F*[n_each_h2D];
  pftracks_eta_n_h  = new TH2F*[n_each_h2D];
  pftracks_eta_smallpt_h  = new TH2F*[n_each_h2D];
  pftracks_eta_eSum_h  = new TH2F*[n_each_h2D];
  //pf_pt_eta_h = new TH2F*[n_each_h2D];
  // neutral
  pf_e_h0_h  = new TH2F*[n_each_h2D];
  pf_eMaxTower_h0_h  = new TH2F*[n_each_h2D];
  pf_e_hch_h = new TH2F*[n_each_h2D];
  pf_eMaxTower_hch_h = new TH2F*[n_each_h2D];
  pf_e_em0_h = new TH2F*[n_each_h2D];
  pf_eMaxTower_em0_h = new TH2F*[n_each_h2D];

  for (unsigned int i = 0; i < n_cuts; i++) {
    title1 = "pfcand_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; N_{pfcand}";
    pfcand_h[i] = new TH1F(title1.Data(), title2.Data(), 201,  -1, 200);
    pfcand_h[i]->SetDirectory(directory);

    title1 = "pf_e_eta_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; #eta;E [GeV]";
    pf_e_eta_h[i] = new TH2F(title1.Data(), title2.Data(), 28, -7, 7, 21000, -100, 2000);
    pf_e_eta_h[i]->SetDirectory(directory);

    title1 = "pf_eMaxTower_eta_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; #eta;E_{Tower} [GeV]";
    pf_eMaxTower_eta_h[i] = new TH2F(title1.Data(), title2.Data(), 28, -7, 7, 21000, -100, 2000);
    pf_eMaxTower_eta_h[i]->SetDirectory(directory);
    /*
    title1 = "pf_pt_eta_h["; title1+=i; title1+="]";
    title2 = title1; title2+="; #eta; P_{t} [GeV/c]";
    pf_pt_eta_h[i] = new TH2F(title1.Data(), title2.Data(), 28, -7,7, 21000, -100, 2000);
    pf_pt_eta_h[i]->SetDirectory(directory);
    */
    title1 = "pf_e_h0_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; #eta;E [GeV]";
    pf_e_h0_h[i] = new TH2F(title1.Data(), title2.Data(), 28, -7, 7, 21000, -100, 2000);
    pf_e_h0_h[i]->SetDirectory(directory);

    title1 = "pf_eMaxTower_h0_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; #eta;E_{Tower} [GeV]";
    pf_eMaxTower_h0_h[i] = new TH2F(title1.Data(), title2.Data(), 28, -7, 7, 21000, -100, 2000);
    pf_eMaxTower_h0_h[i]->SetDirectory(directory);

    title1 = "pf_e_hch_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; #eta;E [GeV]";
    pf_e_hch_h[i] = new TH2F(title1.Data(), title2.Data(), 28, -7, 7, 21000, -100, 2000);
    pf_e_hch_h[i]->SetDirectory(directory);

    title1 = "pf_eMaxTower_hch_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; #eta;E_{Tower} [GeV]";
    pf_eMaxTower_hch_h[i] = new TH2F(title1.Data(), title2.Data(), 28, -7, 7, 21000, -100, 2000);
    pf_eMaxTower_hch_h[i]->SetDirectory(directory);

    title1 = "pf_e_em0_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; #eta;E [GeV]";
    pf_e_em0_h[i] = new TH2F(title1.Data(), title2.Data(), 28, -7, 7, 21000, -100, 2000);
    pf_e_em0_h[i]->SetDirectory(directory);

    title1 = "pf_eMaxTower_em0_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; #eta;E_{Tower} [GeV]";
    pf_eMaxTower_em0_h[i] = new TH2F(title1.Data(), title2.Data(), 28, -7, 7, 21000, -100, 2000);
    pf_eMaxTower_em0_h[i]->SetDirectory(directory);

    /****************************************************************/
    title1      = "pfcharged_eta_smallpt_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; p_T";
    pfcharged_eta_smallpt_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 100, 0, 1);
    pfcharged_eta_smallpt_h[i]->SetDirectory(directory);

    title1      = "pfcharged_eta_n_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; N";
    pfcharged_eta_n_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 100, 0, 100);
    pfcharged_eta_n_h[i]->SetDirectory(directory);

    title1      = "pfcharged_eta_eSum_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; E_{Sum}";
    pfcharged_eta_eSum_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 10000, 0, 100);
    pfcharged_eta_eSum_h[i]->SetDirectory(directory);

    title1      = "pfcharged_eta_ePF_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; E_{PFObject}";
    pfcharged_eta_ePF_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 1000, 0, 10);
    pfcharged_eta_ePF_h[i]->SetDirectory(directory);

    title1      = "pfneut_eta_n_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; N";
    pfneut_eta_n_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 100, 0, 100);
    pfneut_eta_n_h[i]->SetDirectory(directory);

    title1      = "pfneut_eta_eSum_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; E_{Sum}";
    pfneut_eta_eSum_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 10000, 0, 100);
    pfneut_eta_eSum_h[i]->SetDirectory(directory);

    title1      = "pfneut_eta_ePF_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; E_{PFObject}";
    pfneut_eta_ePF_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 1000, 0, 10);
    pfneut_eta_ePF_h[i]->SetDirectory(directory);

    title1      = "pftracks_eta_smallpt_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; p_T";
    pftracks_eta_smallpt_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 100, 0, 1);
    pftracks_eta_smallpt_h[i]->SetDirectory(directory);

    title1      = "pftracks_eta_n_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; N";
    pftracks_eta_n_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 100, 0, 100);
    pftracks_eta_n_h[i]->SetDirectory(directory);

    title1      = "pftracks_eta_eSum_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; E_{Sum}";
    pftracks_eta_eSum_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 10000, 0, 100);
    pftracks_eta_eSum_h[i]->SetDirectory(directory);

  };
  h1D->push_back(pfcand_h);
  h2D->push_back(pf_e_eta_h);
  h2D->push_back(pf_eMaxTower_eta_h);
  //h2D->push_back(pf_pt_eta_h);
  h2D->push_back(pf_e_h0_h);
  h2D->push_back(pf_eMaxTower_h0_h);
  h2D->push_back(pf_e_hch_h);
  h2D->push_back(pf_eMaxTower_hch_h);
  h2D->push_back(pf_e_em0_h);
  h2D->push_back(pf_eMaxTower_em0_h);
  h2D->push_back(pfcharged_eta_smallpt_h);
  h2D->push_back(pfcharged_eta_eSum_h);
  h2D->push_back(pfcharged_eta_ePF_h);
  h2D->push_back(pfneut_eta_eSum_h);
  h2D->push_back(pfneut_eta_ePF_h);
  h2D->push_back(pfcharged_eta_n_h);
  h2D->push_back(pfneut_eta_n_h);
  h2D->push_back(pftracks_eta_n_h);
  h2D->push_back(pftracks_eta_smallpt_h);
  h2D->push_back(pftracks_eta_eSum_h);
}
//
