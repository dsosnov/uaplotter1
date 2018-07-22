/*
 *
 */

#include "uamc.h"
#include "uathresholds.h"
#include "iostream"
#include <map>

ClassImp(uamc)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uamc::uamc(TChain *tree,
           TDirectory *dir,
           const bool cmstotem,
           const bool      pPb,
           const short int MC,
           const short unsigned int Ncuts): uabasecentral(cmstotem, MC, Ncuts, dir), ppb(pPb)
{
  MCthuth = 0;
  if (tree_combined_flag) {
    std::cout << "uamc::uamc, there is no CMS MC in the combined ntuples\n";
  } else {
    tree->SetBranchAddress("genPart",          &MCthuth);
  };
  create_histos();
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uamc::~uamc()
{
  std::cout << "uamc::~uamc(): deleting "
            << h1D->size() << "+" << h2D->size() << " histos" << std::endl;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uamc::FillLastEvent(const short unsigned int cut)
{
  if (cut >= n_cuts) {
    std::cout << "uamc::FillLastEvent: required cut number is larger that possible, do nothing. Please define larger uaforward::n_cut!\n";
    return false;
  };
  double pt2 = protonPt * protonPt;
  proton_pt2_h[cut]->Fill(pt2);
  proton_e_h[cut]->Fill(protonE);
  proton_xi_mc_h[cut]->Fill(protonXi, totalXi);
  proton_pt2_xi_h[cut]->Fill(pt2, protonXi);
  //proton_e_eta_h[cut]->Fill(fabs(protonEta), protonE);
  //proton_e_pt2_h[cut]->Fill(pt2, protonE);

  double energySum_neut[12]; memset(energySum_neut, 0, sizeof(energySum_neut));
  double energySum_charged[12]; memset(energySum_charged, 0, sizeof(energySum_charged));
  double energySum_tracks[12]; memset(energySum_tracks, 0, sizeof(energySum_tracks));
  double n_neut[12]; memset(n_neut, 0, sizeof(n_neut));
  double n_tracks[12]; memset(n_tracks, 0, sizeof(n_tracks));
  double n_charged[12]; memset(n_charged, 0, sizeof(n_charged));
  for (auto part : *MCthuth) {
    auto bin = static_cast<int>(floor((part.eta() + 3.0) / 0.5));
    if(bin < 0 || bin >= 12) continue;
    if(part.charge == 0){
      energySum_neut[bin] += part.E();
      neut_eta_eMC_h[cut]->Fill(part.eta(), part.E());
      n_neut[bin] ++;
    } else {
      energySum_charged[bin] += part.E();
      charged_eta_smallpt_h[cut]->Fill(part.eta(), part.Pt());      
      charged_eta_eMC_h[cut]->Fill(part.eta(), part.E());      
      n_charged[bin] ++;
      if(part.Pt() > TRCK_PT_THR){
        energySum_tracks[bin] += part.E();
        tracks_eta_smallpt_h[cut]->Fill(part.eta(), part.Pt());      
        tracks_eta_eMC_h[cut]->Fill(part.eta(), part.E());      
        n_tracks[bin] ++;
      }
    }
  }
  for (unsigned int i = 0; i < 12; i++) {
    auto etai = i * 0.5 - 3.0 + 0.25;
    tracks_eta_eSum_h[cut]->Fill(etai, energySum_tracks[i]);
    neut_eta_eSum_h[cut]->Fill(etai, energySum_neut[i]);
    charged_eta_eSum_h[cut]->Fill(etai, energySum_charged[i]);
    tracks_eta_n_h[cut]->Fill(etai, n_tracks[i]);
    neut_eta_n_h[cut]->Fill(etai, n_neut[i]);
    charged_eta_n_h[cut]->Fill(etai, n_charged[i]);
  }

  return true;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uamc::PrintEventInfo(const bool detailed)
{
  std::cout << "uamc::PrintEventInfo:\n\t";
  PrintActivity(false);
  if (detailed) {
    std::cout << "\tenergy per rapidity bin:\n";
    for (unsigned int bin = 0; bin < N_ETA_BINS; bin++) {
      std::cout << "\t" << energy[bin];
    }; std::cout << std::endl;
  };
  std::cout << "outer E - : " << outRangeE[0] << "; outer E + " << outRangeE[1] << std::endl;
  std::cout << "\tT2-:" << trksT2[0] << "; T2+:" << trksT2[1] << std::endl;
  std::cout << "\tCastor energy: " << castorE << std::endl;
  std::cout << "\tZDC-: " << (eZDCn[0] + eZDCg[0]) << " (" << eZDCg[0]
            << "); ZDC+: " << (eZDCn[1] + eZDCg[1]) << " (" << eZDCg[1] << ")" << std::endl;
  PrintProtonInfo();

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uamc::PrintProtonInfo()
{
  std::cout << "\tuamc proton: side : " << protonSign << "; xi: " << protonXi
            << "; Pt:" << protonPt << "; E:" << protonE << std::endl;
  std::cout << "\tuamcX: E, Pz, xi: " << totalE << "  " << totalPz << "  " << totalXi << "\n";
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uamc::ProceedEvent(const short unsigned int cut, const bool fill, const bool info)
{
  castorE    = 0;
  castorPz   = 0;

  protonSign = 0;
  protonXi   = 1;
  protonPt   = 0;
  protonE    = 0;
  protonPz   = 0;

  PrepareArrays();
  memset(trksT2,     0, sizeof(trksT2));
  memset(outRangeE,  0, sizeof(outRangeE));
  memset(outRangePz, 0, sizeof(outRangePz));
  memset(eZDCn,      0, sizeof(eZDCn));
  memset(eZDCg,      0, sizeof(eZDCg));
  memset(pzZDCn,     0, sizeof(pzZDCn));
  memset(pzZDCg,     0, sizeof(pzZDCg));
  memset(energyMax,  0, sizeof(energyMax));
  memset(energyH0,   0, sizeof(energyH0));
  memset(energyH0Max,  0, sizeof(energyH0Max));
  memset(ptChargedMax,  0, sizeof(ptChargedMax));
  memset(hfE,        0, sizeof(hfE));
  memset(hfPz,       0, sizeof(hfPz));
  memset(hfE_max,    0, sizeof(hfE_max));
  memset(blindSpot_max, 0, sizeof(blindSpot_max));

  totalE  = 0;
  totalPz = 0;

  unsigned int n = 0;
  bool outer, zdc;

  for (auto part : *MCthuth) {
    double e = part.energy();
    float eta    = part.eta();
    float abseta = fabs(eta);
    outer = false;
    zdc   = false;
    int ind = int(eta > 0);
    //<============================================================= T2 "trigger"
    if ((T2_ABSETA_MIN < abseta) && (abseta < T2_ABSETA_MAX)) {
      if ((part.charge != 0) && part.Pt() > T2_PT_THR) {
        trksT2[ind]++;
      };
    };//========

    if ((CAS_ETA_MIN < eta) && (eta < CAS_ETA_MAX)) { //<=============== Castor
      castorE  += e;
      castorPz += part.Pz();
    };//=========

    if ((HF_ETA_MIN < abseta) && (abseta < HF_ETA_MAX)) { //<=============== HF
      hfE[ind]  += e;
      hfPz[ind] += part.Pz();
      if(e > hfE_max[ind]) hfE_max[ind] = e;
    };//=========

    if((3.0 < abseta) && (abseta < HF_ETA_MIN)){ //<====== Blind Spot
      if(e > blindSpot_max[ind]) blindSpot_max[ind] = e;
    }
    
    //                                             <=============== ZDC
    //if((part.charge==0) && (fabs(part.Px())<ZDC_PXY_THR) && (fabs(part.Py())<ZDC_PXY_THR)){
    if ((part.charge == 0) && (abseta > MIN_ZDC_ETA)) {
      zdc = true;
      if (part.pdgId == 2112) {
        eZDCn[ind] += e;
        pzZDCn[ind] += part.Pz();
      } else {
        eZDCg[ind] += e;
        pzZDCg[ind] += part.Pz();
      };
    };


    bool intact_pcandidate = false;
    if (ppb) {
      intact_pcandidate = (eta < -7);
    } else {
      intact_pcandidate = (eta > 7);
    };

    if ((part.pdgId == 2212) && intact_pcandidate && (e > 1500.)) { //<===== intact proton
      if (eta > 0) {
        protonSign = 1;
      } else {
        protonSign = -1;
      };

      protonXi = (MOMBEAM - fabs(part.Pz())) / MOMBEAM;
      protonPt = part.Pt();
      protonEta = part.Eta();
      protonE  = part.E();
      protonPz = part.Pz();
    } else { //<==================================================== all other particles
      double phi = part.Phi();
      if ( (eta >= -4.363 && eta <= -4.191 && phi >= -M_PI + M_PI / 18 * 8 && phi <= -M_PI + M_PI / 18 * 9) ||
           (eta >=  3.839 && eta <=  4.013 && phi >= -M_PI + M_PI / 18 * 5 && phi <= -M_PI + M_PI / 18 * 6) ){
        continue;
      }
      totalE += e;
      totalPz += part.Pz();

      int bin = find_eta_bin(part.eta());
      if (bin >= 0) {
        energy[bin] += e;
        pz[bin]     += part.Pz();
        pt[bin]     += part.Pt();
        if (e > energyMax[bin]) energyMax[bin] = e;
        if (part.charge != 0 && part.Pt() > ptChargedMax[bin])
          ptChargedMax[bin] = part.Pt();
        if (part.charge == 0 && part.pdgId != 22){
          energyH0[bin] += e;
          if (e > energyH0Max[bin]) energyH0Max[bin] = e;
        }

      } else {
        outer = true;
        outRangeE[ind] += e;
        outRangePz[ind] += part.Pz();
      };
    };

    n++;
    if (info) {
      std::cout << n << "  id: " << part.pdgId << " q: " << part.charge << " stat: " << part.status
                << "\teta: " << part.Eta() << "\tphi: " << part.Phi() << "\tPz: " << part.Pz()
                << "\tPt: " << part.Pt() << "\tE:" << part.energy()
                << "\t" << part.energy() - part.Pz() << "\t" << part.energy() + part.Pz() << "\touter: " << outer << "\tzdc: " << zdc << std::endl;
    };
  }; // end loop

  short int sgn =  1;
  if (ppb)   sgn = -1;
  totalXi = (totalE + sgn * totalPz) / (2 * MOMBEAM); // TBD check it!

  for (unsigned int bin = 0; bin < N_ETA_BINS; bin++){
    if (energy[bin] > 0) activity_loose[bin] = true;
    if (TRCK_ETA_MIN_BIN <= bin && bin <= TRCK_ETA_MAX_BIN){
      if (ptChargedMax[bin] > TRCK_PT_THR) activity_tight[bin] = true;
      if (energy[bin] > THR_PFEN_SUME().at(bin)) activity_tight[bin] = true;
    } else {
      if (energyH0[bin] > THR_PFEN_SUME().at(bin)) activity_tight[bin] = true;
    }
  }
  if (fill)
    FillLastEvent(cut);
  if (info)
    PrintEventInfo(true);
  return true;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uamc::create_histos()
{
  n_each_h1D = n_cuts;
  n_each_h2D = n_cuts;
  proton_pt2_h    = new TH1F * [n_each_h1D];
  proton_e_h      = new TH1F * [n_each_h1D];
  //proton_e_eta_h  = new TH2F * [n_each_h2D];
  //proton_e_pt2_h  = new TH2F * [n_each_h2D];
  proton_pt2_xi_h = new TH2F * [n_each_h2D];
  proton_xi_mc_h  = new TH2F * [n_each_h2D];
  tracks_eta_n_h  = new TH2F*[n_each_h2D];
  tracks_eta_smallpt_h  = new TH2F*[n_each_h2D];
  tracks_eta_eSum_h  = new TH2F*[n_each_h2D];
  tracks_eta_eMC_h  = new TH2F*[n_each_h2D];
  charged_eta_n_h  = new TH2F*[n_each_h2D];
  charged_eta_smallpt_h  = new TH2F*[n_each_h2D];
  charged_eta_eSum_h  = new TH2F*[n_each_h2D];
  charged_eta_eMC_h  = new TH2F*[n_each_h2D];
  neut_eta_eSum_h  = new TH2F*[n_each_h2D];
  neut_eta_eMC_h  = new TH2F*[n_each_h2D];
  neut_eta_n_h  = new TH2F*[n_each_h2D];


  TString title1, title2;
  for (unsigned int i = 0; i < n_cuts; i++) {
    title1 = "proton_pt2_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ;P_{t}^{2} [(GeV/c)^{2}]";
    proton_pt2_h[i] = new TH1F(title1.Data(), title2.Data(), 6100,  -0.1, 6);
    proton_pt2_h[i]->SetDirectory(directory);

    title1 = "proton_e_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ;E [GeV]";
    proton_e_h[i] = new TH1F(title1.Data(), title2.Data(), 5040,  1500, 4200);
    proton_e_h[i]->SetDirectory(directory);

    title1 = "proton_pt2_xi_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ;P_{t}^{2} [(GeV/c)^{2}];#xi_{direct}; ]";
    proton_pt2_xi_h[i] = new TH2F(title1.Data(), title2.Data(), 6100,  -0.1, 6, 102, -0.1, 1.1);
    proton_pt2_xi_h[i]->SetDirectory(directory);

//     title1 = "proton_e_eta_h["; title1+=i; title1+="]";
//     title2 = title1; title2+=" ;|#eta|; E [GeV]";
//     proton_e_eta_h[i] = new TH2F(title1.Data(), title2.Data(), 9,  7, 12, 2520,  1500, 4200);
//     proton_e_eta_h[i]->SetDirectory(directory);

//     title1 = "proton_e_pt2_h["; title1+=i; title1+="]";
//     title2 = title1; title2+=" ;E [GeV]; P_{t}^{2} [(GeV/c)^{2}]";
//     proton_e_pt2_h[i] = new TH2F(title1.Data(), title2.Data(), 5100,  -0.1, 5, 2520,  1500, 4200);
//     proton_e_pt2_h[i]->SetDirectory(directory);

    title1 = "proton_xi_mc_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " ;#xi_{direct}; #xi_{X}]";
    proton_xi_mc_h[i] = new TH2F(title1.Data(), title2.Data(), 102, -0.1, 1.1, 102, -0.1, 1.1);
    proton_xi_mc_h[i]->SetDirectory(directory);

    /****************************************************************/
    title1      = "tracks_eta_smallpt_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; p_T";
    tracks_eta_smallpt_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 100, 0, 1);
    tracks_eta_smallpt_h[i]->SetDirectory(directory);

    title1      = "tracks_eta_n_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; N";
    tracks_eta_n_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 100, 0, 100);
    tracks_eta_n_h[i]->SetDirectory(directory);
    
    title1      = "tracks_eta_eSum_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; E_{Sum}";
    tracks_eta_eSum_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 10000, 0, 100);
    tracks_eta_eSum_h[i]->SetDirectory(directory);

    title1      = "tracks_eta_eMC_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; E_{MCTruth}";
    tracks_eta_eMC_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 1000, 0, 10);
    tracks_eta_eMC_h[i]->SetDirectory(directory);

    title1      = "charged_eta_n_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; N";
    charged_eta_n_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 100, 0, 100);
    charged_eta_n_h[i]->SetDirectory(directory);
    
    title1      = "charged_eta_smallpt_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; p_T";
    charged_eta_smallpt_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 100, 0, 1);
    charged_eta_smallpt_h[i]->SetDirectory(directory);

    title1      = "charged_eta_eSum_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; E_{Sum}";
    charged_eta_eSum_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 10000, 0, 100);
    charged_eta_eSum_h[i]->SetDirectory(directory);

    title1      = "charged_eta_eMC_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; E_{MCTruth}";
    charged_eta_eMC_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 1000, 0, 10);
    charged_eta_eMC_h[i]->SetDirectory(directory);

    title1      = "neut_eta_eSum_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; E_{Sum}";
    neut_eta_eSum_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 10000, 0, 100);
    neut_eta_eSum_h[i]->SetDirectory(directory);

    title1      = "neut_eta_eMC_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; E_{MCTruth}";
    neut_eta_eMC_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 1000, 0, 10);
    neut_eta_eMC_h[i]->SetDirectory(directory);

    title1      = "neut_eta_n_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; N";
    neut_eta_n_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 100, 0, 100);
    neut_eta_n_h[i]->SetDirectory(directory);
  };
  h1D->push_back(proton_pt2_h);
  h1D->push_back(proton_e_h);
  h2D->push_back(proton_pt2_xi_h);
  //h2D->push_back(proton_e_eta_h);
  //h2D->push_back(proton_e_pt2_h);
  h2D->push_back(proton_xi_mc_h);

  h2D->push_back(tracks_eta_smallpt_h);
  h2D->push_back(tracks_eta_eSum_h);
  h2D->push_back(tracks_eta_eMC_h);
  h2D->push_back(tracks_eta_n_h);
  h2D->push_back(charged_eta_n_h);
  h2D->push_back(charged_eta_smallpt_h);
  h2D->push_back(charged_eta_eSum_h);
  h2D->push_back(charged_eta_eMC_h);
  h2D->push_back(neut_eta_eSum_h);
  h2D->push_back(neut_eta_eMC_h);
  h2D->push_back(neut_eta_n_h);
}
