/*
 *
 */

#include "uacalo.h"

#include <map>
#include "iostream"

ClassImp(uacalo)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uacalo::uacalo(TChain *tree,
               TDirectory *dir,
               const bool cmstotem,
               const short int MC,
               const short unsigned int Ncuts): uabasecentral(cmstotem, MC, Ncuts, dir)
{
  HF     = 0;
  Towers = 0;
  if (tree_combined_flag) {
    tree->SetBranchAddress("cmshfRecHitsUA",   &HF);
    tree->SetBranchAddress("cmsCaloTowersUA",  &Towers);
  } else {
    tree->SetBranchAddress("hfRecHits",   &HF);
    tree->SetBranchAddress("caloTowers",  &Towers);
  };

  for(uint side = 0; side<2; side++)
    hf_max_energy_tower[side] = (double*) malloc(nBinsHF*sizeof(double));

  create_histos();
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uacalo::~uacalo()
{
  std::cout << "uacalo::~uacalo(): deleting "
            << h1D->size() << "+" << h2D->size() << " histos" << std::endl;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uacalo::FillLastEvent(const short unsigned int cut)
{
  if (cut >= n_cuts) {
    std::cout << "uatracking::FillLastEvent: required cut number is larger that possible, do nothing. Please define larger uaforward::n_cut!\n";
    return false;
  };
  hf_etotal_minus_h[cut]->Fill(hf_total_energy_tower[0]);
  hf_etotal_plus_h [cut]->Fill(hf_total_energy_tower[1]);
  for (unsigned int bin = 0; bin < N_ETA_BINS; bin++)
    if ((bin >= BINMIN) && (bin <= BINMAX)) {
      calotower_e_eta_h[cut]->Fill(ETA_BIN_L[bin], energy[bin]);
      calotower_eMaxTower_eta_h[cut]->Fill(ETA_BIN_L[bin], energyMax[bin]);
     }

  double MaxHFtow[2] = {0, 0};
  for (short unsigned int side = 0; side < 2; side++)
    for (short unsigned int indx = 0; indx < nBinsHF; indx++) {
      if (MaxHFtow[side] < hf_max_energy_tower[side][indx])
        MaxHFtow[side] = hf_max_energy_tower[side][indx];
    };
  hf_max_towerE_minus_h[cut]->Fill(MaxHFtow[0]);
  hf_max_towerE_plus_h[cut]->Fill(MaxHFtow[1]);
  if (hf_max_energy_tower_ill[0]>0.1) hf_max_towerE_minus_ill_h[cut]->Fill(hf_max_energy_tower_ill[0]);
  if (hf_max_energy_tower_ill[1]>0.1) hf_max_towerE_plus_ill_h[cut]->Fill(hf_max_energy_tower_ill[1]);

  for (auto t: *Towers) {
    calotower_e_eta_phi_h[cut]->Fill(t.eta(), t.phi(), t.energy());
    auto ind = int(t.zside > 0);
    switch (ind) {
      case 0:
        calotower_e_eta_phi_minus_h[cut]->Fill(t.eta(), t.phi(), t.energy());
        break;
      case 1:
        calotower_e_eta_phi_plus_h[cut]->Fill(t.eta(), t.phi(), t.energy());
        break;
    };
  }

  return true;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uacalo::PrintEventInfo(const bool detailed)
{
  std::cout << "uacalo::PrintEventInfo:\n\t";
  std::cout << "tower trigger HF- " << hf_trigger_tower[0]
            << "\t; HF+ "             << hf_trigger_tower[1] << std::endl;
  std::cout << "\tenergy HF-         " << hf_total_energy_tower[0]
            << "\t; HF+ "            << hf_total_energy_tower[1] << std::endl;
  PrintActivity(false);
  if (detailed) {
    std::cout << "\tnumber of calotowers per bin:\n\t";
    for (unsigned short int bin = 0; bin < N_ETA_BINS; bin++) {
      std::cout << calotowers[bin] << "  ";
    }; std::cout << std::endl;
    std::cout << "\tenergy per bin:\n\t";
    for (unsigned short int bin = 0; bin < N_ETA_BINS; bin++) {
      std::cout << energy[bin] << "\t";
    }; std::cout << std::endl;
    std::cout <<  "Energies: ";
    for (auto i: energy) {std::cout << i << " ";} std::cout << std::endl;
    std::cout <<  "Energies (max): ";
    for (auto i: energyMax) {std::cout << i << " ";} std::cout << std::endl;
   }
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool uacalo::ProceedEvent(const short unsigned int cut, const bool fill, const bool info)
{

  memset(hf_total_energy_tower,  0, sizeof(hf_total_energy_tower));
  for(uint side = 0; side<2; side++)
    for(uint bin = 0; bin<nBinsHF; bin++)
      hf_max_energy_tower[side][bin] = 0;
  //memset(hf_max_energy_tower,    0, sizeof(hf_max_energy_tower));
  memset(hf_max_energy_tower_ill,    0, sizeof(hf_max_energy_tower_ill));
  memset(hf_trigger_tower,   false, sizeof(hf_trigger_tower));
  memset(calotowers,             0, sizeof(calotowers));
  PrepareArrays();                                          // cleans eta arrays
  memset(energyMax, 0, sizeof(energyMax));

  for (std::vector<MyCaloTower>::iterator it = Towers->begin(); it != Towers->end(); ++it) {
    if (it->Pt() > 0) {
      int bin = find_eta_bin(it->eta());
      if (bin >= 0) {
        energy[bin] += it->energy();
        pz[bin]     += it->Pz();
        pt[bin]     += it->Pt();
        if (it->hasHF) {
          if (it->energy() > HF_TOWER_THR[it->zside > 0])
            calotowers[bin]++;
        } else {
          if (it->energy() > CENT_TOWER_THR)
            calotowers[bin]++;
        };
      }; // end if bin

      //___________________________________________________________________
      if ((it->hasHF) && (fabs(it->eta()) <= CALO_ETA_ACC) && (fabs(it->eta()) > HF_ETA_MIN)) {

        if (it->eta() >= -4.363 && it->eta() <= -4.191 &&
            it->phi() >= -M_PI + M_PI / 18 * 8 && it->phi() <= -M_PI + M_PI / 18 * 9) {
          if (it->energy() > hf_max_energy_tower_ill[0]) hf_max_energy_tower_ill[0] = it->energy();
        } else if (it->eta() >= 3.839 && it->eta() <= 4.013 &&
                   it->phi() >= -M_PI + M_PI / 18 * 5 && it->phi() <= -M_PI + M_PI / 18 * 6) {
          if (it->energy() > hf_max_energy_tower_ill[1]) hf_max_energy_tower_ill[1] = it->energy();
        } else {
          short unsigned int side = int (it->eta() > 0);
          short unsigned int indx = find_eta_bin (CALO_ETA_ACC) - find_eta_bin (fabs (it->eta()));
          if (it->energy() > hf_max_energy_tower[side][indx]) {
            hf_max_energy_tower[side][indx] = it->energy();
          }
        }

        // here HF "trigger" only
        if (fabs(it->eta()) < HF_ETA_MAX) {
          short unsigned ind = int(it->zside > 0);
          hf_total_energy_tower[ind] += it->energy();
          if (!hf_trigger_tower[ind] && (it->energy() > HF_TOWER_THR[ind]))
            hf_trigger_tower[ind] = true;
        };
      }; // end HF
      if (it->energy() > energyMax[bin]) energyMax[bin] = it->energy();
     //___________________________________________________________________
    }
  }// end tower loop

  // <===================================================compare to threshold
  for (unsigned int bin = 0; bin < N_ETA_BINS; bin++)
    if ((THR_CALO[bin] > 0) && (energy[bin] > THR_CALO[bin])) {
      activity_loose[bin] = true;
      activity_tight[bin] = true;
    };

  if (info)
    PrintEventInfo(true);
  if (fill)
    FillLastEvent(cut);

  return true;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uacalo::create_histos()
{
  TString title1, title2;

  n_each_h1D = n_cuts;
  n_each_h2D = n_cuts;

  hf_etotal_minus_h = new TH1F * [n_each_h1D];
  hf_etotal_plus_h  = new TH1F * [n_each_h1D];
  hf_max_towerE_minus_h = new TH1F * [n_each_h1D];
  hf_max_towerE_minus_ill_h = new TH1F * [n_each_h1D];
  hf_max_towerE_plus_h  = new TH1F * [n_each_h1D];
  hf_max_towerE_plus_ill_h  = new TH1F * [n_each_h1D];

  //hf_towers_vs_rechits_h     = new TH2F * [n_each_h2D];
  calotower_e_eta_h          = new TH2F * [n_each_h2D];
  calotower_e_eta_phi_h = new TH2F * [n_each_h2D];
  calotower_e_eta_phi_minus_h = new TH2F * [n_each_h2D];
  calotower_e_eta_phi_plus_h = new TH2F * [n_each_h2D];
  calotower_eMaxTower_eta_h = new TH2F * [n_each_h2D];

  // rechits are not done yet
  for (unsigned int i = 0; i < n_cuts; i++) {
    title1 = "hf_etotal_minus_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; E_{HFtowers} [GeV]";
    hf_etotal_minus_h[i] = new TH1F(title1.Data(), title2.Data(), 71000, -100, 7000);
    hf_etotal_minus_h[i]->SetDirectory(directory);

    title1 = "hf_max_towerE_minus_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; E^{max}_{HFtower} [GeV]";
    hf_max_towerE_minus_h[i] = new TH1F(title1.Data(), title2.Data(), 71000, -100, 7000);
    hf_max_towerE_minus_h[i]->SetDirectory(directory);

    title1 = "hf_max_towerE_minus_ill_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; E^{max}_{HFtower} [GeV]";
    hf_max_towerE_minus_ill_h[i] = new TH1F(title1.Data(), title2.Data(), 71000, -100, 7000);
    hf_max_towerE_minus_ill_h[i]->SetDirectory(directory);

    title1 = "hf_etotal_plus_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; E_{HFtowers} [GeV]";
    hf_etotal_plus_h[i] = new TH1F(title1.Data(), title2.Data(), 71000, -100, 7000);
    hf_etotal_plus_h[i]->SetDirectory(directory);

    title1 = "hf_max_towerE_plus_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; E^{max}_{HFtower} [GeV]";
    hf_max_towerE_plus_h[i] = new TH1F(title1.Data(), title2.Data(), 71000, -100, 7000);
    hf_max_towerE_plus_h[i]->SetDirectory(directory);

    title1 = "hf_max_towerE_plus_ill_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; E^{max}_{HFtower} [GeV]";
    hf_max_towerE_plus_ill_h[i] = new TH1F(title1.Data(), title2.Data(), 71000, -100, 7000);
    hf_max_towerE_plus_ill_h[i]->SetDirectory(directory);

    //     title1 = "hf_towers_vs_rechits_h["; title1+=i; title1+="]";
    //     title2 = title1; title2+="; E_{HFrechits} [GeV]; E_{HFtowers}";
    //     hf_towers_vs_rechits_h[i] = new TH2F(title1.Data(), title2.Data(), 2100, -100, 2000, 2100, -100, 2000);
    //     hf_towers_vs_rechits_h[i]->SetDirectory(directory);

    title1 = "calotower_e_eta_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; #eta; E_{towers}";
    calotower_e_eta_h[i] = new TH2F(title1.Data(), title2.Data(), 28, -7, 7, 21000, -100, 2000);
    calotower_e_eta_h[i]->SetDirectory(directory);

    double hfMinusBins[] = {-5.191,-4.889,-4.716,-4.538,-4.363,-4.191,-4.013,-3.839,-3.664,-3.489,-3.314,-3.139,-2.964,-2.853};
    title1 = "calotower_e_eta_phi_minus_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; #eta; #phi";
    calotower_e_eta_phi_minus_h[i] = new TH2F(title1.Data(), title2.Data(), 13, &hfMinusBins[0], 36, -M_PI, M_PI);
    calotower_e_eta_phi_minus_h[i]->SetDirectory(directory);

    double hfPlusBins[] = {2.853,2.964,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.191};
    title1 = "calotower_e_eta_phi_plus_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; #eta; #phi";
    calotower_e_eta_phi_plus_h[i] = new TH2F(title1.Data(), title2.Data(), 13, &hfPlusBins[0], 36, -M_PI, M_PI);
    calotower_e_eta_phi_plus_h[i]->SetDirectory(directory);

    double hfBins[] = { -5.191,-4.889,-4.716,-4.538,-4.363,-4.191,-4.013,-3.839,-3.664,-3.489,-3.314,-3.139,-2.964,-2.853, -2.650,-2.500,-2.322,-2.172,-2.043,-1.930,-1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.870, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0.000, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.870, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.740, 1.830,1.930,2.043,2.172,2.322,2.500,2.650, 2.853,2.964,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.191 };
    title1 = "calotower_e_eta_phi_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; #eta; #phi";
    calotower_e_eta_phi_h[i] = new TH2F(title1.Data(), title2.Data(), 82, &hfBins[0], 72, -M_PI, M_PI);
    calotower_e_eta_phi_h[i]->SetDirectory(directory);

    title1 = "calotower_eMaxTower_eta_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += "; #eta; E_{Tower}";
    calotower_eMaxTower_eta_h[i] = new TH2F(title1.Data(), title2.Data(), 28, -7, 7, 21000, -100, 2000);
    calotower_eMaxTower_eta_h[i]->SetDirectory(directory);
  };
  h1D->push_back(hf_etotal_minus_h);
  h1D->push_back(hf_etotal_plus_h);
  h1D->push_back(hf_max_towerE_minus_h);
  h1D->push_back(hf_max_towerE_minus_ill_h);
  h1D->push_back(hf_max_towerE_plus_h);
  h1D->push_back(hf_max_towerE_plus_ill_h);
  //h2D->push_back(hf_towers_vs_rechits_h);
  h2D->push_back(calotower_e_eta_h);
  h2D->push_back(calotower_e_eta_phi_h);
  h2D->push_back(calotower_e_eta_phi_minus_h);
  h2D->push_back(calotower_e_eta_phi_plus_h);
  h2D->push_back(calotower_eMaxTower_eta_h);
}
