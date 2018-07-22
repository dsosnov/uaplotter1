/*
 *
 */

#include "uatracking.h"

#include "iostream"

ClassImp(uatracking)


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uatracking::uatracking(TChain *tree,
                       TDirectory *dir,
                       const bool cmstotem,
                       const short int MC,
                       const short unsigned int Ncuts): uabasecentral(cmstotem, MC, Ncuts, dir)
{
  Vertices = 0;
  Tracks   = 0;

  if (tree_combined_flag) {
    tree->SetBranchAddress("cmsVerticesUA",          &Vertices);
    tree->SetBranchAddress("cmsTracksUA",            &Tracks);
  } else {
    tree->SetBranchAddress("offlinePrimaryVertices", &Vertices);
    tree->SetBranchAddress("generalTracks",          &Tracks);
  };
  create_histos();
  uabasecentral::create_histos();
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uatracking::~uatracking()
{
  std::cout << "uatracking::~uatracking(): deleting "
            << h1D->size() << "+" << h2D->size() << " histos" << std::endl;
}

bool isTrackGood(const MyTracks track){
  return (track.quality[2]);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uatracking::FillLastEvent(const short unsigned int cut)
{
  if (cut >= n_cuts) {
    std::cout << "uatracking::FillLastEvent: required cut number is larger that possible, do nothing. Please define larger uaforward::n_cut!\n";
    return false;
  };
  uabasecentral::FillLastEvent(cut);
  tracks_h[cut]->Fill(nTracks);
  tracks_h[cut + n_cuts]->Fill(nTracksGood);
  vertices_h[cut]->Fill(nVtx);
  vertices_h[cut + n_cuts]->Fill(nVtxGoog);
  for (unsigned int bin = 0; bin < N_ETA_BINS; bin++) {
    tracks_eta_n_h[cut]->Fill(ETA_BIN_L[bin], n_tracks_bin_all[bin]);
    tracks_eta_n_h[cut + n_cuts]->Fill(ETA_BIN_L[bin], n_tracks_bin_good[bin]);
  };
  for (auto t: *Tracks) {
    tracks_pt_h[cut]->Fill(t.pt());
    tracks_eta_pt_h[cut]->Fill(t.eta(), t.pt());
    if (!isTrackGood(t)) continue;
    tracks_pt_h[cut+n_cuts]->Fill(t.pt());
    tracks_eta_pt_h[cut+n_cuts]->Fill(t.eta(), t.pt());
  }
  return true;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uatracking::PrintEventInfo(const bool detailed)
{
  std::cout << "uatracking::PrintEventInfo:\n\t";
  std::cout << "Nvtx = " << nVtx    << " ( " << nVtxGoog    << " good)\n\t";
  std::cout << "Ntrk = " << nTracks << " ( " << nTracksGood << " tight)\n";
  PrintActivity(false);

  if (detailed) {
    std::cout << "\ttracks per eta bin:\n";
    for (unsigned short int bin = 0; bin < N_ETA_BINS; bin++) {
      std::cout << "\t" << n_tracks_bin_all[bin];
    };
    std::cout << "\n\t good tracks per eta bin:\n";
    for (unsigned short int bin = 0; bin < N_ETA_BINS; bin++) {
      std::cout << "\t" << n_tracks_bin_good[bin];
    }; std::cout << std::endl;
  };
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool uatracking::ProceedEvent(const short unsigned int cut, const bool fill, const bool info)
{
  memset(n_tracks_bin_all,  0, sizeof(n_tracks_bin_all));
  memset(n_tracks_bin_good, 0, sizeof(n_tracks_bin_good));
  memset(ptMax_bin_all,  0, sizeof(ptMax_bin_all));
  memset(ptMax_bin_good, 0, sizeof(ptMax_bin_good));
  PrepareArrays(); // cleans eta arrays
  nVtx        = 0;
  nVtxGoog    = 0;
  nTracks     = 0;
  nTracksGood = 0;

  short unsigned int i = 0;
  for (auto v: *Vertices) {
    if (v.ndof == 0) continue;
    nVtx++;
    if ((v.ntracks > 1) && (fabs(v.z) < 15) &&
      (fabs(v.x) < 0.2) && (fabs(v.y) < 0.2)) { // good vertices
      nVtxGoog++;
    }
  } // end vertex loop


  for (auto t: *Tracks) {
    int bin = find_eta_bin(t.eta());
    n_tracks_bin_all[bin]++;
    nTracks++;
    if (t.pt() > ptMax_bin_all[bin]) ptMax_bin_all[bin] = t.pt();
    if (!isTrackGood(t)) continue;
    if (t.pt() > ptMax_bin_good[bin]) ptMax_bin_good[bin] = t.pt();
    n_tracks_bin_good[bin]++;
    nTracksGood++;
    if (info) {
      std::cout << "trck# " << i++ << " ";
      for (int b = 0; b < 5; b++) std::cout << t.quality[b];
      std::cout << " " << t.chi2n << " \t\teta: " << t.Eta() << "\tphi: " << t.Phi()*TMath::RadToDeg() << "\tpt: " << t.Pt() << "\tE:" << t.energy()
                << "\t charge:" << t.charge;
      std::cout << std::endl;
    }
  }

  for (unsigned short int bin = TRCK_ETA_MIN_BIN; bin <= TRCK_ETA_MAX_BIN; bin++) {
    if (ptMax_bin_all[bin]  > TRCK_PT_THR) activity_loose[bin] = true;
    if (ptMax_bin_good[bin] > TRCK_PT_THR) activity_tight[bin] = true;
  };
  if (info)
    PrintEventInfo(true);
  if (fill)
    FillLastEvent(cut);

  return true;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uatracking::create_histos()
{
  n_each_h1D = 2 * n_cuts;
  n_each_h2D = 2 * n_cuts;
  TString title1, title2;

  // vertices
  vertices_h   = new TH1F * [2 * n_cuts];
  //vertices_X_h = new TH1F * [2*n_cuts];
  //vertices_Y_h = new TH1F * [2*n_cuts];
  //vertices_Z_h = new TH1F * [2*n_cuts];

  for (unsigned int i = 0; i < n_cuts; i++) {
    title1 = "vertices_h["; title1 += i; title1 += "]";
    title2 = title1; title2 += " all; nvtx";
    vertices_h[i] = new TH1F(title1.Data(), title2.Data(), 10,  -1, 9);
    vertices_h[i]->SetDirectory(directory);
    title1 = "vertices_h["; title1 += (i + n_cuts); title1 += "]";
    title2 = title1; title2 += " good; nvtx";
    vertices_h[i + n_cuts] = new TH1F(title1.Data(), title2.Data(), 10,  -1, 9);
    vertices_h[i + n_cuts]->SetDirectory(directory);

    /*  if these histos are added, FillLastEvent should be redone
        title1 = "vertices_X_h["; title1+=i; title1+="]";
        title2 = title1; title2+=" all; X [cm]";
        vertices_X_h[i] = new TH1F(title1.Data(), title2.Data(), 100,  -1, 1);
        vertices_X_h[i]->SetDirectory(directory);
        title1 = "vertices_X_h["; title1+=(i+n_cuts); title1+="]";
        title2 = title1; title2+=" good; X [cm]";
        vertices_X_h[i+n_cuts] = new TH1F(title1.Data(), title2.Data(), 100,  -1, 1);
        vertices_h[i+n_cuts]->SetDirectory(directory);

        title1 = "vertices_Y_h["; title1+=i; title1+="]";
        title2 = title1; title2+=" all; Y [cm]";
        vertices_Y_h[i] = new TH1F(title1.Data(), title2.Data(), 100,  -1, 1);
        vertices_Y_h[i]->SetDirectory(directory);
        title1 = "vertices_Y_h["; title1+=(i+n_cuts); title1+="]";
        title2 = title1; title2+=" good; Y [cm]";
        vertices_Y_h[i+n_cuts] = new TH1F(title1.Data(), title2.Data(), 100,  -1, 1);
        vertices_Y_h[i+n_cuts]->SetDirectory(directory);

        title1 = "vertices_Z_h["; title1+=i; title1+="]";
        title2 = title1; title2+=" all; Z [cm]";
        vertices_Z_h[i] = new TH1F(title1.Data(), title2.Data(), 400,  -20, 20);
        vertices_Z_h[i]->SetDirectory(directory);
        title1 = "vertices_Z_h["; title1+=(i+n_cuts); title1+="]";
        title2 = title1; title2+=" good; Z [cm]";
        vertices_Z_h[i+n_cuts] = new TH1F(title1.Data(), title2.Data(), 400,  -20, 20);
        vertices_Z_h[i+n_cuts]->SetDirectory(directory);
      */

  };
  h1D->push_back(vertices_h);
  //h1D->push_back(vertices_X_h);
  //h1D->push_back(vertices_Y_h);
  //h1D->push_back(vertices_Z_h);
  /*
  title1 = "vertices_beamX_h";
  title2 = "vertices_beamX_h; X [cm]";
  vertices_beamX_h = new TH1F(title1.Data(), title2.Data(), 100,  -0.2, 0.2);

  title1 = "vertices_beamY_h";
  title2 = "vertices_beamY_h; Y [cm]";
  vertices_beamY_h = new TH1F(title1.Data(), title2.Data(), 100,  -0.2, 0.2);

  title1 = "vertices_beamZ_h";
  title2 = "vertices_beamZ_h; Z [cm]";
  vertices_beamZ_h = new TH1F(title1.Data(), title2.Data(), 400,  -5, 5);
  */

  tracks_h     = new TH1F * [n_each_h1D];
  tracks_pt_h = new TH1F * [n_each_h1D];
  tracks_eta_n_h = new TH2F * [n_each_h2D];
  tracks_eta_pt_h = new TH2F * [n_each_h2D];
  for (unsigned int i = 0; i < n_cuts; i++) {
    title1      = "tracks_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; ntrcks";
    tracks_h[i] = new TH1F(title1.Data(), title2.Data(), 201, -1, 200);
    tracks_h[i]->SetDirectory(directory);

    title1      = "tracks_h["; title1 += (i + n_cuts); title1 += "]";
    title2      = title1; title2 += " good; ntrcks";
    tracks_h[i + n_cuts] = new TH1F(title1.Data(), title2.Data(), 201, -1, 200);
    tracks_h[i + n_cuts]->SetDirectory(directory);

    title1      = "tracks_pt_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; p_T;";
    tracks_pt_h[i] = new TH1F(title1.Data(), title2.Data(), 7001, -1, 7000);
    tracks_pt_h[i]->SetDirectory(directory);

    title1      = "tracks_pt_h["; title1 += (i + n_cuts); title1 += "]";
    title2      = title1; title2 += " good; p_T";
    tracks_pt_h[i + n_cuts] = new TH1F(title1.Data(), title2.Data(), 7001, -1, 7000);
    tracks_pt_h[i + n_cuts]->SetDirectory(directory);

    title1      = "tracks_eta_n_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; N_{all}";
    tracks_eta_n_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 201, -1, 200);
    tracks_eta_n_h[i]->SetDirectory(directory);

    title1      = "tracks_eta_n_h["; title1 += (i + n_cuts); title1 += "]";
    title2      = title1; title2 += " good; #eta; N_{good}";
    tracks_eta_n_h[i + n_cuts] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 201, -1, 200);
    tracks_eta_n_h[i + n_cuts]->SetDirectory(directory);

    title1      = "tracks_eta_pt_h["; title1 += i; title1 += "]";
    title2      = title1; title2 += "; #eta; p_T";
    tracks_eta_pt_h[i] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 70010, -1, 7000);
    tracks_eta_pt_h[i]->SetDirectory(directory);

    title1      = "tracks_eta_pt_h["; title1 += (i + n_cuts); title1 += "]";
    title2      = title1; title2 += " good; #eta; p_T";
    tracks_eta_pt_h[i + n_cuts] = new TH2F(title1.Data(), title2.Data(), 12, -3, 3, 70010, -1, 7000);
    tracks_eta_pt_h[i + n_cuts]->SetDirectory(directory);
  };
  h1D->push_back(tracks_h);
  h1D->push_back(tracks_pt_h);
  h2D->push_back(tracks_eta_n_h);
  h2D->push_back(tracks_eta_pt_h);
}
