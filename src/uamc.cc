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
uamc::uamc(TChain* tree, 
           TDirectory* dir, 
           const bool cmstotem,
	   const bool      pPb,
           const short int MC, 
           const short unsigned int Ncuts): uabasecentral(cmstotem, MC, Ncuts, dir), ppb(pPb)
{
  MCthuth = 0;
  if(tree_combined_flag){
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
  if( cut>=n_cuts ){
    std::cout << "uamc::FillLastEvent: required cut number is larger that possible, do nothing. Please define larger uaforward::n_cut!\n"; 
    return false;
  }; 
  double pt2 = protonPt*protonPt;
  proton_pt2_h[cut]->Fill(pt2);
  proton_e_h[cut]->Fill(protonE);
  proton_xi_mc_h[cut]->Fill(protonXi, totalXi);
  proton_pt2_xi_h[cut]->Fill(pt2, protonXi);
  //proton_e_eta_h[cut]->Fill(fabs(protonEta), protonE);
  //proton_e_pt2_h[cut]->Fill(pt2, protonE);
  return true;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void uamc::PrintEventInfo(const bool detailed)
{
    std::cout << "uamc::PrintEventInfo:\n\t";
    PrintActivity(false);
    if(detailed){
      std::cout << "\tenergy per rapidity bin:\n";
      for(unsigned int bin=0; bin<N_ETA_BINS; bin++){
	std::cout << "\t"<< energy[bin];
      };std::cout << std::endl;  
    };
    std::cout << "outer E - : " << outRangeE[0] << "; outer E + " << outRangeE[1] << std::endl;
    std::cout << "\tT2-:" << trksT2[0] << "; T2+:" << trksT2[1] << std::endl;
    std::cout << "\tCastor energy: " << castorE << std::endl;
    std::cout << "\tZDC-: " << (eZDCn[0]+eZDCg[0]) << " (" << eZDCg[0] 
	    << "); ZDC+: " << (eZDCn[1]+eZDCg[1]) << " (" << eZDCg[1] << ")"<< std::endl;
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
  memset(trksT2,     0,     sizeof(trksT2));
  memset(outRangeE,  0,     sizeof(outRangeE));
  memset(outRangePz, 0,     sizeof(outRangePz));
  memset(eZDCn,      0,     sizeof(eZDCn));
  memset(eZDCg,      0,     sizeof(eZDCg));
  memset(pzZDCn,     0,     sizeof(pzZDCn));
  memset(pzZDCg,     0,     sizeof(pzZDCg));
  
  totalE  = 0;
  totalPz = 0;
  
  unsigned int n=0;
  bool outer, zdc;
  
  for(std::vector<MyGenPart>::iterator part=MCthuth->begin(); part!=MCthuth->end(); ++part){
    double e = (*part).Energy();
    float eta    = (*part).PseudoRapidity();
    float abseta = fabs(eta);
    outer = false;
    zdc   = false;
    //<============================================================= T2 "trigger"
    if( (abseta>T2_ABSETA_MIN) && (abseta<T2_ABSETA_MAX) ){
      if ( ((*part).charge!=0) && (*part).Pt()>T2_PT_THR ){
	int ind=int(eta>0);
	trksT2[ind]++;
      };
    };//========
    
    if( (eta<CAS_ETA_MAX) && (eta>CAS_ETA_MIN) ){//<=============== Castor 
      castorE  +=e;
      castorPz +=(*part).Pz();
    };//=========

    //                                             <=============== ZDC
    //if(((*part).charge==0) && (fabs((*part).Px())<ZDC_PXY_THR) && (fabs((*part).Py())<ZDC_PXY_THR)){
    if(((*part).charge==0) && (abseta>MIN_ZDC_ETA)){
      zdc = true;
      int ind=int(eta>0);
      if( (*part).pdgId==2112 ){
	eZDCn[ind]+=e; 
	pzZDCn[ind]+=(*part).Pz();
      }else{
	eZDCg[ind]+=e; 
	pzZDCg[ind]+=(*part).Pz();
      };
    };  
    
    
    bool intact_pcandidate = false;
    if(ppb){
      intact_pcandidate = (eta<-7);
    }else{
      intact_pcandidate = (eta>7);
    };
    
    if( ((*part).pdgId==2212) && intact_pcandidate && (e>1500.) ){ //<===== intact proton
      if(eta>0){
	protonSign = 1;
      }else{
	protonSign = -1;
      };
    
      protonXi = ( MOMBEAM - fabs((*part).Pz()) ) / MOMBEAM;
      protonPt = (*part).Pt();
      protonEta= (*part).Eta();
      protonE  = (*part).E();
      protonPz = (*part).Pz();
    }else{ //<==================================================== all other particles                                              
      totalE +=e;
      totalPz+=(*part).Pz();
      
      int bin = find_eta_bin((*part).PseudoRapidity());
      if(bin>=0){
	energy[bin] += e;
	pz[bin]     += (*part).Pz();
	pt[bin]     += (*part).Pt();
      }else{
	outer = true;
	int ind = int(eta>0);
	outRangeE[ind]+=e;
	outRangePz[ind]+=(*part).Pz();
      };
    };
    
    n++;
    if(info){
      std::cout << n << "  id: " << (*part).pdgId << " q: " << (*part).charge << " stat: " << (*part).status 
		<< "\teta: " << (*part).Eta() << "\tphi: " << (*part).Phi() << "\tPz: " << (*part).Pz() 
		<< "\tPt: " << (*part).Pt() << "\tE:" << (*part).Energy() 
		<< "\t" << (*part).Energy()-(*part).Pz() << "\t" << (*part).Energy()+(*part).Pz() << "\touter: " << outer << "\tzdc: " << zdc << std::endl;
    };
  }; // end loop
  
  short int sgn =  1;
  if(ppb)   sgn = -1;
  totalXi = (totalE + sgn*totalPz)/(2*MOMBEAM); // TBD check it!
  
  for(unsigned int bin=0; bin<N_ETA_BINS; bin++)
    if(energy[bin]>0) {
      activity_loose[bin]=true;
      activity_tight[bin]=true;
    };
  if(fill)
    FillLastEvent(cut);
  if(info)
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
  
  TString title1, title2;
  for(unsigned int i=0; i<n_cuts; i++){
    title1 = "proton_pt2_h["; title1+=i; title1+="]";
    title2 = title1; title2+=" ;P_{t}^{2} [(GeV/c)^{2}]";
    proton_pt2_h[i] = new TH1F(title1.Data(), title2.Data(), 6100,  -0.1, 6);
    proton_pt2_h[i]->SetDirectory(directory); 

    title1 = "proton_e_h["; title1+=i; title1+="]";
    title2 = title1; title2+=" ;E [GeV]";
    proton_e_h[i] = new TH1F(title1.Data(), title2.Data(), 5040,  1500, 4200);
    proton_e_h[i]->SetDirectory(directory); 

    title1 = "proton_pt2_xi_h["; title1+=i; title1+="]";
    title2 = title1; title2+=" ;P_{t}^{2} [(GeV/c)^{2}];#xi_{direct}; ]";
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

    title1 = "proton_xi_mc_h["; title1+=i; title1+="]";
    title2 = title1; title2+=" ;#xi_{direct}; #xi_{X}]";
    proton_xi_mc_h[i] = new TH2F(title1.Data(), title2.Data(), 102, -0.1, 1.1, 102, -0.1, 1.1);
    proton_xi_mc_h[i]->SetDirectory(directory); 
  };
  h1D->push_back(proton_pt2_h);
  h1D->push_back(proton_e_h);
  h2D->push_back(proton_pt2_xi_h);
  //h2D->push_back(proton_e_eta_h);
  //h2D->push_back(proton_e_pt2_h);
  h2D->push_back(proton_xi_mc_h);
}
