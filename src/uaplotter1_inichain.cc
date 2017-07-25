#include "uaplotter1.h"

#include "iostream"
#include "stdlib.h"

// also sets zdc56
TString uaplotter1::initializeChain(){

  TString str="";

  switch (mc){
    // NOISE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    case -2:{
      str="dataPARandomNew";
      zdc56 = true;
      chainTree->Add("/data/HI13/CmsData/rnd2015/UAcms_210885_RND_000.root");
      chainTree->Add("/data/HI13/CmsData/rnd2015/UAcms_210885_RND_001.root");
      chainTree->Add("/data/HI13/CmsData/rnd2015/UAcms_210885_RND_002.root");
      chainTree->Add("/data/HI13/CmsData/rnd2015/UAcms_210885_RND_003.root");
      chainTree->Add("/data/HI13/CmsData/rnd2015/UAcms_210885_RND_004.root");
      chainTree->Add("/data/HI13/CmsData/rnd2015/UAcms_210885_RND_006.root");
      chainTree->Add("/data/HI13/CmsData/rnd2015/UAcms_210885_RND_008.root");
      chainTree->Add("/data/HI13/CmsData/rnd2015/UAcms_210885_RND_020.root");      
      break;
    };
    case -1:{
      str="data_PARandom";
      zdc56 = true;
      chainTree->Add("/home/kkuzn/pA13/noise/210534/UABaseTree_PAMinBiasUPC_reco_210534_HLT_PARandom_v1_HLT_PAL1Tech53_MB_v1_1.root");
      chainTree->Add("/home/kkuzn/pA13/noise/210534/UABaseTree_PAMinBiasUPC_reco_210534_HLT_PARandom_v1_HLT_PAL1Tech53_MB_v1_2.root");
      chainTree->Add("/home/kkuzn/pA13/noise/210534/UABaseTree_PAMinBiasUPC_reco_210534_HLT_PARandom_v1_HLT_PAL1Tech53_MB_v1_3.root");
      break;
    };
    case 0:{ // DATA ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zdc56 = false;
      if(ppb){
	  str="data_pPb";
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8869_10000001_20000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8869_1_10000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8869_20000001_30000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8869_30000001_40000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8870_30000001_40000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8870_40000001_50000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8870_50000001_60000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8870_60000001_70000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8871_100000001_110000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8871_110000001_120000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8871_70000001_80000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8871_80000001_90000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8871_90000001_100000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8872_120000001_130000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8872_130000001_140000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8873_130000001_140000000.root");
	  //chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8873_140000001_150000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8873_150000001_160000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8873_160000001_170000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8873_170000001_180000000.root");
	  //chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8874_170000001_180000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8874_180000001_190000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8875_180000001_190000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8875_190000001_200000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8876_200000001_210000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8876_210000001_220000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8876_220000001_230000000.root");
	  //chainTree->Add("/data/HI13/CmsTotemData/210885/cms_tot_8876_230000001_240000000.root");
      }else{
	  str="data_Pbp";
	  chainTree->Add("/data/HI13/CmsTotemData/211328/cms_tot_8933_100000001_110000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211328/cms_tot_8933_10000001_20000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211328/cms_tot_8933_110000001_122000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211328/cms_tot_8933_20000001_30000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211328/cms_tot_8933_30000001_40000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211328/cms_tot_8933_40000001_50000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211328/cms_tot_8933_50000001_60000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211328/cms_tot_8933_60000001_70000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211328/cms_tot_8933_70000001_80000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211328/cms_tot_8933_80000001_90000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211328/cms_tot_8933_90000001_100000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211421/cms_tot_8946_30000001_40000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211421/cms_tot_8946_40000001_50000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211421/cms_tot_8946_50000001_60000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211421/cms_tot_8946_60000001_70000000.root");
	  chainTree->Add("/data/HI13/CmsTotemData/211421/cms_tot_8946_70000001_80000000.root");
      };
      break;
    };

    case 1:{ // Hijing ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zdc56 = false;
      if(ppb){
	  str="Hijing_pPb";
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_Hijing_PPb502_MB_pa_STARTHI53_V25v1__04043DE5-DA75-E211-95B9-00266CFAE1EC.root");  
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_Hijing_PPb502_MB_pa_STARTHI53_V25v1__0674A051-D075-E211-A67E-848F69FD290A.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_Hijing_PPb502_MB_pa_STARTHI53_V25v1__06AB7BC9-E175-E211-8001-848F69FD294F.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_Hijing_PPb502_MB_pa_STARTHI53_V25v1__06CE93DE-CD75-E211-9F44-008CFA0016A4.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_Hijing_PPb502_MB_pa_STARTHI53_V25v1__0881AC11-1D76-E211-A11D-00266CF256CC.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_Hijing_PPb502_MB_pa_STARTHI53_V25v1__0A45972D-D375-E211-BF77-008CFA0004B0.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_Hijing_PPb502_MB_pa_STARTHI53_V25v1__0C6030FD-D975-E211-AF0B-0024E8768D68.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_Hijing_PPb502_MB_pa_STARTHI53_V25v1__0E085C45-D075-E211-A918-00266CF9B630.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_Hijing_PPb502_MB_pa_STARTHI53_V25v1__105358BC-CE75-E211-B32E-00266CF279F8.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_Hijing_PPb502_MB_pa_STARTHI53_V25v1__1055B5DD-D075-E211-8BFF-008CFA001D7C.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_Hijing_PPb502_MB_pa_STARTHI53_V25v1__10D2E05E-CF75-E211-A34D-008CFA001444.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_Hijing_PPb502_MB_pa_STARTHI53_V25v1__10D9A79C-D275-E211-9022-0024E87699F4.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_Hijing_PPb502_MB_pa_STARTHI53_V25v1__10DCAF62-CF75-E211-8575-0024E8768CBF.root");
      }else{
	  str="Hijing_Pbp";
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_HijingPbp100k_dNdetaCentralityStudy_v1_job10-19.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_HijingPbp100k_dNdetaCentralityStudy_v1_job20-29.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_HijingPbp100k_dNdetaCentralityStudy_v1_job30-39.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_HijingPbp100k_dNdetaCentralityStudy_v1_job40-49.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_HijingPbp100k_dNdetaCentralityStudy_v1_job50-59.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_HijingPbp100k_dNdetaCentralityStudy_v1_job60-69.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_HijingPbp100k_dNdetaCentralityStudy_v1_job70-79.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_HijingPbp100k_dNdetaCentralityStudy_v1_job80-89.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_HijingPbp100k_dNdetaCentralityStudy_v1_job90-99.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_HijingPbp100k_dNdetaCentralityStudy_v1_job00-09.root");
      };
      break;
    };

    case 2:{ // EPOS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zdc56 = false;
      if(ppb){
	  str="EPOS_pPb";
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__004C2C1B-0777-E211-B514-0024E86E8D0B.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__005F009C-F476-E211-9E70-848F69FD2484.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__021E903B-4176-E211-990A-00266CFAE050.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__045598A5-EF76-E211-ADA7-180373FF8446.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__0676A8A4-EF76-E211-B587-00266CF9AF00.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__069A9303-E176-E211-B8D4-00266CF9B874.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__0867EC3F-DC76-E211-9AE6-00A0D1EE8C64.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__0887049B-1177-E211-9BBA-008CFA0004B0.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__08A7AFCA-2A77-E211-B4C4-008CFA001DB8.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__08A95CA0-F676-E211-A40A-008CFA008D0C.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__08C34118-3A76-E211-8FE0-00266CF25708.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__0A16086C-4376-E211-97B8-00A0D1EEF364.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__0A169B21-E176-E211-8182-00A0D1EE9238.root");
	  chainTree->Add("/data/HI13/MC/pPb/UABaseTree__HI_MC_EPOS_PPb502_MB_pa_STARTHI53_V25v1__0AF1E8E7-C476-E211-BE7D-848F69FD29CA.root");
      }else{
	  str="EPOS_Pbp";
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_EposPbp100k_dNdetaCentralityStudy_v3_job00-09.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_EposPbp100k_dNdetaCentralityStudy_v3_job10-19.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_EposPbp100k_dNdetaCentralityStudy_v3_job20-29.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_EposPbp100k_dNdetaCentralityStudy_v3_job30-39.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_EposPbp100k_dNdetaCentralityStudy_v3_job40-49.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_EposPbp100k_dNdetaCentralityStudy_v3_job50-59.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_EposPbp100k_dNdetaCentralityStudy_v3_job60-69.root");
	  chainTree->Add("/data/HI13/MC/Pbp/UABaseTree__tuos_EposPbp100k_dNdetaCentralityStudy_v3_job70-79.root");
      };
      break;
    };

    default:{
      zdc56 = false;
      std::cout << "uaplotter1::initializeChain(): unknown value for MC type = " << mc << std::endl;
      std::cout << "uaplotter1::initializeChain(): the TChain will not be initialized; exiting.\n";
      exit(10);
    };
  };
  return str;
};