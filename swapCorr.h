#include "TFile.h"
#include "TF1.h"
#include "/net/hisrv0001/home/cfmcginn/DijetMPt/CMSSW_5_3_20/src/DijetInitialSkim/sType.h"

const Int_t nFileSWAPPbPb = 4;
const Int_t nFileSWAPPP = 4;

TFile* SWAP_VsCaloFile_p[nFileSWAPPbPb];
TF1* SWAP12_VsCaloCorr_010_f[nFileSWAPPbPb];
TF1* SWAP12_VsCaloCorr_1030_f[nFileSWAPPbPb];
TF1* SWAP12_VsCaloCorr_3050_f[nFileSWAPPbPb];
TF1* SWAP12_VsCaloCorr_5070_f[nFileSWAPPbPb];
TF1* SWAP12_VsCaloCorr_70100_f[nFileSWAPPbPb];
TF1* SWAP23_VsCaloCorr_010_f[nFileSWAPPbPb];
TF1* SWAP23_VsCaloCorr_1030_f[nFileSWAPPbPb];
TF1* SWAP23_VsCaloCorr_3050_f[nFileSWAPPbPb];
TF1* SWAP23_VsCaloCorr_5070_f[nFileSWAPPbPb];
TF1* SWAP23_VsCaloCorr_70100_f[nFileSWAPPbPb];

TFile* SWAP_CaloFile_p[nFileSWAPPP];
TF1* SWAP12_CaloCorr_PP_f[nFileSWAPPP];
TF1* SWAP23_CaloCorr_PP_f[nFileSWAPPP];

void InitSWAPCorrFiles(sampleType sType = kHIDATA)
{
  if(sType == kHIDATA || sType == kHIMC){
    // SWAP_VsCaloFile_p[0] = new TFile("swapPlot_v6_NPFEffCor_R2.root");
    // SWAP_VsCaloFile_p[1] = new TFile("swapPlot_v6_NPFEffCor_R3.root");
    // SWAP_VsCaloFile_p[2] = new TFile("swapPlot_v6_NPFEffCor_R4.root");
    // SWAP_VsCaloFile_p[3] = new TFile("swapPlot_v6_NPFEffCor_R5.root");
    SWAP_VsCaloFile_p[0] = new TFile("swapPlot_v5_R2.root");
    SWAP_VsCaloFile_p[1] = new TFile("swapPlot_v5_R3.root");
    SWAP_VsCaloFile_p[2] = new TFile("swapPlot_v5_R4.root");
    SWAP_VsCaloFile_p[3] = new TFile("swapPlot_v5_R5.root");
  } 
  else if(sType == kPPDATA || sType == kPPMC){
    SWAP_CaloFile_p[0] = new TFile("swapPlot_PP_v5_R2.root");
    SWAP_CaloFile_p[1] = new TFile("swapPlot_PP_v5_R3.root");
    SWAP_CaloFile_p[2] = new TFile("swapPlot_PP_v5_R4.root");
    SWAP_CaloFile_p[3] = new TFile("swapPlot_PP_v5_R5.root");
  }  

  return;
}

void InitSWAPCorrFits(sampleType sType = kHIDATA)
{
  if(sType == kHIDATA || sType == kHIMC){
    for(Int_t iter = 0; iter < 4; iter++){
      SWAP12_VsCaloCorr_010_f[iter] = (TF1*)SWAP_VsCaloFile_p[iter]->Get("fPLeading_0_20");
      SWAP12_VsCaloCorr_1030_f[iter] = (TF1*)SWAP_VsCaloFile_p[iter]->Get("fPLeading_20_60");
      SWAP12_VsCaloCorr_3050_f[iter] = (TF1*)SWAP_VsCaloFile_p[iter]->Get("fPLeading_60_100");
      SWAP12_VsCaloCorr_5070_f[iter] = (TF1*)SWAP_VsCaloFile_p[iter]->Get("fPLeading_100_140");
      SWAP12_VsCaloCorr_70100_f[iter] = (TF1*)SWAP_VsCaloFile_p[iter]->Get("fPLeading_140_200");

      SWAP23_VsCaloCorr_010_f[iter] = (TF1*)SWAP_VsCaloFile_p[iter]->Get("fPSubleading_0_20");
      SWAP23_VsCaloCorr_1030_f[iter] = (TF1*)SWAP_VsCaloFile_p[iter]->Get("fPSubleading_20_60");
      SWAP23_VsCaloCorr_3050_f[iter] = (TF1*)SWAP_VsCaloFile_p[iter]->Get("fPSubleading_60_100");
      SWAP23_VsCaloCorr_5070_f[iter] = (TF1*)SWAP_VsCaloFile_p[iter]->Get("fPSubleading_100_140");
      SWAP23_VsCaloCorr_70100_f[iter] = (TF1*)SWAP_VsCaloFile_p[iter]->Get("fPSubleading_140_200");
    }
  }
  else if(sType == kPPDATA || sType == kPPMC){
    for(Int_t iter = 0; iter < 4; iter++){
      SWAP12_CaloCorr_PP_f[iter] = (TF1*)SWAP_CaloFile_p[iter]->Get("fPLeading_0_20");
      SWAP23_CaloCorr_PP_f[iter] = (TF1*)SWAP_CaloFile_p[iter]->Get("fPSubleading_0_20");
    }
  }

  return;
}


Float_t GetJtSWAP12Prob(sampleType sType, Int_t jtRBin, Int_t hiBin, Float_t Asymm12)
{
  Float_t prob = 1.0;
  const Int_t hiBinCut[5] = {19, 59, 99, 139, 199};


  if(sType == kHIDATA || sType == kHIMC){

    if(hiBin <= hiBinCut[0]) prob = SWAP12_VsCaloCorr_010_f[jtRBin]->Eval(Asymm12);
    else if(hiBin <= hiBinCut[1]) prob = SWAP12_VsCaloCorr_1030_f[jtRBin]->Eval(Asymm12);
    else if(hiBin <= hiBinCut[2]) prob = SWAP12_VsCaloCorr_3050_f[jtRBin]->Eval(Asymm12);
    else if(hiBin <= hiBinCut[3]) prob = SWAP12_VsCaloCorr_5070_f[jtRBin]->Eval(Asymm12);
    else prob = SWAP12_VsCaloCorr_70100_f[jtRBin]->Eval(Asymm12);
  }
  else if(sType == kPPDATA || sType == kPPMC) 
  prob = SWAP12_CaloCorr_PP_f[jtRBin]->Eval(Asymm12);
  
  return prob;
}

Float_t GetJtSWAP23Prob(sampleType sType, Int_t jtRBin, Int_t hiBin, Float_t Asymm23)
{
  Float_t prob = 1.0;
  const Int_t hiBinCut[5] = {19, 59, 99, 139, 199};


  if(sType == kHIDATA || sType == kHIMC){

    if(hiBin <= hiBinCut[0]) prob = SWAP23_VsCaloCorr_010_f[jtRBin]->Eval(Asymm23);
    else if(hiBin <= hiBinCut[1]) prob = SWAP23_VsCaloCorr_1030_f[jtRBin]->Eval(Asymm23);
    else if(hiBin <= hiBinCut[2]) prob = SWAP23_VsCaloCorr_3050_f[jtRBin]->Eval(Asymm23);
    else if(hiBin <= hiBinCut[3]) prob = SWAP23_VsCaloCorr_5070_f[jtRBin]->Eval(Asymm23);
    else prob = SWAP23_VsCaloCorr_70100_f[jtRBin]->Eval(Asymm23);
  }
  else if(sType == kPPDATA || sType == kPPMC) 
  prob = SWAP23_CaloCorr_PP_f[jtRBin]->Eval(Asymm23);
  
  return prob;
}
