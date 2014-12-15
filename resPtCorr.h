//=============================================
// Author: Chris McGinn
//
// JEC Res. Corr. Class (MC)
//
//=============================================     

#include "TFile.h"
#include "TF1.h"
#include "sType.h"

const Int_t nFileRESPbPb = 4;
const Int_t nFileRESPP = 4;

TFile* RES_VsCaloFile_p[nFileRESPbPb];
TF1* RES_VsCaloCorr_010_f[nFileRESPbPb];
TF1* RES_VsCaloCorr_1030_f[nFileRESPbPb];
TF1* RES_VsCaloCorr_3050_f[nFileRESPbPb];
TF1* RES_VsCaloCorr_50100_f[nFileRESPbPb];

TFile* RES_CaloFile_p[nFileRESPP];
TF1* RES_CaloCorr_PP_f[nFileRESPP];

void InitRESCorrFiles(sampleType sType = kHIDATA)
{
  if(sType == kHIDATA || sType == kHIMC){
    RES_VsCaloFile_p[0] = new TFile("residualcorr_akVs2Calo.root");
    RES_VsCaloFile_p[1] = new TFile("residualcorr_akVs3Calo.root");
    RES_VsCaloFile_p[2] = new TFile("residualcorr_akVs4Calo.root");
    RES_VsCaloFile_p[3] = new TFile("residualcorr_akVs5Calo.root");
  }
  else if(sType == kPPDATA || sType == kPPMC){
    RES_CaloFile_p[0] = new TFile("residualcorr_ak2Calo.root");
    RES_CaloFile_p[1] = new TFile("residualcorr_ak3Calo.root");
    RES_CaloFile_p[2] = new TFile("residualcorr_ak4Calo.root");
    RES_CaloFile_p[3] = new TFile("residualcorr_ak5Calo.root");
  }  

  return;
}

void InitRESCorrFits(sampleType sType = kHIDATA)
{
  if(sType == kHIDATA || sType == kHIMC){
    for(Int_t iter = 0; iter < 4; iter++){
      RES_VsCaloCorr_010_f[iter] = (TF1*)RES_VsCaloFile_p[iter]->Get("fit0");
      RES_VsCaloCorr_1030_f[iter] = (TF1*)RES_VsCaloFile_p[iter]->Get("fit1");
      RES_VsCaloCorr_3050_f[iter] = (TF1*)RES_VsCaloFile_p[iter]->Get("fit2");
      RES_VsCaloCorr_50100_f[iter] = (TF1*)RES_VsCaloFile_p[iter]->Get("fit3");
    }
  }
  else if(sType == kPPDATA || sType == kPPMC){
    for(Int_t iter = 0; iter < 4; iter++){
      RES_CaloCorr_PP_f[iter] = (TF1*)RES_CaloFile_p[iter]->Get("fit0");
    }
  }

  return;
}


Float_t GetJtRESCorrPt(sampleType sType, Int_t jtRBin, Int_t hiBin, Float_t jtPt)
{
  Float_t corrPt = 1.0;
  const Int_t hiBinCut[4] = {19, 59, 99, 199};

  if(jtPt > 20 && jtPt < 700){
    if(sType == kHIDATA || sType == kHIMC){
      if(hiBin <= hiBinCut[0]) corrPt = RES_VsCaloCorr_010_f[jtRBin]->Eval(jtPt);
      else if(hiBin <= hiBinCut[1]) corrPt = RES_VsCaloCorr_1030_f[jtRBin]->Eval(jtPt);
      else if(hiBin <= hiBinCut[2]) corrPt = RES_VsCaloCorr_3050_f[jtRBin]->Eval(jtPt);
      else corrPt = RES_VsCaloCorr_50100_f[jtRBin]->Eval(jtPt);
    }
    else if(sType == kPPDATA || sType == kPPMC){
      corrPt = RES_CaloCorr_PP_f[jtRBin]->Eval(jtPt);
    }
  }

  return 1/corrPt;
}
