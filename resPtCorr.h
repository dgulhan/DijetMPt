//=============================================
// Author: Chris McGinn
//
// JEC Res. Corr. Class (MC)
//
//=============================================     

#include "TFile.h"
#include "TF1.h"
#include "sType.h"

const Int_t nFilePbPb = 4;
const Int_t nFilePP = 4;

TFile* RES_VsCaloFile_p[nFilePbPb];
TF1* RES_VsCaloCorr_010_f[nFilePbPb];
TF1* RES_VsCaloCorr_1030_f[nFilePbPb];
TF1* RES_VsCaloCorr_3050_f[nFilePbPb];
TF1* RES_VsCaloCorr_50100_f[nFilePbPb];


void InitRESCorrFiles(sampleType sType = kHIDATA)
{
  if(sType == kHIDATA || sType == kHIMC){
    RES_VsCaloFile_p[0] = new TFile("residualcorr_akVs2Calo.root");
    RES_VsCaloFile_p[1] = new TFile("residualcorr_akVs3Calo.root");
    RES_VsCaloFile_p[2] = new TFile("residualcorr_akVs4Calo.root");
    RES_VsCaloFile_p[3] = new TFile("residualcorr_akVs5Calo.root");
  }
  
  return;
}

void InitRESCorrHists(sampleType sType = kHIDATA)
{
  if(sType == kHIDATA || sType == kHIMC){
    for(Int_t iter = 0; iter < 4; iter++){
      RES_VsCaloCorr_010_f[iter] = (TF1*)RES_VsCaloFile_p[iter]->Get("fit0");
      RES_VsCaloCorr_1030_f[iter] = (TF1*)RES_VsCaloFile_p[iter]->Get("fit1");
      RES_VsCaloCorr_3050_f[iter] = (TF1*)RES_VsCaloFile_p[iter]->Get("fit2");
      RES_VsCaloCorr_50100_f[iter] = (TF1*)RES_VsCaloFile_p[iter]->Get("fit3");
    }
  }

  return;
}


Float_t GetJtRESCorrPt(sampleType sType, Int_t jtRBin, Int_t hiBin, Float_t jtPt)
{
  Float_t corrPt = 1.0;
  const Int_t hiBinCut[4] = {19, 59, 99, 199};

  if(jtPt > 20){
    if(sType == kHIDATA || sType == kHIMC){
      if(hiBin <= hiBinCut[0]) corrPt = RES_VsCaloCorr_010_f[jtRBin]->Eval(jtPt);
      else if(hiBin <= hiBinCut[1]) corrPt = RES_VsCaloCorr_1030_f[jtRBin]->Eval(jtPt);
      else if(hiBin <= hiBinCut[2]) corrPt = RES_VsCaloCorr_3050_f[jtRBin]->Eval(jtPt);
      else corrPt = RES_VsCaloCorr_50100_f[jtRBin]->Eval(jtPt);
    }
  }

  return 1/corrPt;
}
