//=============================================
// Author: Chris McGinn
//
// DiJet Analysis Class (MC)
//
//=============================================     


#include "TFile.h"
#include "TH2D.h"
#include "etaPhiFunc.h"
#include "sType.h"

const Int_t nFilePbPb = 4;
const Int_t nFilePP = 4;

TFile* JEC_VsCaloFile_p[nFilePbPb];
TH2D* JEC_VsCaloCorr_010_h[nFilePbPb];
TH2D* JEC_VsCaloCorr_1030_h[nFilePbPb];
TH2D* JEC_VsCaloCorr_3050_h[nFilePbPb];
TH2D* JEC_VsCaloCorr_50100_h[nFilePbPb];


TFile* JEC_CaloFile_p[nFilePP];
TH2D* JEC_CaloCorr_PP_h[nFilePP];

void InitJECCorrFiles(sampleType sType = kHIDATA)
{
  if(sType == kHIDATA || sType == kHIMC){
    JEC_VsCaloFile_p[0] = new TFile("FFJEC_correction_PF_akVs2Calo_pt2.root");
    JEC_VsCaloFile_p[1] = new TFile("FFJEC_correction_PF_akVs3Calo_pt2.root");
    JEC_VsCaloFile_p[2] = new TFile("FFJEC_correction_PF_akVs4Calo_pt2.root");
    JEC_VsCaloFile_p[3] = new TFile("FFJEC_correction_PF_akVs5Calo_pt2.root");
  }
  else if(sType == kPPDATA || sType == kPPMC){
    JEC_CaloFile_p[0] = new TFile("FFJEC_correction_PF_ak2Calo_pt2.root");
    JEC_CaloFile_p[1] = new TFile("FFJEC_correction_PF_ak3Calo_pt2.root");
    JEC_CaloFile_p[2] = new TFile("FFJEC_correction_PF_ak4Calo_pt2.root");
    JEC_CaloFile_p[3] = new TFile("FFJEC_correction_PF_ak5Calo_pt2.root");
  }
  
  return;
}

void InitJECCorrHists(sampleType sType = kHIDATA)
{
  if(sType == kHIDATA || sType == kHIMC){
    for(Int_t iter = 0; iter < 4; iter++){
      JEC_VsCaloCorr_010_h[iter] = (TH2D*)JEC_VsCaloFile_p[iter]->Get("pNtrk_pt0");
      JEC_VsCaloCorr_1030_h[iter] = (TH2D*)JEC_VsCaloFile_p[iter]->Get("pNtrk_pt1");
      JEC_VsCaloCorr_3050_h[iter] = (TH2D*)JEC_VsCaloFile_p[iter]->Get("pNtrk_pt2");
      JEC_VsCaloCorr_50100_h[iter] = (TH2D*)JEC_VsCaloFile_p[iter]->Get("pNtrk_pt3");
    }
  }
  else if(sType == kPPDATA || sType == kPPMC){
    for(Int_t iter = 0; iter < 4; iter++){
      JEC_CaloCorr_PP_h[iter] = (TH2D*)JEC_CaloFile_p[iter]->Get("pNtrk_pt");
    }
  }

  return;
}


Float_t Get2PFCand(Float_t jtR, Float_t jtPhi, Float_t jtEta, Int_t nPFpart, Float_t pfPt[], Int_t pfID[], Float_t pfPhi[], Float_t pfEta[])
{
  Int_t nInJt = 0;
  const Float_t pfPtCut = 2.0;

  for(Int_t iter = 0; iter < nPFpart; iter++){
    if(pfPt[iter] < pfPtCut) continue;
    if(pfID[iter] != 1) continue;

    if(getDR(jtEta, jtPhi, pfEta[iter], pfPhi[iter]) < jtR) nInJt++;
  }

  return nInJt;
}


Float_t GetJtCorrPt(sampleType sType, Int_t jtRBin, Int_t hiBin, Float_t jtPt, Float_t nPF)
{
  Float_t corrPt = 1.0;
  const Int_t hiBinCut[4] = {19, 59, 99, 199};

  if(jtPt > 20){
    if(sType == kHIDATA || sType == kHIMC){
      if(hiBin <= hiBinCut[0]) corrPt = JEC_VsCaloCorr_010_h[jtRBin]->GetBinContent(JEC_VsCaloCorr_010_h[jtRBin]->GetXaxis()->FindBin(nPF), JEC_VsCaloCorr_010_h[jtRBin]->GetYaxis()->FindBin(jtPt));
      else if(hiBin <= hiBinCut[1]) corrPt = JEC_VsCaloCorr_1030_h[jtRBin]->GetBinContent(JEC_VsCaloCorr_1030_h[jtRBin]->GetXaxis()->FindBin(nPF), JEC_VsCaloCorr_1030_h[jtRBin]->GetYaxis()->FindBin(jtPt));
      else if(hiBin <= hiBinCut[2]) corrPt = JEC_VsCaloCorr_3050_h[jtRBin]->GetBinContent(JEC_VsCaloCorr_3050_h[jtRBin]->GetXaxis()->FindBin(nPF), JEC_VsCaloCorr_3050_h[jtRBin]->GetYaxis()->FindBin(jtPt));
      else corrPt = JEC_VsCaloCorr_50100_h[jtRBin]->GetBinContent(JEC_VsCaloCorr_50100_h[jtRBin]->GetXaxis()->FindBin(nPF), JEC_VsCaloCorr_50100_h[jtRBin]->GetYaxis()->FindBin(jtPt));
    }
    else if(sType == kPPDATA || sType == kPPMC){
      corrPt = JEC_CaloCorr_PP_h[jtRBin]->GetBinContent(JEC_CaloCorr_PP_h[jtRBin]->GetXaxis()->FindBin(nPF), JEC_CaloCorr_PP_h[jtRBin]->GetYaxis()->FindBin(jtPt));
    }
  }

  return 1/corrPt;
}
