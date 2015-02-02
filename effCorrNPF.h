//=============================================
// Author: Chris McGinn
//
// JEC EFF. Corr Class (MC)
//
//=============================================     


#include "TFile.h"
#include "TH2D.h"
#include "etaPhiFunc.h"
#include "sType.h"

const Int_t nFileEFFPbPb = 4;

TFile* EFF_VsCaloFile_p[nFileEFFPbPb];
TH1D* EFF_Eta_VsCaloCorr_010_h[nFileEFFPbPb];
TH1D* EFF_Eta_VsCaloCorr_1030_h[nFileEFFPbPb];
TH1D* EFF_Eta_VsCaloCorr_3050_h[nFileEFFPbPb];
TH1D* EFF_Eta_VsCaloCorr_50100_h[nFileEFFPbPb];
TH1D* EFF_Pt_VsCaloCorr_010_h[nFileEFFPbPb];
TH1D* EFF_Pt_VsCaloCorr_1030_h[nFileEFFPbPb];
TH1D* EFF_Pt_VsCaloCorr_3050_h[nFileEFFPbPb];
TH1D* EFF_Pt_VsCaloCorr_50100_h[nFileEFFPbPb];

void InitEFFCorrFiles()
{
  EFF_VsCaloFile_p[0] = new TFile("eff_NPF_PH_jtpt_eta_akVs2Calo.root");
  EFF_VsCaloFile_p[1] = new TFile("eff_NPF_PH_jtpt_eta_akVs3Calo.root");
  EFF_VsCaloFile_p[2] = new TFile("eff_NPF_PH_jtpt_eta_akVs4Calo.root");
  EFF_VsCaloFile_p[3] = new TFile("eff_NPF_PH_jtpt_eta_akVs5Calo.root");
  
  return;
}

void InitEFFCorrHists()
{
  for(Int_t iter = 0; iter < 4; iter++){
    EFF_Eta_VsCaloCorr_010_h[iter] = (TH1D*)EFF_VsCaloFile_p[iter]->Get("heta_cent0");
    EFF_Eta_VsCaloCorr_1030_h[iter] = (TH1D*)EFF_VsCaloFile_p[iter]->Get("heta_cent1");
    EFF_Eta_VsCaloCorr_3050_h[iter] = (TH1D*)EFF_VsCaloFile_p[iter]->Get("heta_cent2");
    EFF_Eta_VsCaloCorr_50100_h[iter] = (TH1D*)EFF_VsCaloFile_p[iter]->Get("heta_cent3");
    EFF_Pt_VsCaloCorr_010_h[iter] = (TH1D*)EFF_VsCaloFile_p[iter]->Get("hpt_cent0");
    EFF_Pt_VsCaloCorr_1030_h[iter] = (TH1D*)EFF_VsCaloFile_p[iter]->Get("hpt_cent1");
    EFF_Pt_VsCaloCorr_3050_h[iter] = (TH1D*)EFF_VsCaloFile_p[iter]->Get("hpt_cent2");
    EFF_Pt_VsCaloCorr_50100_h[iter] = (TH1D*)EFF_VsCaloFile_p[iter]->Get("hpt_cent3");
  }
  return;
}


Float_t GetEFFCorr(Int_t jtRBin, Float_t jtPt, Float_t jtEta, Int_t hiBin)
{
  Float_t effCorr = 1.0;
  const Int_t hiBinCut[4] = {19, 59, 99, 199};

  Float_t jtPtCheck = jtPt;
  if(jtPtCheck < 15) jtPtCheck = 15;

  if(hiBin <= hiBinCut[0]){
    effCorr = EFF_Eta_VsCaloCorr_010_h[jtRBin]->GetBinContent(EFF_Eta_VsCaloCorr_010_h[jtRBin]->FindBin(jtEta));
    effCorr *= EFF_Pt_VsCaloCorr_010_h[jtRBin]->GetBinContent(EFF_Pt_VsCaloCorr_010_h[jtRBin]->FindBin(jtPtCheck));
  }else if(hiBin <= hiBinCut[1]){
    effCorr = EFF_Eta_VsCaloCorr_1030_h[jtRBin]->GetBinContent(EFF_Eta_VsCaloCorr_1030_h[jtRBin]->FindBin(jtEta));
    effCorr *= EFF_Pt_VsCaloCorr_1030_h[jtRBin]->GetBinContent(EFF_Pt_VsCaloCorr_1030_h[jtRBin]->FindBin(jtPtCheck));
  }else if(hiBin <= hiBinCut[2]){ 
    effCorr = EFF_Eta_VsCaloCorr_3050_h[jtRBin]->GetBinContent(EFF_Eta_VsCaloCorr_3050_h[jtRBin]->FindBin(jtEta));
    effCorr *= EFF_Pt_VsCaloCorr_3050_h[jtRBin]->GetBinContent(EFF_Pt_VsCaloCorr_3050_h[jtRBin]->FindBin(jtPtCheck));
  }else{
    effCorr = EFF_Eta_VsCaloCorr_50100_h[jtRBin]->GetBinContent(EFF_Eta_VsCaloCorr_50100_h[jtRBin]->FindBin(jtEta));
    effCorr *= EFF_Pt_VsCaloCorr_50100_h[jtRBin]->GetBinContent(EFF_Pt_VsCaloCorr_50100_h[jtRBin]->FindBin(jtPtCheck));
  }
  
  return 1/effCorr;
}