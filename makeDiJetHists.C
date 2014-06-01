//=============================================                                                
// Author: Chris McGinn                                                                        
//                                                                                             
// DiJet Initial Skim Class (MC)                                                               
//                                                                                             
// !!NOTE: Written for jets sorted by pt, tracks unsorted!!                                    
//                                                                                             
//=============================================     

#include "commonUtility.h"
#include "TTree.h"
#include "TDatime.h"
#include "TFile.h"
#include "TCut.h"
#include "TProfile.h"
#include "TH1F.h"
#include "diJetFileTag.h"


Float_t getDPHI( Float_t phi1, Float_t phi2){
  Float_t dphi = phi1 - phi2;

  if(dphi > TMath::Pi())
    dphi = dphi - 2.*(TMath::Pi());
  if(dphi <= -(TMath::Pi()) )
    dphi = dphi + 2.*(TMath::Pi());

  if(TMath::Abs(dphi) > TMath::Pi())
    std::cout << " commonUtility::getDPHI error!!! dphi is bigger than TMath::Pi() " << std::endl;

  return dphi;
}


Float_t getAbsDphi(Float_t phi1, Float_t phi2){
  return TMath::Abs(getDPHI(phi1, phi2));
}


Bool_t sameSign(Double_t num1, Double_t num2){
  if((num1 > 0 && num2 > 0) || (num1 < 0 && num2 < 0)) return true;

  return false;
}


void niceTH1(TH1F* uglyTH1, float max , float min, float ndivX, float ndivY)
{
  handsomeTH1N(uglyTH1);
  uglyTH1->SetMaximum(max);
  uglyTH1->SetMinimum(min);
  uglyTH1->SetNdivisions(ndivX);
  uglyTH1->SetNdivisions(ndivY, "Y");
}


Bool_t checkSetRange(Int_t setNum)
{
  if(setNum > 2 || setNum < 0){
    std::cout << "checkSetRange: setNum must be between 0-2, empty cut returned" << std::endl;
    return false;
  }

  return true;
}


TCut makeSetCut(Int_t setNum)
{
  if(!checkSetRange(setNum))
    return "";

  return Form("eventSet[%d]", setNum);
}


TCut makeCentCut(Int_t centLow, Int_t centHi)
{
  if(centLow >= 0 && centHi >= centLow && centHi <= 199)
    return Form("hiBin >= %d && hiBin <= %d", centLow, centHi);
  else{
    std::cout << "makeCentCut: centLow/centHi incorrectly specified, empty cut returned" << std::endl;
    return "";
  }
}


TCut makeAsymmCut(Int_t setNum, Float_t asymmLow, Float_t asymmHi)
{
  if(!checkSetRange(setNum))
    return "";

  const char* cutVar = Form("AlgJtAsymm[%d]", setNum);

  if(asymmLow >= .00 && asymmHi >= asymmLow && asymmHi <= 1.)
    return Form("%s > %f && %s < %f ", cutVar, asymmLow, cutVar, asymmHi);
  else{
    std::cout << "makeAsymmCut: asymmLow/asymmHi incorrectly specified, empty cut returned" << std::endl;
    return "";
  }
}


TCut makeEtaCut(Int_t setNum, Float_t overallCut = 2.0)
{
  if(!checkSetRange(setNum))
    return "";

  const char* leadJt = Form("AlgLeadJtEta[%d]", setNum);
  const char* subLeadJt = Form("AlgSubLeadJtEta[%d]", setNum);

  return Form("TMath::Abs(%s) < %f && TMath::Abs(%s) < %f", leadJt, overallCut, subLeadJt, overallCut);
}


TCut makeDelPhiCut(Int_t setNum, Float_t delPhiLow = 0)
{
  if(!checkSetRange(setNum))
    return "";

  const char* jtDelPhi = Form("AlgJtDelPhi[%d]", setNum);

  return Form("%s > %f", jtDelPhi, delPhiLow);
}


void makeImbAsymmHist(TTree* anaTree_p, const char* outName, const char* gorr, Int_t setNum, const char* CNC, const char* FPT, Int_t centLow, Int_t centHi, Int_t histLow, Int_t histHi, const char* Corr = "", Bool_t montecarlo = false)
{
  Int_t setCorrNum = setNum;
  if(!strcmp("", Corr))
    setCorrNum = setNum + 3;

  const char* title = Form();
}
