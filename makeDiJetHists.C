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

const char* algType[5] = {"PuCalo", "VsCalo", "T"};

TFile* inFile_p = 0;
TFile* outFile_p = 0;

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


void niceTH1(TH1F* uglyTH1, float max , float min)
{
  handsomeTH1N(uglyTH1);
  uglyTH1->SetMaximum(max);
  uglyTH1->SetMinimum(min);
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
  inFile_p->cd();

  Int_t setCorrNum = setNum;
  if(!strcmp("", Corr))
    setCorrNum = setNum + 3;

  const char* title = Form("%s%sImbAsymmProjA%s%s%s_%d%d_%s_h", gorr, algType[setNum], CNC, FPT, Corr, (Int_t)(centLow*.5), (Int_t)((centHi + 1)*.5), fileTag);

  Float_t xArr[5] = {.0, .11, .22, .33, .50};

  TH1F* imbAsymmHist_p = new TH1F("imbAsymmHist_p", "imbAsymmHist_p", 4, xArr);
  imbAsymmHist_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTH1(imbAsymmHist_p, histHi, histLow);

  TH1F* getHist_p;

  TString var = Form("%sAlgImbProjA%s%s[%d]", gorr, CNC, FPT, setCorrNum);

  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut etaCut = makeEtaCut(setNum, 1.6);
  TCut phiCut = makeDelPhiCut(setNum, 5*TMath::Pi()/6);

  TCut jetLCut = Form("AlgLeadJtPt[%d] > 120", setNum);
  TCut jetSLCut = Form("AlgSubLeadJtPt[%d] > 120", setNum);
  TCut pthat = "pthat > 80";

  const char* name1[4] = {"0_1(10000, -10000, 10000)", "1_2(10000, -10000, 10000)", "2_3(10000, -10000, 10000)", "3_5(10000, -10000, 10000)"};
  const char* name2[4] = {"0_1", "1_2", "2_3", "3_5"};
  Float_t asymmBins[5] = {.00, .11, .22, .33, 1.00};

  for(Int_t binIter = 0; binIter < 4; binIter++){
    TCut asymmCut = makeAsymmCut(setNum, asymmBins[binIter], asymmBins[binIter + 1]);

    if(montecarlo)
      anaTree_p->Project(name1[binIter], var, setCut && centCut && etaCut && phiCut && jetLCut && jetSLCut && asymmCut, "");
    else
      anaTree_p->Project(name1[binIter], var, setCut && centCut && etaCut && phiCut && jetLCut && jetSLCut && asymmCut, "");

    getHist_p = (TH1F*)inFile_p->Get(name2[binIter]);

    imbAsymmHist_p->SetBinContent(binIter + 1, getHist_p->GetMean());
    imbAsymmHist_p->SetBinError(binIter + 1, getHist_p->GetMeanError());
  }

  outFile_p = new TFile(outName, "UPDATE");
  imbAsymmHist_p->Write(title);
  outFile_p->Close();

  delete outFile_p;
  delete imbAsymmHist_p;

  return;
}


void makeDiJetHists(const char* inName, const char* outName, Bool_t montecarlo = false)
{
  TH1::SetDefaultSumw2();

  setFileTag(inName);

  inFile_p = new TFile(inName, "READ");
  TTree* anaTree_p = (TTree*)inFile_p->Get("jetTreeAna");
  anaTree_p->AddFriend("trackTreeAna");

  Int_t jetAlgMax = 2;

  if(montecarlo){
    anaTree_p->AddFriend("genTreeAna");
    jetAlgMax = 3;
  }

  const char* corr[2] = {"", "Corr"};
  const char* CNC[3] = {"", "C", "NC"};
  const char* FPT[6] = {"F", "0_1", "1_2", "2_4", "4_8", "8_100"};

  Int_t centLow[6] = {0, 20, 60, 100, 0, 60};
  Int_t centHi[6] = {19, 59, 99, 199, 59, 199};

  for(Int_t algIter = 0; algIter < jetAlgMax; algIter++){
    std::cout << "Algorithm: " << algType[algIter] << std::endl;

    for(Int_t corrIter = 0; corrIter < 2; corrIter++){
      for(Int_t CNCIter = 0; CNCIter < 3; CNCIter++){
	for(Int_t centIter = 0; centIter < 6; centIter++){
	  for(Int_t FPTIter = 0; FPTIter < 6; FPTIter++){

	    makeImbAsymmHist(anaTree_p, outName, "r", algIter, CNC[CNCIter], FPT[FPTIter], centLow[centIter], centHi[centIter], -60, 60, corr[corrIter], montecarlo);

	    if(montecarlo && corrIter == 0)
	      makeImbAsymmHist(anaTree_p, outName, "g", algIter, CNC[CNCIter], FPT[FPTIter], centLow[centIter], centHi[centIter], -60, 60, corr[corrIter], montecarlo);

	  }
	}
      }
    }

  }

  inFile_p->Close();
  delete inFile_p;

  return;
}
