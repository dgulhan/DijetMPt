//=============================================                                                
// Author: Chris McGinn                                                                        
//                                                                                             
// DiJet Histogram Maker, Jet Properties                                                              
//                                                                                             
//=============================================     

#include "TTree.h"
#include "TDatime.h"
#include "TFile.h"
#include "diJetFileTag.h"

TFile* inFile_p = 0;
TFile* outFile_p = 0;

void makeAsymmHist(TTree* anaTree_p, const char* outName, Int_t setNum, Int_t nBins, Float_t histLow, Float_t histHi, Int_t centLow, Int_t centHi, sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  inFile_p->cd();

  Bool_t montecarlo = false;
  if(sType == kPPMC || sType == kHIMC)
    montecarlo = true;

  const char* title;

  if(sType == kPPMC || sType == kPPDATA)
    title = Form("%sAj_PP_%s", algType[setNum], fileTag);
  else
    title = Form("%sAj_%d%d_%s", algType[setNum], (Int_t)(centLow*.5), (Int_t)((centHi + 1)*.5), fileTag);

  TH1F* ajHist_p;

  TString name = Form("%s_h(%d, %f, %f)", title, nBins, histLow, histHi);
  InitCuts();
  SetCuts(setNum, sType, centLow, centHi, isHighPtTrk);

  if(montecarlo && setNum != 2)
    anaTree_p->Project(name, Form("AlgJtAsymm[%d]", setNum), Form("centWeight_80[%d]", setNum)*(centCut && setCut && etaCut && phiCut && jetLCut && jetSLCut && pthat && trkCut));
  else
    anaTree_p->Project(name, Form("AlgJtAsymm[%d]", setNum), centCut && setCut && etaCut && phiCut && jetLCut && jetSLCut && pthat && trkCut);

  ajHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));
  if(sType == kPPMC || sType == kPPDATA)
    niceTH1N(ajHist_p, .2999, 0., 405, 506, kBlue, 1, 25);
  else{
    if(sType == kHIDATA)
      niceTH1N(ajHist_p, .2999, 0., 405, 506, kRed, 1, 28);
    else{
      niceTH1N(ajHist_p, .2999, 0., 405, 506, 1, 0, 28);
      ajHist_p->SetFillColor(16);
    }
  }

  ajHist_p->SetYTitle("Event Fraction");
  ajHist_p->SetXTitle("A_{J}");


  outFile_p = new TFile(outName, "UPDATE");
  ajHist_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
  return;
}


void makeDiJetHists_Jet(const char* inName, const char* outName, sampleType sType = kHIDATA, Bool_t isPercent = false, Bool_t isHighPtTrk = false)
{
  TH1::SetDefaultSumw2();

  Bool_t montecarlo = false;
  if(sType == kHIMC || sType == kPPMC || sType == kPAMC)
    montecarlo = true;

  setFileTag(inName);

  inFile_p = new TFile(inName, "READ");
  TTree* anaTree_p = (TTree*)inFile_p->Get("jetTreeAna");

  Int_t jetAlgMax = 2;

  if(montecarlo)
    jetAlgMax = 3;

  Int_t centBins = 4;
  Int_t centLow[4] = {0, 20, 60, 100};
  Int_t centHi[4] = {19, 59, 99, 199};

  if(sType == kPPMC || sType == kPPDATA) centBins = 1;

  for(Int_t algIter = 1; algIter < jetAlgMax; algIter++){
    std::cout << "Algorithm: " << algType[algIter] << std::endl;

    for(Int_t centIter = 0; centIter < centBins; centIter++){
      makeAsymmHist(anaTree_p, outName, algIter, 20, 0.0001, 0.9999, centLow[centIter], centHi[centIter], sType, isHighPtTrk);
    }
  }

  inFile_p->Close();
  delete inFile_p;

  return;
}
