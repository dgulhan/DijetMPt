//=============================================                                                
// Author: Chris McGinn                                                                        
//                                                                                             
// DiJet Histogram Maker                                                              
//                                                                                             
//=============================================     

#include "TTree.h"
#include "TDatime.h"
#include "TFile.h"
#include "diJetFileTag.h"

TFile* inFile_p = 0;
TFile* outFile_p = 0;

Int_t leadJtCut = 120;
Int_t subLeadJtCut = 50;

const char* FPT[6] = {"0_1", "1_2", "2_4", "4_8", "8_100", "F"};

void makeImbAsymmHist(TTree* anaTree_p, const char* outName, const char* gorr, Int_t setNum, const char* CNCR, Int_t FPTNum, Int_t centLow, Int_t centHi, Float_t histLow, Float_t histHi, const char* Corr = "", const char* Tight = "", sampleType sType = kHIDATA, Bool_t isPercent = false, Bool_t isHighPtTrk = false)
{
  inFile_p->cd();

  Bool_t montecarlo = false;
  if(sType == kHIMC || sType == kPPMC || sType == kPAMC)
    montecarlo = true;

  Int_t setCorrNum = setNum;
  if(!strcmp("Corr", Corr))
    setCorrNum = setNum + 3;

  const char* title;

  if(sType == kHIDATA || sType == kHIMC)
    title = Form("%s%sImbProjA%s%s%s%s_%d%d_%s_h", gorr, algType[setNum], CNCR, FPT[FPTNum], Corr, Tight, (Int_t)(centLow*.5), (Int_t)((centHi + 1)*.5), fileTag);
  else
    title = Form("%s%sImbProjA%s%s%s%s_PP_%s_h", gorr, algType[setNum], CNCR, FPT[FPTNum], Corr, Tight, fileTag);

  Float_t xArr[5] = {.0001, .11, .22, .33, .4999};
  Float_t xArrPerc[5] = {.0001, 25., 50., 75., 99.9999};
  Float_t xArrTight[9] = {.0001, .055, .110, .165, .220, .275, .33, .415, .4999};
  Float_t xArrTightPerc[9] = {.0001, 12.5, 25., 37.5, 50., 62.5, 75., 87.5, 99.9999};

  TH1F* imbAsymmHist_p;
  if(!isPercent){
    if(!strcmp(Tight, ""))
      imbAsymmHist_p = new TH1F("imbAsymmHist_p", "imbAsymmHist_p", 4, xArr);
    else
      imbAsymmHist_p = new TH1F("imbAsymmHist_p", "imbAsymmHist_p", 8, xArrTight);
  }
  else{
    if(!strcmp(Tight, ""))
      imbAsymmHist_p = new TH1F("imbAsymmHist_p", "imbAsymmHist_p", 4, xArrPerc);
    else
      imbAsymmHist_p = new TH1F("imbAsymmHist_p", "imbAsymmHist_p", 8, xArrTightPerc);
  }

  if(!isPercent){
    imbAsymmHist_p->GetXaxis()->SetLimits(0.00, 0.50);
    niceTH1(imbAsymmHist_p, histHi, histLow, 505, 406);
  }
  else{
    imbAsymmHist_p->GetXaxis()->SetLimits(0.00, 100.00);
    niceTH1(imbAsymmHist_p, histHi, histLow, 504, 406);
  }

  TH1F* getHist_p;

  TString var = Form("%sAlgImbProjA%s[%d][%d]", gorr, CNCR, setCorrNum, FPTNum);

  TCut setCut = makeSetCut(setNum);
  TCut centCut = "";
  if(sType == kHIDATA || sType == kHIMC) centCut = makeCentCut(centLow, centHi);
  TCut etaCut = makeEtaCut(setNum, 1.6);
  TCut phiCut = makeDelPhiCut(setNum, 5*TMath::Pi()/6);

  TCut jetLCut = Form("AlgJtPt[%d][0] > %d", setNum, leadJtCut);
  TCut jetSLCut = Form("AlgJtPt[%d][1] > %d", setNum, subLeadJtCut);
  TCut pthat = "pthat > 80";

  TCut trkCut = "";
  if(isHighPtTrk)
    trkCut = Form("AlgJtTrkMax[%d][0] > 12 && AlgJtTrkMax[%d][1] > 12", setNum, setNum);

  const char* name1[8] = {"0(10000, -10000, 10000)", "1(10000, -10000, 10000)", "2(10000, -10000, 10000)", "3(10000, -10000, 10000)", "4(10000, -10000, 10000)", "5(10000, -10000, 10000)", "6(10000, -10000, 10000)", "7(10000, -10000, 10000)"};
  const char* name2[8] = {"0", "1", "2", "3", "4", "5", "6", "7"};

  Float_t asymmBins[5] = {.00, .11, .22, .33, 1.00};
  Float_t asymmBinsTight[9] = {.00, .055, .11, .165, .22, .275, .33, .415, 1.00};

  if(isPercent){
    if(!strcmp(Tight, ""))
      setAsymmPercBins(asymmBins, centLow, centHi, sType, Tight);
    else
      setAsymmPercBins(asymmBinsTight, centLow, centHi, sType, Tight);
  }

  Int_t nBins = 4;
  if(!strcmp(Tight, "Tight"))
    nBins = 8;

  for(Int_t binIter = 0; binIter < nBins; binIter++){
    TCut asymmCut = "";

    if(!strcmp(Tight, ""))
      asymmCut = makeAsymmCut(setNum, asymmBins[binIter], asymmBins[binIter + 1]);
    else
      asymmCut = makeAsymmCut(setNum, asymmBinsTight[binIter], asymmBinsTight[binIter + 1]);

    if(montecarlo)
      anaTree_p->Project(name1[binIter], var, setCut && centCut && etaCut && phiCut && jetLCut && jetSLCut && asymmCut && pthat && trkCut, "");
    else
      anaTree_p->Project(name1[binIter], var, setCut && centCut && etaCut && phiCut && jetLCut && jetSLCut && asymmCut && trkCut, "");

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


void makeImbDelRHist(TTree* anaTree_p, const char* outName, const char* gorr, Int_t setNum, const char* CNCR, Int_t FPTNum, Int_t centLow, Int_t centHi, Float_t histLow, Float_t histHi, const char* Corr = "", sampleType sType = kHIDATA, Bool_t isPercent = false, Bool_t isHighPtTrk = false)
{
  inFile_p->cd();

  Bool_t montecarlo = false;
  if(sType == kHIMC || sType == kPPMC || sType == kPAMC)
    montecarlo = true;

  Int_t setCorrNum = setNum;
  if(!strcmp("Corr", Corr))
    setCorrNum = setNum + 3;

  const char* title;

  if(sType == kHIDATA || sType == kHIMC)
    title = Form("%s%sImbProjA%s%s%s_%d%d_%s_h", gorr, algType[setNum], CNCR, FPT[FPTNum], Corr, (Int_t)(centLow*.5), (Int_t)((centHi + 1)*.5), fileTag);
  else
    title = Form("%s%sImbProjA%s%s%s_PP_%s_h", gorr, algType[setNum], CNCR, FPT[FPTNum], Corr, fileTag);

  Float_t xArr[11] = {0.0001, 0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 1.999};

  //edit here

  TH1F* imbDelRHist_p = new TH1F("imbDelRHist_p", "imbDelRHist_p", 10, xArr);

  imbDelRHist_p->GetXaxis()->SetLimits(0.00, 2.00);
  niceTH1(imbDelRHist_p, histHi, histLow, 505, 403);

  TH1F* getHist_p;


  //Throw dis shit in da loop

  TCut setCut = makeSetCut(setNum);
  TCut centCut = "";
  if(sType == kHIDATA || sType == kHIMC) centCut = makeCentCut(centLow, centHi);
  TCut etaCut = makeEtaCut(setNum, 0.5);
  TCut phiCut = makeDelPhiCut(setNum, 5*TMath::Pi()/6);
  TCut asymmCut = makeAsymmCut(setNum, 0.0, 1.0);
  if(!strcmp(CNCR, "RU")){
    asymmCut = makeAsymmCut(setNum, 0.22, 1.0);
    if(isPercent){
      Float_t dummArr[5] = {0, 0, 0, 0, 0};
      setAsymmPercBins(dummArr, centLow, centHi, sType);
      asymmCut = makeAsymmCut(setNum, dummArr[2], 1.0);
    }
  }
  else if(!strcmp(CNCR, "RD")){
    asymmCut = makeAsymmCut(setNum, 0.0, 0.22);
    if(isPercent){
      Float_t dummArr[5] = {0, 0, 0, 0, 0};
      setAsymmPercBins(dummArr, centLow, centHi, sType);
      asymmCut = makeAsymmCut(setNum, 0.0, dummArr[2]);
    }
  }

  TCut jetLCut = Form("AlgJtPt[%d][0] > %d", setNum, leadJtCut);
  TCut jetSLCut = Form("AlgJtPt[%d][1] > %d", setNum, subLeadJtCut);
  TCut pthat = "pthat > 80";

  TCut trkCut = "";
  if(isHighPtTrk)
    trkCut = Form("AlgJtTrkMax[%d][0] > 12 && AlgJtTrkMax[%d][1] > 12", setNum, setNum);

  const char* name1[10] = {"0(10000, -10000, 10000)", "1(10000, -10000, 10000)", "2(10000, -10000, 10000)", "3(10000, -10000, 10000)", "4(10000, -10000, 10000)", "5(10000, -10000, 10000)", "6(10000, -10000, 10000)", "7(10000, -1000, 10000)", "8(10000, -1000, 10000)", "9(10000, -1000, 10000)"};
  const char* name2[10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};

  for(Int_t rIter = 0; rIter < 10; rIter++){
    TString var = Form("%sAlgImbProjAR[%d][%d][%d]", gorr, setCorrNum, FPTNum, rIter);

    if(montecarlo)
      anaTree_p->Project(name1[rIter], var, setCut && centCut && etaCut && phiCut && jetLCut && jetSLCut && asymmCut && pthat && trkCut, "");
    else
      anaTree_p->Project(name1[rIter], var, setCut && centCut && etaCut && phiCut && jetLCut && jetSLCut && asymmCut && trkCut, "");

    getHist_p = (TH1F*)inFile_p->Get(name2[rIter]);

    imbDelRHist_p->SetBinContent(rIter+1, getHist_p->GetMean());
    imbDelRHist_p->SetBinError(rIter+1, getHist_p->GetMeanError());
  }

  outFile_p = new TFile(outName, "UPDATE");
  imbDelRHist_p->Write(title);
  outFile_p->Close();

  delete outFile_p;
  delete imbDelRHist_p;

  return;
}

void makeMultDiffHist(TTree* anaTree_p, const char* outName, Int_t setNum, Int_t centLow, Int_t centHi, Float_t histLow, Float_t histHi, sampleType sType = kHIDATA, const char* Tight = "", Bool_t isHighPtTrk = false)
{
  inFile_p->cd();

  Bool_t montecarlo = false;
  if(sType == kHIMC || sType == kPPMC || sType == kPAMC)
    montecarlo = true;

  const char* title;
  
  if(sType == kHIDATA || sType == kHIMC)
    title = Form("%sMultA%s_%d%d_%s_h", algType[setNum], Tight, (Int_t)(centLow*.5), (Int_t)((centHi + 1)*.5), fileTag);
  else
    title = Form("%sMultA%s_PP_%s_h", algType[setNum], Tight, fileTag);

  Float_t xArr[5] = {.0001, .11, .22, .33, .4999};
  Float_t xArrTight[9] = {.0001, .055, .11, .165, .22, .275, .33, .415, .4999};
  TH1F* multHist_p;

  if(!strcmp(Tight, ""))
    multHist_p = new TH1F("multHist_p", "multHist_p", 4, xArr);
  else
    multHist_p = new TH1F("multHist_p", "multHist_p", 8, xArrTight);

  multHist_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTH1(multHist_p, histHi, histLow, 505, 508);

  TH1F* getHist_p;

  TString var = Form("AlgJtMult[%d][1] - AlgJtMult[%d][0]", setNum, setNum);

  TCut setCut = makeSetCut(setNum);
  TCut centCut = "";
  if(sType == kHIDATA || sType == kHIMC) centCut = makeCentCut(centLow, centHi);
  TCut etaCut = makeEtaCut(setNum, 1.6);
  TCut phiCut = makeDelPhiCut(setNum, 5*TMath::Pi()/6);

  TCut jetLCut = Form("AlgJtPt[%d][0] > %d", setNum, leadJtCut);
  TCut jetSLCut = Form("AlgJtPt[%d][1] > %d", setNum, subLeadJtCut);
  TCut pthat = "";
  if(montecarlo) pthat = "pthat > 80";

  TCut trkCut = "";
  if(isHighPtTrk)
    trkCut = Form("AlgJtTrkMax[%d][0] > 12 && AlgJtTrkMax[%d][1] > 12", setNum, setNum);

  const char* name1[8] = {"0(10000, -10000, 10000)", "1(10000, -10000, 10000)", "2(10000, -10000, 10000)", "3(10000, -10000, 10000)", "4(10000, -10000, 10000)", "5(10000, -10000, 10000)", "6(10000, -10000, 10000)", "7(10000, -10000, 10000)"};
  const char* name2[8] = {"0", "1", "2", "3", "4", "5", "6", "7"};

  Float_t asymmBins[5] = {.00, .11, .22, .33, 1.00};
  Float_t asymmBinsTight[9] = {.00, .055, .11, .165, .22, .275, .33, .415, 1.00};

  Int_t nBins = 4;
  if(!strcmp(Tight, "Tight"))
    nBins = 8;

  for(Int_t binIter = 0; binIter < nBins; binIter++){
    TCut asymmCut = "";

    if(!strcmp(Tight, ""))
      asymmCut = makeAsymmCut(setNum, asymmBins[binIter], asymmBins[binIter + 1]);
    else
      asymmCut = makeAsymmCut(setNum, asymmBinsTight[binIter], asymmBinsTight[binIter + 1]);

      anaTree_p->Project(name1[binIter], var, setCut && centCut && etaCut && phiCut && jetLCut && jetSLCut && asymmCut && pthat && trkCut, "");

      getHist_p = (TH1F*)inFile_p->Get(name2[binIter]);

      multHist_p->SetBinContent(binIter + 1, getHist_p->GetMean());
      multHist_p->SetBinError(binIter + 1, getHist_p->GetMeanError());
  }

  outFile_p = new TFile(outName, "UPDATE");
  multHist_p->Write(title);
  outFile_p->Close();

  delete outFile_p;
  delete multHist_p;

  return;
}

void makeDiJetHists(const char* inName, const char* outName, sampleType sType = kHIDATA, Bool_t isPercent = false, Bool_t isHighPtTrk = false)
{
  TH1::SetDefaultSumw2();

  Bool_t montecarlo = false;
  if(sType == kHIMC || sType == kPPMC || sType == kPAMC)
    montecarlo = true;

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
  const char* CNCR[6] = {"", "C", "NC", "R", "RD", "RU"};
  const char* Tight[2] = {"", "Tight"};

  Int_t centBins = 6;
  Int_t centLow[6] = {0, 20, 60, 100, 0, 60};
  Int_t centHi[6] = {19, 59, 99, 199, 59, 199};

  if(sType == kPPMC || sType == kPPDATA) centBins = 1;

  for(Int_t algIter = 1; algIter < jetAlgMax; algIter++){
    std::cout << "Algorithm: " << algType[algIter] << std::endl;

    for(Int_t centIter = 0; centIter < centBins; centIter++){
      for(Int_t corrIter = 0; corrIter < 2; corrIter++){
	for(Int_t FPTIter = 0; FPTIter < 6; FPTIter++){
	  for(Int_t CNCRIter = 0; CNCRIter < 6; CNCRIter++){

	    if(CNCRIter < 3){

	      for(Int_t tightIter = 0; tightIter < 2; tightIter++){

		if((CNCRIter == 0 && centIter < 4) || (CNCRIter != 0 && centIter >= 4) || sType == kPPMC || sType == kPPDATA){
		  makeImbAsymmHist(anaTree_p, outName, "r", algIter, CNCR[CNCRIter], FPTIter, centLow[centIter], centHi[centIter], -60., 59.999, corr[corrIter], Tight[tightIter], sType, isPercent, isHighPtTrk);
		  
		  if(corrIter == 0 && FPTIter == 0 && CNCRIter == 0)
		    makeMultDiffHist(anaTree_p, outName, algIter, centLow[centIter], centHi[centIter], .0001, 39.9999, sType, Tight[tightIter], false);      

		  if(montecarlo && corrIter == 0)
		    makeImbAsymmHist(anaTree_p, outName, "g", algIter, CNCR[CNCRIter], FPTIter, centLow[centIter], centHi[centIter], -60., 59.999, corr[corrIter], Tight[tightIter], sType, isPercent, isHighPtTrk);
		}
	      }
	    }
	    else{
	      if(centIter >= 4 || sType == kPPMC || sType == kPPDATA){
		makeImbDelRHist(anaTree_p, outName, "r", algIter, CNCR[CNCRIter], FPTIter, centLow[centIter], centHi[centIter], -40., 19.999, corr[corrIter], sType, isPercent, isHighPtTrk);
		
		if(montecarlo && corrIter == 0)
		  makeImbDelRHist(anaTree_p, outName, "r", algIter, CNCR[CNCRIter], FPTIter, centLow[centIter], centHi[centIter], -40., 19.999, corr[corrIter], sType, isPercent, isHighPtTrk);
	      }
	    }
	  }
	}
      }
    }
  }
  
  inFile_p->Close();
  delete inFile_p;

  return;
}
