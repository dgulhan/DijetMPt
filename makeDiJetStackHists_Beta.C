//=============================================                      
// Author: Chris McGinn
//                     
// DiJet Histogram Maker, Missing Pt
//                                                                            
//=============================================     

#include <string>

#include "TTree.h"
#include "TDatime.h"
#include "TMath.h"
#include "TFile.h"
#include "diJetFileTag.h"
#include "/net/hisrv0001/home/cfmcginn/DijetMPt/CMSSW_5_3_20/src/DijetHists/cfmDijetMeanTree.h"

#include <fstream>

TFile* inFile_p = 0;
TFile* outFile_p = 0;

const std::string FPT[6] = {"0_1", "1_2", "2_4", "4_8", "8_100", "F"};

const Float_t trkMaxCut = 8.;

const Float_t loose[5] = {0.00, .11, .22, .33, 1.00};
const Float_t tight[9] = {0.00, .055, .11, .165, .22, .275, .33, .415, 1.00};
const Float_t niceNumCNC[4] = {59.999, -60, 505, 406};
const Float_t niceNumR[4] = {19.999, -40, 505, 403};
const Float_t niceNumMult[4] = {9.999, -10, 505, 403};

const Int_t evtCutPos[7] = {2, 2, 2, 6, 6, 6, 6};
const Bool_t inVenn = false;
const Bool_t outVenn = false;

const Float_t midRapCut = 0.5;
const Float_t midRapCut2 = 0.8;

Int_t retBinNumber(const std::string Tight)
{
  if(!strcmp(Tight.c_str(), "")) return 4;
  else return 8;
}


void getBinArr(const Int_t nBins, Float_t bins[], const Float_t initArray[])
{
  for(Int_t iter = 0; iter < nBins + 1; iter++){
    bins[iter] = initArray[iter];
  }
 
  bins[0] = .0001;
  bins[nBins] = .4999;

  return;
}


void CleanHist(TH1* inHist_p)
{
  delete inHist_p;
  inHist_p = 0;
  
  return;
}


std::string getCentString(sampleType sType, Int_t centLow, Int_t centHi)
{
  if(isHI(sType)) return Form("%d%d", (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5));
  else return "PP";
}


void InitHist(TH1F* inHist_p[6])
{
  for(Int_t iter = 0; iter < 6; iter++){
    inHist_p[iter] = 0;
  }
  return;
}


void BookHist(TH1F* inHist_p[6], const std::string gR, const std::string inAlgType, const std::string multProj, const std::string CNCR, const std::string Corr, const std::string centString, const Int_t nBins, Float_t xArr[], Float_t upBound, const Float_t niceNum[4])
{
  for(Int_t iter = 0; iter < 6; iter++){
    const std::string title = Form("%s%s%sA%s%s%s_%s_%s_h", gR.c_str(), inAlgType.c_str(), multProj.c_str(), CNCR.c_str(), FPT[iter].c_str(), Corr.c_str(), centString.c_str(), fileTag.c_str());
    inHist_p[iter] = new TH1F(title.c_str(), title.c_str(), nBins, xArr);

    inHist_p[iter]->GetXaxis()->SetLimits(0.00, upBound);
    niceTH1(inHist_p[iter], niceNum[0], niceNum[1], niceNum[2], niceNum[3]);
  }
  return;
}

void makeImbAHist(TTree* inTree_p, const std::string outName, const std::string gR, Int_t setNum, Int_t centLow, Int_t centHi, const std::string Corr = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

  Float_t xArr[nAjBins+1];
  getBinArr(nAjBins, xArr, loose);
  
  Int_t setCorrNum = setNum;
  if(!strcmp("Corr", Corr.c_str())) setCorrNum = setNum + nSumAlg;

  const std::string centString = getCentString(sType, centLow, centHi);
  Int_t ajCentNum = -1;
  if(hi){
    for(Int_t ajCentIter = 0; ajCentIter < nAjCentBins; ajCentIter++){
      if(centLow < ajCentBins[ajCentIter]){
	ajCentNum = ajCentIter;
	break;
      }
    }
  }
  else ajCentNum = 0;

  TH1F* rImbAHist_p[6];  
  InitHist(rImbAHist_p);
  BookHist(rImbAHist_p, gR, algType[setNum], "Proj", "", Corr, centString, nAjBins, xArr, 0.50, niceNumCNC);

  Int_t outMult[nPtBins][nAjBins];
  Float_t outMean[nPtBins][nAjBins];
  Float_t outWeight[nPtBins][nAjBins];
  Float_t outWeight2[nPtBins][nAjBins];
  Float_t outSig[nPtBins][nAjBins];

  for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
    for(Int_t ajIter = 0; ajIter < nAjBins; ajIter++){
      outMult[ptIter][ajIter] = 0;
      outMean[ptIter][ajIter] = 0;
      outWeight[ptIter][ajIter] = 0;
      outWeight2[ptIter][ajIter] = 0;
      outSig[ptIter][ajIter] = 0;
    }
  }

  Int_t nentries = inTree_p->GetEntries();

  for(Int_t jEntry = 0; jEntry < nentries; jEntry++){
    inTree_p->GetEntry(jEntry);

    for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
      for(Int_t ajIter = 0; ajIter < nAjBins; ajIter++){
	if(!strcmp(gR.c_str(), "r")){
	  outMult[ptIter][ajIter] += rProjA_MULT_[setCorrNum][ptIter][ajCentNum][ajIter];
	  outMean[ptIter][ajIter] += rProjA_MEAN_[setCorrNum][ptIter][ajCentNum][ajIter]*rProjA_WEIGHT_[setCorrNum][ptIter][ajCentNum][ajIter];
	  outWeight[ptIter][ajIter] += rProjA_WEIGHT_[setCorrNum][ptIter][ajCentNum][ajIter];
	  outWeight2[ptIter][ajIter] += rProjA_WEIGHT2_[setCorrNum][ptIter][ajCentNum][ajIter];
	}
	else{
	  outMult[ptIter][ajIter] += gProjA_MULT_[setNum][ptIter][ajCentNum][ajIter];
	  outMean[ptIter][ajIter] += gProjA_MEAN_[setNum][ptIter][ajCentNum][ajIter]*gProjA_WEIGHT_[setNum][ptIter][ajCentNum][ajIter];
	  outWeight[ptIter][ajIter] += gProjA_WEIGHT_[setNum][ptIter][ajCentNum][ajIter];
	  outWeight2[ptIter][ajIter] += gProjA_WEIGHT2_[setNum][ptIter][ajCentNum][ajIter];
	}
      }
    }
  }

  for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
    for(Int_t ajIter = 0; ajIter < nAjBins; ajIter++){
      outMean[ptIter][ajIter] /= outWeight[ptIter][ajIter];
    }
  }

  for(Int_t jEntry = 0; jEntry < nentries; jEntry++){
    inTree_p->GetEntry(jEntry);

    for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
      for(Int_t ajIter = 0; ajIter < nAjBins; ajIter++){

	if(!strcmp(gR.c_str(), "r")) outSig[ptIter][ajIter] += rProjA_WEIGHT_[setCorrNum][ptIter][ajCentNum][ajIter]*(rProjA_SIG_[setCorrNum][ptIter][ajCentNum][ajIter] + TMath::Power(rProjA_MEAN_[setCorrNum][ptIter][ajCentNum][ajIter], 2));
	else outSig[ptIter][ajIter] += gProjA_WEIGHT_[setNum][ptIter][ajCentNum][ajIter]*(gProjA_SIG_[setNum][ptIter][ajCentNum][ajIter] + TMath::Power(gProjA_MEAN_[setNum][ptIter][ajCentNum][ajIter], 2));
      }
    }
  }

  for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
    for(Int_t ajIter = 0; ajIter < nAjBins; ajIter++){
      outSig[ptIter][ajIter] /= outWeight[ptIter][ajIter];
      outSig[ptIter][ajIter] -= TMath::Power(outMean[ptIter][ajIter], 2);
      outSig[ptIter][ajIter] = TMath::Sqrt(outSig[ptIter][ajIter]/(outWeight[ptIter][ajIter]*outWeight[ptIter][ajIter]/outWeight2[ptIter][ajIter]));
    }
  }

  for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
    for(Int_t ajIter = 0; ajIter < nAjBins; ajIter++){
      rImbAHist_p[ptIter]->SetBinContent(ajIter+1, outMean[ptIter][ajIter]);
      rImbAHist_p[ptIter]->SetBinError(ajIter+1, outSig[ptIter][ajIter]);
    }
  }

  outFile_p = new TFile(Form("%s.root", outName.c_str()), "UPDATE");
  std::cout << outName << std::endl;

  for(Int_t iter = 0; iter < 6; iter++){
    rImbAHist_p[iter]->Write("", TObject::kOverwrite);
  }

  outFile_p->Close();
  delete outFile_p;
  outFile_p = 0;
  
  for(Int_t iter = 0; iter < 6; iter++){
    CleanHist(rImbAHist_p[iter]);
  }

  return;
}




void makeImbARHist(TTree* inTree_p, const std::string outName, const std::string gR, Int_t setNum, Int_t centLow, Int_t centHi, const std::string Corr = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

  const Int_t nBins = 10;
  Float_t xArr[nBins+1] = {0.0001, 0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 1.9999};
  
  Int_t setCorrNum = setNum;
  if(!strcmp("Corr", Corr.c_str())) setCorrNum = setNum + nSumAlg;

  const std::string centString = getCentString(sType, centLow, centHi);
  Int_t rCentNum = 0;
  if(hi){
    for(Int_t rCentIter = 0; rCentIter < nRCentBins; rCentIter++){
      if(centLow < rCentBins[rCentIter]){
	rCentNum = rCentIter;
	break;
      }
    }
  }

  TH1F* rImbARHist_p[6];  
  InitHist(rImbARHist_p);
  BookHist(rImbARHist_p, gR, algType[setNum], "Proj", "R", Corr, centString, nBins, xArr, 2.00, niceNumR);

  TH1F* rImbARDHist_p[6];  
  InitHist(rImbARDHist_p);
  BookHist(rImbARDHist_p, gR, algType[setNum], "Proj", "RD", Corr, centString, nBins, xArr, 2.00, niceNumR);

  TH1F* rImbARUHist_p[6];  
  InitHist(rImbARUHist_p);
  BookHist(rImbARUHist_p, gR, algType[setNum], "Proj", "RU", Corr, centString, nBins, xArr, 2.00, niceNumR);


  Int_t outRMult[nPtBins][nRBins];
  Float_t outRMean[nPtBins][nRBins];
  Float_t outRWeight[nPtBins][nRBins];
  Float_t outRWeight2[nPtBins][nRBins];
  Float_t outRSig[nPtBins][nRBins];

  Int_t outRDMult[nPtBins][nRBins];
  Float_t outRDMean[nPtBins][nRBins];
  Float_t outRDWeight[nPtBins][nRBins];
  Float_t outRDWeight2[nPtBins][nRBins];
  Float_t outRDSig[nPtBins][nRBins];

  Int_t outRUMult[nPtBins][nRBins];
  Float_t outRUMean[nPtBins][nRBins];
  Float_t outRUWeight[nPtBins][nRBins];
  Float_t outRUWeight2[nPtBins][nRBins];
  Float_t outRUSig[nPtBins][nRBins];

  for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
    for(Int_t rIter = 0; rIter < nRBins; rIter++){
      outRMult[ptIter][rIter] = 0;
      outRMean[ptIter][rIter] = 0;
      outRWeight[ptIter][rIter] = 0;
      outRWeight2[ptIter][rIter] = 0;
      outRSig[ptIter][rIter] = 0;

      outRDMult[ptIter][rIter] = 0;
      outRDMean[ptIter][rIter] = 0;
      outRDWeight[ptIter][rIter] = 0;
      outRDWeight2[ptIter][rIter] = 0;
      outRDSig[ptIter][rIter] = 0;

      outRUMult[ptIter][rIter] = 0;
      outRUMean[ptIter][rIter] = 0;
      outRUWeight[ptIter][rIter] = 0;
      outRUWeight2[ptIter][rIter] = 0;
      outRUSig[ptIter][rIter] = 0;
    }
  }

  Int_t nentries = inTree_p->GetEntries();

  for(Int_t jEntry = 0; jEntry < nentries; jEntry++){
    inTree_p->GetEntry(jEntry);

    for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
      for(Int_t rIter = 0; rIter < nRBins; rIter++){
	if(!strcmp(gR.c_str(), "r")){
	  outRMult[ptIter][rIter] += rProjAR_MULT_[setCorrNum][ptIter][rCentNum][rIter];
	  outRMean[ptIter][rIter] += rProjAR_MEAN_[setCorrNum][ptIter][rCentNum][rIter]*rProjAR_WEIGHT_[setCorrNum][ptIter][rCentNum][rIter];
	  outRWeight[ptIter][rIter] += rProjAR_WEIGHT_[setCorrNum][ptIter][rCentNum][rIter];
	  outRWeight2[ptIter][rIter] += rProjAR_WEIGHT2_[setCorrNum][ptIter][rCentNum][rIter];

	  outRDMult[ptIter][rIter] += rProjARD_MULT_[setCorrNum][ptIter][rCentNum][rIter];
	  outRDMean[ptIter][rIter] += rProjARD_MEAN_[setCorrNum][ptIter][rCentNum][rIter]*rProjARD_WEIGHT_[setCorrNum][ptIter][rCentNum][rIter];
	  outRDWeight[ptIter][rIter] += rProjARD_WEIGHT_[setCorrNum][ptIter][rCentNum][rIter];
	  outRDWeight2[ptIter][rIter] += rProjARD_WEIGHT2_[setCorrNum][ptIter][rCentNum][rIter];

	  outRUMult[ptIter][rIter] += rProjARU_MULT_[setCorrNum][ptIter][rCentNum][rIter];
	  outRUMean[ptIter][rIter] += rProjARU_MEAN_[setCorrNum][ptIter][rCentNum][rIter]*rProjARU_WEIGHT_[setCorrNum][ptIter][rCentNum][rIter];
	  outRUWeight[ptIter][rIter] += rProjARU_WEIGHT_[setCorrNum][ptIter][rCentNum][rIter];
	  outRUWeight2[ptIter][rIter] += rProjARU_WEIGHT2_[setCorrNum][ptIter][rCentNum][rIter];
	}
	else{
	  outRMult[ptIter][rIter] += gProjAR_MULT_[setNum][ptIter][rCentNum][rIter];
	  outRMean[ptIter][rIter] += gProjAR_MEAN_[setNum][ptIter][rCentNum][rIter]*gProjAR_WEIGHT_[setNum][ptIter][rCentNum][rIter];
	  outRWeight[ptIter][rIter] += gProjAR_WEIGHT_[setNum][ptIter][rCentNum][rIter];
	  outRWeight2[ptIter][rIter] += gProjAR_WEIGHT2_[setNum][ptIter][rCentNum][rIter];

	  outRDMult[ptIter][rIter] += gProjARD_MULT_[setNum][ptIter][rCentNum][rIter];
	  outRDMean[ptIter][rIter] += gProjARD_MEAN_[setNum][ptIter][rCentNum][rIter]*gProjARD_WEIGHT_[setNum][ptIter][rCentNum][rIter];
	  outRDWeight[ptIter][rIter] += gProjARD_WEIGHT_[setNum][ptIter][rCentNum][rIter];
	  outRDWeight2[ptIter][rIter] += gProjARD_WEIGHT2_[setNum][ptIter][rCentNum][rIter];

	  outRUMult[ptIter][rIter] += gProjARU_MULT_[setNum][ptIter][rCentNum][rIter];
	  outRUMean[ptIter][rIter] += gProjARU_MEAN_[setNum][ptIter][rCentNum][rIter]*gProjARU_WEIGHT_[setNum][ptIter][rCentNum][rIter];
	  outRUWeight[ptIter][rIter] += gProjARU_WEIGHT_[setNum][ptIter][rCentNum][rIter];
	  outRUWeight2[ptIter][rIter] += gProjARU_WEIGHT2_[setNum][ptIter][rCentNum][rIter];
	}
      }
    }
  }

  for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
    for(Int_t rIter = 0; rIter < nRBins; rIter++){
      outRMean[ptIter][rIter] /= outRWeight[ptIter][rIter];

      outRDMean[ptIter][rIter] /= outRDWeight[ptIter][rIter];

      outRUMean[ptIter][rIter] /= outRUWeight[ptIter][rIter];
    }
  }

  for(Int_t jEntry = 0; jEntry < nentries; jEntry++){
    inTree_p->GetEntry(jEntry);

    for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
      for(Int_t rIter = 0; rIter < nRBins; rIter++){
	outRSig[ptIter][rIter] += rProjAR_WEIGHT_[setCorrNum][ptIter][rCentNum][rIter]*(rProjAR_SIG_[setCorrNum][ptIter][rCentNum][rIter] + TMath::Power(rProjAR_MEAN_[setCorrNum][ptIter][rCentNum][rIter], 2));

	outRDSig[ptIter][rIter] += rProjARD_WEIGHT_[setCorrNum][ptIter][rCentNum][rIter]*(rProjARD_SIG_[setCorrNum][ptIter][rCentNum][rIter] + TMath::Power(rProjARD_MEAN_[setCorrNum][ptIter][rCentNum][rIter], 2));

	outRUSig[ptIter][rIter] += rProjARU_WEIGHT_[setCorrNum][ptIter][rCentNum][rIter]*(rProjARU_SIG_[setCorrNum][ptIter][rCentNum][rIter] + TMath::Power(rProjARU_MEAN_[setCorrNum][ptIter][rCentNum][rIter], 2));
      }
    }
  }

  for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
    for(Int_t rIter = 0; rIter < nRBins; rIter++){
      outRSig[ptIter][rIter] /= outRWeight[ptIter][rIter];
      outRSig[ptIter][rIter] -= TMath::Power(outRMean[ptIter][rIter], 2);
      outRSig[ptIter][rIter] = TMath::Sqrt(outRSig[ptIter][rIter]/(outRWeight[ptIter][rIter]*outRWeight[ptIter][rIter]/outRWeight2[ptIter][rIter]));

      outRDSig[ptIter][rIter] /= outRDWeight[ptIter][rIter];
      outRDSig[ptIter][rIter] -= TMath::Power(outRDMean[ptIter][rIter], 2);
      outRDSig[ptIter][rIter] = TMath::Sqrt(outRDSig[ptIter][rIter]/(outRDWeight[ptIter][rIter]*outRDWeight[ptIter][rIter]/outRDWeight2[ptIter][rIter]));

      outRUSig[ptIter][rIter] /= outRUWeight[ptIter][rIter];
      outRUSig[ptIter][rIter] -= TMath::Power(outRUMean[ptIter][rIter], 2);
      outRUSig[ptIter][rIter] = TMath::Sqrt(outRUSig[ptIter][rIter]/(outRUWeight[ptIter][rIter]*outRUWeight[ptIter][rIter]/outRUWeight2[ptIter][rIter]));
    }
  }

  for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
    for(Int_t rIter = 0; rIter < nRBins; rIter++){
      rImbARHist_p[ptIter]->SetBinContent(rIter+1, outRMean[ptIter][rIter]);
      rImbARHist_p[ptIter]->SetBinError(rIter+1, outRSig[ptIter][rIter]);

      rImbARDHist_p[ptIter]->SetBinContent(rIter+1, outRDMean[ptIter][rIter]);
      rImbARDHist_p[ptIter]->SetBinError(rIter+1, outRDSig[ptIter][rIter]);

      rImbARUHist_p[ptIter]->SetBinContent(rIter+1, outRUMean[ptIter][rIter]);
      rImbARUHist_p[ptIter]->SetBinError(rIter+1, outRUSig[ptIter][rIter]);
    }
  }

  outFile_p = new TFile(Form("%s.root", outName.c_str()), "UPDATE");
  std::cout << outName << std::endl;

  for(Int_t iter = 0; iter < 6; iter++){
    rImbARHist_p[iter]->Write("", TObject::kOverwrite);

    rImbARDHist_p[iter]->Write("", TObject::kOverwrite);

    rImbARUHist_p[iter]->Write("", TObject::kOverwrite);
  }

  outFile_p->Close();
  delete outFile_p;
  outFile_p = 0;
  
  for(Int_t iter = 0; iter < 6; iter++){
    CleanHist(rImbARHist_p[iter]);
    CleanHist(rImbARDHist_p[iter]);
    CleanHist(rImbARUHist_p[iter]);
  }

  return;
}




int makeDiJetStackHists(const std::string inName, sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  TH1::SetDefaultSumw2();
  Bool_t montecarlo = isMonteCarlo(sType);

  setFileTag(inName);
 
  std::string outName = inName;
  const std::string cutString[2] = {"MeanTree", ".root"};
  const std::string repString[2] = {"StackHists", ""};

  std::cout << "Replace string" << std::endl;

  for(Int_t iter = 0; iter < 2; iter++){
    std::size_t strIndex = outName.find(cutString[iter]);
    if(!(strIndex == std::string::npos)){
      outName.replace(strIndex, cutString[iter].length(), repString[iter]);
    }
  }

  inFile_p = new TFile(inName.c_str(), "READ");
  GetDiJetMeanTree(inFile_p, sType);

  std::cout << "YOLO" << std::endl;

  const std::string Corr[2] = {"", "Corr"};

  const char* rStr[39] = {"R", "RD", "RU", "Eta", "EtaD", "EtaU", "Phi", "PhiD", "PhiU", "RCut", "RCutD", "RCutU", "RCutEta", "RCutEtaD", "RCutEtaU", "RCutPhi", "RCutPhiD", "RCutPhiU", "EtaCut", "EtaCutD", "EtaCutU", "PhiCut", "PhiCutD", "PhiCutU", "RFOR", "RFORD", "RFORU", "RFORMID", "RFORMIDD", "RFORMIDU", "RFORFOR", "RFORFORD", "RFORFORU", "RMID", "RMIDD", "RMIDU", "R2", "R2U", "R2D"};

  if(montecarlo) trackTreeMean_p->AddFriend(genTreeMean_p);

  for(Int_t algIter = 0; algIter < nSumAlg; algIter++){
    for(Int_t corrIter = 1; corrIter < 2; corrIter++){

      std::cout << algType[algIter] << std::endl;

      if(algIter < nSumAlg - 4){
	makeImbAHist(trackTreeMean_p, outName, "r", algIter, 0, 19, Corr[corrIter], sType, isHighPtTrk);
	
	if(isHI(sType)){
	  makeImbAHist(trackTreeMean_p, outName, "r", algIter, 20, 59, Corr[corrIter], sType, isHighPtTrk);
	  makeImbAHist(trackTreeMean_p, outName, "r", algIter, 60, 99, Corr[corrIter], sType, isHighPtTrk);
	  makeImbAHist(trackTreeMean_p, outName, "r", algIter, 100, 199, Corr[corrIter], sType, isHighPtTrk);
	}

	makeImbARHist(trackTreeMean_p, outName, "r", algIter, 0, 59, Corr[corrIter], sType, isHighPtTrk);
	if(isHI(sType)) makeImbARHist(trackTreeMean_p, outName, "r", algIter, 60, 199, Corr[corrIter], sType, isHighPtTrk);
      }

      if(montecarlo && algIter < nSumAlg){
	makeImbAHist(trackTreeMean_p, outName, "g", algIter, 0, 19, Corr[corrIter], sType, isHighPtTrk);
	
	if(isHI(sType)){
	  makeImbAHist(trackTreeMean_p, outName, "g", algIter, 20, 59, Corr[corrIter], sType, isHighPtTrk);
	  makeImbAHist(trackTreeMean_p, outName, "g", algIter, 60, 99, Corr[corrIter], sType, isHighPtTrk);
	  makeImbAHist(trackTreeMean_p, outName, "g", algIter, 100, 199, Corr[corrIter], sType, isHighPtTrk);
	}

	makeImbARHist(trackTreeMean_p, outName, "g", algIter, 0, 59, Corr[corrIter], sType, isHighPtTrk);
	if(isHI(sType)) makeImbARHist(trackTreeMean_p, outName, "g", algIter, 60, 199, Corr[corrIter], sType, isHighPtTrk);
      }

    }
  }

  return 0;
}
