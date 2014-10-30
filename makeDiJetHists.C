//=============================================                      
// Author: Chris McGinn
//                     
// DiJet Histogram Maker, Missing Pt
//                                                                            
//=============================================     

#include <string>

#include "TTree.h"
#include "TDatime.h"
#include "TFile.h"
#include "diJetFileTag.h"
#include "/net/hisrv0001/home/cfmcginn/DijetMPt/CMSSW_5_3_20/src/DijetAnalysisSkim/cfmDiJetAnaSkim.h"
#include "/net/hisrv0001/home/cfmcginn/DijetMPt/CMSSW_5_3_20/src/DijetInitialSkim/cfmVectFunc.h"

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


void getBinArr(const Int_t nBins, Float_t bins[], Float_t cuts[], const Float_t initArray[])
{
  for(Int_t iter = 0; iter < nBins + 1; iter++){
    bins[iter] = initArray[iter];
    cuts[iter] = initArray[iter];
  }
 
  bins[0] = .0001;
  bins[nBins] = .4999;

  return;
}


void CleanHist(TH1* inHist_p)
{
  if(inHist_p!=0){
    delete inHist_p;
    inHist_p = 0;
  }
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


Bool_t isEventCut(Int_t setNum, sampleType sType, Int_t centLow, Int_t centHi, Bool_t isHighPtTrk = false)
{
  if(!eventSet_[setNum]) return true;
  if(AlgJtPt_[setNum][0] < leadJtPtCut[setNum] || AlgJtPt_[setNum][1] < subLeadJtPtCut[setNum]) return true;
  if(TMath::Abs(AlgJtEta_[setNum][0]) > 1.6 || TMath::Abs(AlgJtEta_[setNum][1]) > 1.6) return true;
  if(AlgJtDelPhi12_[setNum] < 5.0*TMath::Pi()/6.0) return true;


  if(inVenn){
    if(!eventSet_[evtCutPos[setNum]]) return true;
    if(AlgJtPt_[evtCutPos[setNum]][0] < leadJtPtCut[setNum] || AlgJtPt_[evtCutPos[setNum]][1] < subLeadJtPtCut[setNum]) return true;
    if(TMath::Abs(AlgJtEta_[evtCutPos[setNum]][0]) > 1.6 || TMath::Abs(AlgJtEta_[evtCutPos[setNum]][1]) > 1.6) return true;
    if(AlgJtDelPhi12_[evtCutPos[setNum]] < 5.0*TMath::Pi()/6.0) return true;
  }
  else if(outVenn){
    if(eventSet_[evtCutPos[setNum]]){
      if(AlgJtPt_[evtCutPos[setNum]][0] > leadJtPtCut[setNum] && AlgJtPt_[evtCutPos[setNum]][1] > subLeadJtPtCut[setNum]){
	if(TMath::Abs(AlgJtEta_[evtCutPos[setNum]][0]) < 1.6 && TMath::Abs(AlgJtEta_[evtCutPos[setNum]][1]) < 1.6){
	  if(AlgJtDelPhi12_[evtCutPos[setNum]] > 5.0*TMath::Pi()/6.0) return true;
	}
      }
    }
  }

  if(isHI(sType)){
     if(hiBin_ < centLow || hiBin_ > centHi) return true;
  }

  if(isHighPtTrk && (AlgJtTrkMax_[setNum][0] < trkMaxCut || AlgJtTrkMax_[setNum][1] < trkMaxCut)) return true;

  if(isMonteCarlo(sType) && pthat_ < 80) return true;

  return false;
}

/*
void makeMultAHist(TTree* anaTree_p, const std::string outName, Int_t setNum, Int_t centLow, Int_t centHi, const std::string Corr = "", const std::string Tight = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  Int_t setCorrNum = setNum;
  if(!strcmp("Corr", Corr.c_str()))
    setCorrNum = setNum + 3;

  const Int_t nBins = retBinNumber(Tight);
  Float_t xArr[nBins+1];
  Float_t xArrCut[nBins+1];

  if(nBins == 4) getBinArr(nBins, xArr, xArrCut, loose);
  else if(nBins == 8) getBinArr(nBins, xArr, xArrCut, tight);

  const std::string centString = getCentString(sType, centLow, centHi);

  std::string title = Form("r%sMultA%s%s_%s_%s_h", algType[setNum], Corr.c_str(), Tight.c_str(), centString.c_str(), fileTag);
  TH1F* rMultAHist_p = 0;
  rMultAHist_p = new TH1F(Form("rMultAHist_p"), Form("rMultAHist_p"), nBins, xArr);
  rMultAHist_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTH1(rMultAHist_p, niceNumCNC[0], niceNumCNC[1], niceNumCNC[2], niceNumCNC[3]);
  std::vector<Float_t>* mean_rMultA_p[nBins];
 
  for(Int_t iter = 0; iter < nBins; iter++){
    mean_rMultA_p[iter] = new std::vector<Float_t>;
  }

  return;
}
*/

void makeImbAHist(TTree* anaTree_p, const std::string outName, Int_t setNum, Int_t centLow, Int_t centHi, const std::string Corr = "", const std::string Tight = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);

  Int_t setCorrNum = setNum;
  if(!strcmp("Corr", Corr.c_str()))
    setCorrNum = setNum + 8;

  const Int_t nBins = retBinNumber(Tight);
  Float_t xArr[nBins+1];
  Float_t xArrCut[nBins+1];

  Int_t evtsPass = 0;

  if(nBins == 4) getBinArr(nBins, xArr, xArrCut, loose);
  else if(nBins == 8) getBinArr(nBins, xArr, xArrCut, tight);

  const std::string centString = getCentString(sType, centLow, centHi);

  TH1F* rImbAHist_p[6];
  TH1F* mean_rProjA_p[6][nBins];
  
  InitHist(rImbAHist_p);
  BookHist(rImbAHist_p, "r", algType[setNum], "Proj", "", Corr, centString, nBins, xArr, 0.50, niceNumCNC);

  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      mean_rProjA_p[iter][iter2] = new TH1F(Form("mean_rProjA_%d_%d_h", iter, iter2), Form("mean_rProjA_%d_%d_h", iter, iter2), 100, -2000, 2000);
    }
  }    

  for(Int_t jEntry = 0; jEntry < (Int_t)anaTree_p->GetEntries(); jEntry++){
    anaTree_p->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(isEventCut(setNum, sType, centLow, centHi, isHighPtTrk)) continue;

    evtsPass++;

    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      if(AlgJtAsymm_[setNum] < xArrCut[iter2+1]){

	for(Int_t iter = 0; iter < 6; iter++){
	  if(rAlgImbProjA_[setCorrNum][iter] == 0) continue;

	  Float_t weight = 1.0;
	  if(montecarlo){
	    weight = pthatWeight_;
	    if(isHI(sType)) weight *= centWeight_[setNum];
	  }

	  mean_rProjA_p[iter][iter2]->Fill(rAlgImbProjA_[setCorrNum][iter], weight);
	}

	break;	
      }
    }
  }

  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      if(mean_rProjA_p[iter][iter2]->GetEntries() != 0){
	rImbAHist_p[iter]->SetBinContent(iter2+1, mean_rProjA_p[iter][iter2]->GetMean());
	rImbAHist_p[iter]->SetBinError(iter2+1, mean_rProjA_p[iter][iter2]->GetMeanError());
      }
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
  
  std::ofstream txtFile;
  txtFile.open(Form("%s.txt", outName.c_str()), std::ofstream::out | std::ofstream::app);
  if(centLow == 0) txtFile << Form("ProjA %s,\n", algType[setNum].c_str());
  txtFile << Form("   %d-%d: %d.\n", centLow, centHi, evtsPass);
  txtFile.close();
  
  for(Int_t iter = 0; iter < 6; iter++){
    delete rImbAHist_p[iter];
    rImbAHist_p[iter] = 0;

    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      delete mean_rProjA_p[iter][iter2];
      mean_rProjA_p[iter][iter2] = 0;
    }
  }

}


void makeImbACNCHist(TTree* anaTree_p, const std::string outName, Int_t setNum, Int_t centLow, Int_t centHi, const std::string Corr = "", const std::string Tight = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);

  Int_t setCorrNum = setNum;
  if(!strcmp("Corr", Corr.c_str()))
    setCorrNum = setNum + 8;

  const Int_t nBins = retBinNumber(Tight);
  Float_t xArr[nBins+1];
  Float_t xArrCut[nBins+1];

  Int_t evtsPass = 0;

  if(nBins == 4) getBinArr(nBins, xArr, xArrCut, loose);
  else if(nBins == 8) getBinArr(nBins, xArr, xArrCut, tight);

  const std::string centString = getCentString(sType, centLow, centHi);

  TH1F* rImbACHist_p[6];
  TH1F* mean_rProjAC_p[6][nBins];
  TH1F* rImbANCHist_p[6];
  TH1F* mean_rProjANC_p[6][nBins];
  TH1F* rImbAC0Hist_p[6];
  TH1F* mean_rProjAC0_p[6][nBins];
  TH1F* rImbAC1Hist_p[6];
  TH1F* mean_rProjAC1_p[6][nBins];
  TH1F* rImbAC2Hist_p[6];
  TH1F* mean_rProjAC2_p[6][nBins];
  TH1F* rImbAC3Hist_p[6];
  TH1F* mean_rProjAC3_p[6][nBins];

  InitHist(rImbACHist_p);
  InitHist(rImbANCHist_p);
  InitHist(rImbAC0Hist_p);
  InitHist(rImbAC1Hist_p);
  InitHist(rImbAC2Hist_p);
  InitHist(rImbAC3Hist_p);

  BookHist(rImbACHist_p, "r", algType[setNum], "Proj", "C", Corr, centString, nBins, xArr, 0.50, niceNumCNC);
  BookHist(rImbANCHist_p, "r", algType[setNum], "Proj", "NC", Corr, centString, nBins, xArr, 0.50, niceNumCNC);
  BookHist(rImbAC0Hist_p, "r", algType[setNum], "Proj", "C0", Corr, centString, nBins, xArr, 0.50, niceNumCNC);
  BookHist(rImbAC1Hist_p, "r", algType[setNum], "Proj", "C1", Corr, centString, nBins, xArr, 0.50, niceNumCNC);
  BookHist(rImbAC2Hist_p, "r", algType[setNum], "Proj", "C2", Corr, centString, nBins, xArr, 0.50, niceNumCNC);
  BookHist(rImbAC3Hist_p, "r", algType[setNum], "Proj", "C3", Corr, centString, nBins, xArr, 0.50, niceNumCNC);

  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      mean_rProjAC_p[iter][iter2] = new TH1F(Form("mean_rProjAC_%d_%d_h", iter, iter2), Form("mean_rProjAC_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_rProjANC_p[iter][iter2] = new TH1F(Form("mean_rProjANC_%d_%d_h", iter, iter2), Form("mean_rProjANC_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_rProjAC0_p[iter][iter2] = new TH1F(Form("mean_rProjAC0_%d_%d_h", iter, iter2), Form("mean_rProjAC0_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_rProjAC1_p[iter][iter2] = new TH1F(Form("mean_rProjAC1_%d_%d_h", iter, iter2), Form("mean_rProjAC1_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_rProjAC2_p[iter][iter2] = new TH1F(Form("mean_rProjAC2_%d_%d_h", iter, iter2), Form("mean_rProjAC2_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_rProjAC3_p[iter][iter2] = new TH1F(Form("mean_rProjAC3_%d_%d_h", iter, iter2), Form("mean_rProjAC3_%d_%d_h", iter, iter2), 100, -2000, 2000);
    }
  }    

  for(Int_t jEntry = 0; jEntry < (Int_t)anaTree_p->GetEntries(); jEntry++){
    anaTree_p->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(isEventCut(setNum, sType, centLow, centHi, isHighPtTrk)) continue;

    evtsPass++;

    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      if(AlgJtAsymm_[setNum] < xArrCut[iter2+1]){

	for(Int_t iter = 0; iter < 6; iter++){
	  if(rAlgImbProjAC_[setCorrNum][iter] == 0) continue;

	  Float_t weight = 1.0;
	  if(montecarlo){
	    weight = pthatWeight_;
	    if(isHI(sType)) weight *= centWeight_[setNum];
	  }

	  mean_rProjAC_p[iter][iter2]->Fill(rAlgImbProjAC_[setCorrNum][iter], weight);
	  mean_rProjANC_p[iter][iter2]->Fill(rAlgImbProjANC_[setCorrNum][iter], weight);

	  if(TMath::Abs(AlgJtEta_[setNum][0]) < midRapCut && TMath::Abs(AlgJtEta_[setNum][1]) < midRapCut){
	    mean_rProjAC0_p[iter][iter2]->Fill(rAlgImbProjAC0_[setCorrNum][iter], weight);
	    mean_rProjAC1_p[iter][iter2]->Fill(rAlgImbProjAC1_[setCorrNum][iter], weight);
	    mean_rProjAC2_p[iter][iter2]->Fill(rAlgImbProjAC2_[setCorrNum][iter], weight);
	    mean_rProjAC3_p[iter][iter2]->Fill(rAlgImbProjAC3_[setCorrNum][iter], weight);
	  }
	}

	break;	
      }
    }
  }

  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){

      if(mean_rProjAC_p[iter][iter2]->GetEntries() != 0){
	rImbACHist_p[iter]->SetBinContent(iter2+1, mean_rProjAC_p[iter][iter2]->GetMean());
	rImbACHist_p[iter]->SetBinError(iter2+1, mean_rProjAC_p[iter][iter2]->GetMeanError());
      }

      if(mean_rProjANC_p[iter][iter2]->GetEntries() != 0){
	rImbANCHist_p[iter]->SetBinContent(iter2+1, mean_rProjANC_p[iter][iter2]->GetMean());
	rImbANCHist_p[iter]->SetBinError(iter2+1, mean_rProjANC_p[iter][iter2]->GetMeanError());
      }

      if(mean_rProjAC0_p[iter][iter2]->GetEntries() != 0){
	rImbAC0Hist_p[iter]->SetBinContent(iter2+1, mean_rProjAC0_p[iter][iter2]->GetMean());
	rImbAC0Hist_p[iter]->SetBinError(iter2+1, mean_rProjAC0_p[iter][iter2]->GetMeanError());
      }

      if(mean_rProjAC1_p[iter][iter2]->GetEntries() != 0){
	rImbAC1Hist_p[iter]->SetBinContent(iter2+1, mean_rProjAC1_p[iter][iter2]->GetMean());
	rImbAC1Hist_p[iter]->SetBinError(iter2+1, mean_rProjAC1_p[iter][iter2]->GetMeanError());
      }

      if(mean_rProjAC2_p[iter][iter2]->GetEntries() != 0){
	rImbAC2Hist_p[iter]->SetBinContent(iter2+1, mean_rProjAC2_p[iter][iter2]->GetMean());
	rImbAC2Hist_p[iter]->SetBinError(iter2+1, mean_rProjAC2_p[iter][iter2]->GetMeanError());
      }

      if(mean_rProjAC3_p[iter][iter2]->GetEntries() != 0){
	rImbAC3Hist_p[iter]->SetBinContent(iter2+1, mean_rProjAC3_p[iter][iter2]->GetMean());
	rImbAC3Hist_p[iter]->SetBinError(iter2+1, mean_rProjAC3_p[iter][iter2]->GetMeanError());
      }

    }
  }
  
  outFile_p = new TFile(Form("%s.root", outName.c_str()), "UPDATE");
  std::cout << outName << std::endl;

  for(Int_t iter = 0; iter < 6; iter++){
    rImbACHist_p[iter]->Write("", TObject::kOverwrite);
    rImbANCHist_p[iter]->Write("", TObject::kOverwrite);

    rImbAC0Hist_p[iter]->Write("", TObject::kOverwrite);
    rImbAC1Hist_p[iter]->Write("", TObject::kOverwrite);
    rImbAC2Hist_p[iter]->Write("", TObject::kOverwrite);
    rImbAC3Hist_p[iter]->Write("", TObject::kOverwrite);
  }

  outFile_p->Close();
  delete outFile_p;
  outFile_p = 0;

  std::ofstream txtFile;
  txtFile.open(Form("%s.txt", outName.c_str()), std::ofstream::out | std::ofstream::app);
  if(centLow == 0) txtFile << Form("ProjACNC %s,\n", algType[setNum].c_str());
  txtFile << Form("   %d-%d: %d.\n", centLow, centHi, evtsPass);
  txtFile.close();

  for(Int_t iter = 0; iter < 6; iter++){
    delete rImbACHist_p[iter];
    rImbACHist_p[iter] = 0;
    delete rImbANCHist_p[iter];
    rImbANCHist_p[iter] = 0;

    delete rImbAC0Hist_p[iter];
    rImbAC0Hist_p[iter] = 0;
    delete rImbAC1Hist_p[iter];
    rImbAC1Hist_p[iter] = 0;
    delete rImbAC2Hist_p[iter];
    rImbAC2Hist_p[iter] = 0;
    delete rImbAC3Hist_p[iter];
    rImbAC3Hist_p[iter] = 0;


    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      delete mean_rProjAC_p[iter][iter2];
      mean_rProjAC_p[iter][iter2] = 0;
      delete mean_rProjANC_p[iter][iter2];
      mean_rProjANC_p[iter][iter2] = 0;

      delete mean_rProjAC0_p[iter][iter2];
      mean_rProjAC0_p[iter][iter2] = 0;
      delete mean_rProjAC1_p[iter][iter2];
      mean_rProjAC1_p[iter][iter2] = 0;
      delete mean_rProjAC2_p[iter][iter2];
      mean_rProjAC2_p[iter][iter2] = 0;
      delete mean_rProjAC3_p[iter][iter2];
      mean_rProjAC3_p[iter][iter2] = 0;
    }
  }

}



void makeImbARHist(TTree* anaTree_p, const std::string outName, const std::string gR, Int_t setNum, Int_t centLow, Int_t centHi, const std::string Corr = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);                                                    

  Int_t setCorrNum = setNum;
  if(!strcmp("Corr", Corr.c_str())) setCorrNum = setNum + 8;

  Int_t evtsPass = 0;
  Int_t evtsPassU = 0;
  Int_t evtsPassD = 0;

  const Int_t nBins = 10;
  Float_t rBins[11] = {0.0001, 0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 1.9999};
  const std::string centString = getCentString(sType, centLow, centHi);

  const Int_t nBins2 = 9;
  Float_t rBins2[10] = {0.0001, 0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.7999};

  TH1F* imbARHist_p[6];
  TH1F* imbARDHist_p[6];
  TH1F* imbARUHist_p[6];

  TH1F* imbAR2Hist_p[6];
  TH1F* imbAR2DHist_p[6];
  TH1F* imbAR2UHist_p[6];

  TH1F* imbARCutHist_p[6];
  TH1F* imbARCutDHist_p[6];
  TH1F* imbARCutUHist_p[6];
  TH1F* imbARCutEtaHist_p[6];
  TH1F* imbARCutEtaDHist_p[6];
  TH1F* imbARCutEtaUHist_p[6];
  TH1F* imbARCutPhiHist_p[6];
  TH1F* imbARCutPhiDHist_p[6];
  TH1F* imbARCutPhiUHist_p[6];
  TH1F* imbAEtaHist_p[6];
  TH1F* imbAEtaDHist_p[6];
  TH1F* imbAEtaUHist_p[6];
  TH1F* imbAEtaCutHist_p[6];
  TH1F* imbAEtaCutDHist_p[6];
  TH1F* imbAEtaCutUHist_p[6];
  TH1F* imbAPhiHist_p[6];
  TH1F* imbAPhiDHist_p[6];
  TH1F* imbAPhiUHist_p[6];
  TH1F* imbAPhiCutHist_p[6];
  TH1F* imbAPhiCutDHist_p[6];
  TH1F* imbAPhiCutUHist_p[6];

  TH1F* imbARFORHist_p[6];
  TH1F* imbARFORUHist_p[6];
  TH1F* imbARFORDHist_p[6];
  TH1F* imbARFORMIDHist_p[6];
  TH1F* imbARFORMIDUHist_p[6];
  TH1F* imbARFORMIDDHist_p[6];
  TH1F* imbARFORFORHist_p[6];
  TH1F* imbARFORFORUHist_p[6];
  TH1F* imbARFORFORDHist_p[6];
  TH1F* imbARMIDHist_p[6];
  TH1F* imbARMIDUHist_p[6];
  TH1F* imbARMIDDHist_p[6];

  TH1F* multARHist_p[6];
  TH1F* multARDHist_p[6];
  TH1F* multARUHist_p[6];
  TH1F* multARCutHist_p[6];
  TH1F* multARCutDHist_p[6];
  TH1F* multARCutUHist_p[6];
  TH1F* multARCutEtaHist_p[6];
  TH1F* multARCutEtaDHist_p[6];
  TH1F* multARCutEtaUHist_p[6];
  TH1F* multARCutPhiHist_p[6];
  TH1F* multARCutPhiDHist_p[6];
  TH1F* multARCutPhiUHist_p[6];
  TH1F* multAEtaHist_p[6];
  TH1F* multAEtaDHist_p[6];
  TH1F* multAEtaUHist_p[6];
  TH1F* multAEtaCutHist_p[6];
  TH1F* multAEtaCutDHist_p[6];
  TH1F* multAEtaCutUHist_p[6];
  TH1F* multAPhiHist_p[6];
  TH1F* multAPhiDHist_p[6];
  TH1F* multAPhiUHist_p[6];
  TH1F* multAPhiCutHist_p[6];
  TH1F* multAPhiCutDHist_p[6];
  TH1F* multAPhiCutUHist_p[6];

  TH1F* mean_projAR_p[6][nBins];
  TH1F* mean_projARD_p[6][nBins];
  TH1F* mean_projARU_p[6][nBins];

  TH1F* mean_projAR2_p[6][nBins2];
  TH1F* mean_projAR2D_p[6][nBins2];
  TH1F* mean_projAR2U_p[6][nBins2];

  TH1F* mean_projARCut_p[6][nBins];
  TH1F* mean_projARCutD_p[6][nBins];
  TH1F* mean_projARCutU_p[6][nBins];
  TH1F* mean_projARCutEta_p[6][nBins];
  TH1F* mean_projARCutEtaD_p[6][nBins];
  TH1F* mean_projARCutEtaU_p[6][nBins];
  TH1F* mean_projARCutPhi_p[6][nBins];
  TH1F* mean_projARCutPhiD_p[6][nBins];
  TH1F* mean_projARCutPhiU_p[6][nBins];
  TH1F* mean_projAEta_p[6][nBins];
  TH1F* mean_projAEtaD_p[6][nBins];
  TH1F* mean_projAEtaU_p[6][nBins];
  TH1F* mean_projAEtaCut_p[6][nBins];
  TH1F* mean_projAEtaCutD_p[6][nBins];
  TH1F* mean_projAEtaCutU_p[6][nBins];
  TH1F* mean_projAPhi_p[6][nBins];
  TH1F* mean_projAPhiD_p[6][nBins];
  TH1F* mean_projAPhiU_p[6][nBins];
  TH1F* mean_projAPhiCut_p[6][nBins];
  TH1F* mean_projAPhiCutD_p[6][nBins];
  TH1F* mean_projAPhiCutU_p[6][nBins];

  TH1F* mean_projARFOR_p[6][nBins];
  TH1F* mean_projARFORD_p[6][nBins];
  TH1F* mean_projARFORU_p[6][nBins];
  TH1F* mean_projARFORMID_p[6][nBins];
  TH1F* mean_projARFORMIDD_p[6][nBins];
  TH1F* mean_projARFORMIDU_p[6][nBins];
  TH1F* mean_projARFORFOR_p[6][nBins];
  TH1F* mean_projARFORFORD_p[6][nBins];
  TH1F* mean_projARFORFORU_p[6][nBins];
  TH1F* mean_projARMID_p[6][nBins];
  TH1F* mean_projARMIDD_p[6][nBins];
  TH1F* mean_projARMIDU_p[6][nBins];

  TH1F* mean_multAR_p[6][nBins];
  TH1F* mean_multARD_p[6][nBins];
  TH1F* mean_multARU_p[6][nBins];
  TH1F* mean_multARCut_p[6][nBins];
  TH1F* mean_multARCutD_p[6][nBins];
  TH1F* mean_multARCutU_p[6][nBins];
  TH1F* mean_multARCutEta_p[6][nBins];
  TH1F* mean_multARCutEtaD_p[6][nBins];
  TH1F* mean_multARCutEtaU_p[6][nBins];
  TH1F* mean_multARCutPhi_p[6][nBins];
  TH1F* mean_multARCutPhiD_p[6][nBins];
  TH1F* mean_multARCutPhiU_p[6][nBins];
  TH1F* mean_multAEta_p[6][nBins];
  TH1F* mean_multAEtaD_p[6][nBins];
  TH1F* mean_multAEtaU_p[6][nBins];
  TH1F* mean_multAEtaCut_p[6][nBins];
  TH1F* mean_multAEtaCutD_p[6][nBins];
  TH1F* mean_multAEtaCutU_p[6][nBins];
  TH1F* mean_multAPhi_p[6][nBins];
  TH1F* mean_multAPhiD_p[6][nBins];
  TH1F* mean_multAPhiU_p[6][nBins];
  TH1F* mean_multAPhiCut_p[6][nBins];
  TH1F* mean_multAPhiCutD_p[6][nBins];
  TH1F* mean_multAPhiCutU_p[6][nBins];

  InitHist(imbARHist_p);
  InitHist(imbARDHist_p);
  InitHist(imbARUHist_p);

  InitHist(imbAR2Hist_p);
  InitHist(imbAR2DHist_p);
  InitHist(imbAR2UHist_p);

  InitHist(imbARCutHist_p);
  InitHist(imbARCutDHist_p);
  InitHist(imbARCutUHist_p);
  InitHist(imbARCutEtaHist_p);
  InitHist(imbARCutEtaDHist_p);
  InitHist(imbARCutEtaUHist_p);
  InitHist(imbARCutPhiHist_p);
  InitHist(imbARCutPhiDHist_p);
  InitHist(imbARCutPhiUHist_p);
  InitHist(imbAEtaHist_p);
  InitHist(imbAEtaDHist_p);
  InitHist(imbAEtaUHist_p);
  InitHist(imbAEtaCutHist_p);
  InitHist(imbAEtaCutDHist_p);
  InitHist(imbAEtaCutUHist_p);
  InitHist(imbAPhiHist_p);
  InitHist(imbAPhiDHist_p);
  InitHist(imbAPhiUHist_p);
  InitHist(imbAPhiCutHist_p);
  InitHist(imbAPhiCutDHist_p);
  InitHist(imbAPhiCutUHist_p);

  InitHist(imbARFORHist_p);
  InitHist(imbARFORDHist_p);
  InitHist(imbARFORUHist_p);
  InitHist(imbARFORMIDHist_p);
  InitHist(imbARFORMIDDHist_p);
  InitHist(imbARFORMIDUHist_p);
  InitHist(imbARFORFORHist_p);
  InitHist(imbARFORFORDHist_p);
  InitHist(imbARFORFORUHist_p);
  InitHist(imbARMIDHist_p);
  InitHist(imbARMIDDHist_p);
  InitHist(imbARMIDUHist_p);

  InitHist(multARHist_p);
  InitHist(multARDHist_p);
  InitHist(multARUHist_p);
  InitHist(multARCutHist_p);
  InitHist(multARCutDHist_p);
  InitHist(multARCutUHist_p);
  InitHist(multARCutEtaHist_p);
  InitHist(multARCutEtaDHist_p);
  InitHist(multARCutEtaUHist_p);
  InitHist(multARCutPhiHist_p);
  InitHist(multARCutPhiDHist_p);
  InitHist(multARCutPhiUHist_p);
  InitHist(multAEtaHist_p);
  InitHist(multAEtaDHist_p);
  InitHist(multAEtaUHist_p);
  InitHist(multAEtaCutHist_p);
  InitHist(multAEtaCutDHist_p);
  InitHist(multAEtaCutUHist_p);
  InitHist(multAPhiHist_p);
  InitHist(multAPhiDHist_p);
  InitHist(multAPhiUHist_p);
  InitHist(multAPhiCutHist_p);
  InitHist(multAPhiCutDHist_p);
  InitHist(multAPhiCutUHist_p);

  BookHist(imbARHist_p, gR, algType[setNum], "Proj", "R", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARDHist_p, gR, algType[setNum], "Proj", "RD", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARUHist_p, gR, algType[setNum], "Proj", "RU", Corr, centString, nBins, rBins, 2.00, niceNumR);

  BookHist(imbAR2Hist_p, gR, algType[setNum], "Proj", "R2", Corr, centString, nBins2, rBins2, 2.00, niceNumR);
  BookHist(imbAR2DHist_p, gR, algType[setNum], "Proj", "R2D", Corr, centString, nBins2, rBins2, 2.00, niceNumR);
  BookHist(imbAR2UHist_p, gR, algType[setNum], "Proj", "R2U", Corr, centString, nBins2, rBins2, 2.00, niceNumR);

  BookHist(imbARCutHist_p, gR, algType[setNum], "Proj", "RCut", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARCutDHist_p, gR, algType[setNum], "Proj", "RCutD", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARCutUHist_p, gR, algType[setNum], "Proj", "RCutU", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARCutEtaHist_p, gR, algType[setNum], "Proj", "RCutEta", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARCutEtaDHist_p, gR, algType[setNum], "Proj", "RCutEtaD", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARCutEtaUHist_p, gR, algType[setNum], "Proj", "RCutEtaU", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARCutPhiHist_p, gR, algType[setNum], "Proj", "RCutPhi", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARCutPhiDHist_p, gR, algType[setNum], "Proj", "RCutPhiD", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARCutPhiUHist_p, gR, algType[setNum], "Proj", "RCutPhiU", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbAEtaHist_p, gR, algType[setNum], "Proj", "Eta", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbAEtaDHist_p, gR, algType[setNum], "Proj", "EtaD", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbAEtaUHist_p, gR, algType[setNum], "Proj", "EtaU", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbAEtaCutHist_p, gR, algType[setNum], "Proj", "EtaCut", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbAEtaCutDHist_p, gR, algType[setNum], "Proj", "EtaCutD", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbAEtaCutUHist_p, gR, algType[setNum], "Proj", "EtaCutU", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbAPhiHist_p, gR, algType[setNum], "Proj", "Phi", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbAPhiDHist_p, gR, algType[setNum], "Proj", "PhiD", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbAPhiUHist_p, gR, algType[setNum], "Proj", "PhiU", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbAPhiCutHist_p, gR, algType[setNum], "Proj", "PhiCut", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbAPhiCutDHist_p, gR, algType[setNum], "Proj", "PhiCutD", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbAPhiCutUHist_p, gR, algType[setNum], "Proj", "PhiCutU", Corr, centString, nBins, rBins, 2.00, niceNumR);

  BookHist(imbARFORHist_p, gR, algType[setNum], "Proj", "RFOR", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARFORDHist_p, gR, algType[setNum], "Proj", "RFORD", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARFORUHist_p, gR, algType[setNum], "Proj", "RFORU", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARFORMIDHist_p, gR, algType[setNum], "Proj", "RFORMID", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARFORMIDDHist_p, gR, algType[setNum], "Proj", "RFORMIDD", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARFORMIDUHist_p, gR, algType[setNum], "Proj", "RFORMIDU", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARFORFORHist_p, gR, algType[setNum], "Proj", "RFORFOR", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARFORFORDHist_p, gR, algType[setNum], "Proj", "RFORFORD", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARFORFORUHist_p, gR, algType[setNum], "Proj", "RFORFORU", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARMIDHist_p, gR, algType[setNum], "Proj", "RMID", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARMIDDHist_p, gR, algType[setNum], "Proj", "RMIDD", Corr, centString, nBins, rBins, 2.00, niceNumR);
  BookHist(imbARMIDUHist_p, gR, algType[setNum], "Proj", "RMIDU", Corr, centString, nBins, rBins, 2.00, niceNumR);

  BookHist(multARHist_p, gR, algType[setNum], "Mult", "R", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multARDHist_p, gR, algType[setNum], "Mult", "RD", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multARUHist_p, gR, algType[setNum], "Mult", "RU", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multARCutHist_p, gR, algType[setNum], "Mult", "RCut", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multARCutDHist_p, gR, algType[setNum], "Mult", "RCutD", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multARCutUHist_p, gR, algType[setNum], "Mult", "RCutU", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multARCutEtaHist_p, gR, algType[setNum], "Mult", "RCutEta", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multARCutEtaDHist_p, gR, algType[setNum], "Mult", "RCutEtaD", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multARCutEtaUHist_p, gR, algType[setNum], "Mult", "RCutEtaU", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multARCutPhiHist_p, gR, algType[setNum], "Mult", "RCutPhi", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multARCutPhiDHist_p, gR, algType[setNum], "Mult", "RCutPhiD", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multARCutPhiUHist_p, gR, algType[setNum], "Mult", "RCutPhiU", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multAEtaHist_p, gR, algType[setNum], "Mult", "Eta", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multAEtaDHist_p, gR, algType[setNum], "Mult", "EtaD", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multAEtaUHist_p, gR, algType[setNum], "Mult", "EtaU", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multAEtaCutHist_p, gR, algType[setNum], "Mult", "EtaCut", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multAEtaCutDHist_p, gR, algType[setNum], "Mult", "EtaCutD", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multAEtaCutUHist_p, gR, algType[setNum], "Mult", "EtaCutU", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multAPhiHist_p, gR, algType[setNum], "Mult", "Phi", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multAPhiDHist_p, gR, algType[setNum], "Mult", "PhiD", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multAPhiUHist_p, gR, algType[setNum], "Mult", "PhiU", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multAPhiCutHist_p, gR, algType[setNum], "Mult", "PhiCut", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multAPhiCutDHist_p, gR, algType[setNum], "Mult", "PhiCutD", Corr, centString, nBins, rBins, 2.00, niceNumMult);
  BookHist(multAPhiCutUHist_p, gR, algType[setNum], "Mult", "PhiCutU", Corr, centString, nBins, rBins, 2.00, niceNumMult);


  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      mean_projAR_p[iter][iter2] = new TH1F(Form("mean_projAR_%d_%d_h", iter, iter2), Form("mean_projAR_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARD_p[iter][iter2] = new TH1F(Form("mean_projARD_%d_%d_h", iter, iter2), Form("mean_projARD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARU_p[iter][iter2] = new TH1F(Form("mean_projARU_%d_%d_h", iter, iter2), Form("mean_projARU_%d_%d_h", iter, iter2), 100, -2000, 2000);

      if(iter2 < nBins2){
	mean_projAR2_p[iter][iter2] = new TH1F(Form("mean_projAR2_%d_%d_h", iter, iter2), Form("mean_projAR2_%d_%d_h", iter, iter2), 100, -2000, 2000);
	mean_projAR2D_p[iter][iter2] = new TH1F(Form("mean_projAR2D_%d_%d_h", iter, iter2), Form("mean_projAR2D_%d_%d_h", iter, iter2), 100, -2000, 2000);
	mean_projAR2U_p[iter][iter2] = new TH1F(Form("mean_projAR2U_%d_%d_h", iter, iter2), Form("mean_projAR2U_%d_%d_h", iter, iter2), 100, -2000, 2000);
      }

      mean_projARCut_p[iter][iter2] = new TH1F(Form("mean_projARCut_%d_%d_h", iter, iter2), Form("mean_projARCut_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARCutD_p[iter][iter2] = new TH1F(Form("mean_projARCutD_%d_%d_h", iter, iter2), Form("mean_projARCutD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARCutU_p[iter][iter2] = new TH1F(Form("mean_projARCutU_%d_%d_h", iter, iter2), Form("mean_projARCutU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARCutEta_p[iter][iter2] = new TH1F(Form("mean_projARCutEta_%d_%d_h", iter, iter2), Form("mean_projARCutEta_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARCutEtaD_p[iter][iter2] = new TH1F(Form("mean_projARCutEtaD_%d_%d_h", iter, iter2), Form("mean_projARCutEtaD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARCutEtaU_p[iter][iter2] = new TH1F(Form("mean_projARCutEtaU_%d_%d_h", iter, iter2), Form("mean_projARCutEtaU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARCutPhi_p[iter][iter2] = new TH1F(Form("mean_projARCutPhi_%d_%d_h", iter, iter2), Form("mean_projARCutPhi_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARCutPhiD_p[iter][iter2] = new TH1F(Form("mean_projARCutPhiD_%d_%d_h", iter, iter2), Form("mean_projARCutPhiD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARCutPhiU_p[iter][iter2] = new TH1F(Form("mean_projARCutPhiU_%d_%d_h", iter, iter2), Form("mean_projARCutPhiU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projAEta_p[iter][iter2] = new TH1F(Form("mean_projAEta_%d_%d_h", iter, iter2), Form("mean_projAEta_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projAEtaD_p[iter][iter2] = new TH1F(Form("mean_projAEtaD_%d_%d_h", iter, iter2), Form("mean_projAEtaD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projAEtaU_p[iter][iter2] = new TH1F(Form("mean_projAEtaU_%d_%d_h", iter, iter2), Form("mean_projAEtaU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projAEtaCut_p[iter][iter2] = new TH1F(Form("mean_projAEtaCut_%d_%d_h", iter, iter2), Form("mean_projAEtaCut_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projAEtaCutD_p[iter][iter2] = new TH1F(Form("mean_projAEtaCutD_%d_%d_h", iter, iter2), Form("mean_projAEtaCutD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projAEtaCutU_p[iter][iter2] = new TH1F(Form("mean_projAEtaCutU_%d_%d_h", iter, iter2), Form("mean_projAEtaCutU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projAPhi_p[iter][iter2] = new TH1F(Form("mean_projAPhi_%d_%d_h", iter, iter2), Form("mean_projAPhi_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projAPhiD_p[iter][iter2] = new TH1F(Form("mean_projAPhiD_%d_%d_h", iter, iter2), Form("mean_projAPhiD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projAPhiU_p[iter][iter2] = new TH1F(Form("mean_projAPhiU_%d_%d_h", iter, iter2), Form("mean_projAPhiU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projAPhiCut_p[iter][iter2] = new TH1F(Form("mean_projAPhiCut_%d_%d_h", iter, iter2), Form("mean_projAPhiCut_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projAPhiCutD_p[iter][iter2] = new TH1F(Form("mean_projAPhiCutD_%d_%d_h", iter, iter2), Form("mean_projAPhiCutD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projAPhiCutU_p[iter][iter2] = new TH1F(Form("mean_projAPhiCutU_%d_%d_h", iter, iter2), Form("mean_projAPhiCutU_%d_%d_h", iter, iter2), 100, -2000, 2000);


      mean_projARFOR_p[iter][iter2] = new TH1F(Form("mean_projARFOR_%d_%d_h", iter, iter2), Form("mean_projARFOR_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARFORD_p[iter][iter2] = new TH1F(Form("mean_projARFORD_%d_%d_h", iter, iter2), Form("mean_projARFORD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARFORU_p[iter][iter2] = new TH1F(Form("mean_projARFORU_%d_%d_h", iter, iter2), Form("mean_projARFORU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARFORMID_p[iter][iter2] = new TH1F(Form("mean_projARFORMID_%d_%d_h", iter, iter2), Form("mean_projARFORMID_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARFORMIDD_p[iter][iter2] = new TH1F(Form("mean_projARFORMIDD_%d_%d_h", iter, iter2), Form("mean_projARFORMIDD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARFORMIDU_p[iter][iter2] = new TH1F(Form("mean_projARFORMIDU_%d_%d_h", iter, iter2), Form("mean_projARFORMIDU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARFORFOR_p[iter][iter2] = new TH1F(Form("mean_projARFORFOR_%d_%d_h", iter, iter2), Form("mean_projARFORFOR_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARFORFORD_p[iter][iter2] = new TH1F(Form("mean_projARFORFORD_%d_%d_h", iter, iter2), Form("mean_projARFORFORD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARFORFORU_p[iter][iter2] = new TH1F(Form("mean_projARFORFORU_%d_%d_h", iter, iter2), Form("mean_projARFORFORU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARMID_p[iter][iter2] = new TH1F(Form("mean_projARMID_%d_%d_h", iter, iter2), Form("mean_projARMID_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARMIDD_p[iter][iter2] = new TH1F(Form("mean_projARMIDD_%d_%d_h", iter, iter2), Form("mean_projARMIDD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_projARMIDU_p[iter][iter2] = new TH1F(Form("mean_projARMIDU_%d_%d_h", iter, iter2), Form("mean_projARMIDU_%d_%d_h", iter, iter2), 100, -2000, 2000);


      mean_multAR_p[iter][iter2] = new TH1F(Form("mean_multAR_%d_%d_h", iter, iter2), Form("mean_multAR_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multARD_p[iter][iter2] = new TH1F(Form("mean_multARD_%d_%d_h", iter, iter2), Form("mean_multARD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multARU_p[iter][iter2] = new TH1F(Form("mean_multARU_%d_%d_h", iter, iter2), Form("mean_multARU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multARCut_p[iter][iter2] = new TH1F(Form("mean_multARCut_%d_%d_h", iter, iter2), Form("mean_multARCut_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multARCutD_p[iter][iter2] = new TH1F(Form("mean_multARCutD_%d_%d_h", iter, iter2), Form("mean_multARCutD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multARCutU_p[iter][iter2] = new TH1F(Form("mean_multARCutU_%d_%d_h", iter, iter2), Form("mean_multARCutU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multARCutEta_p[iter][iter2] = new TH1F(Form("mean_multARCutEta_%d_%d_h", iter, iter2), Form("mean_multARCutEta_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multARCutEtaD_p[iter][iter2] = new TH1F(Form("mean_multARCutEtaD_%d_%d_h", iter, iter2), Form("mean_multARCutEtaD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multARCutEtaU_p[iter][iter2] = new TH1F(Form("mean_multARCutEtaU_%d_%d_h", iter, iter2), Form("mean_multARCutEtaU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multARCutPhi_p[iter][iter2] = new TH1F(Form("mean_multARCutPhi_%d_%d_h", iter, iter2), Form("mean_multARCutPhi_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multARCutPhiD_p[iter][iter2] = new TH1F(Form("mean_multARCutPhiD_%d_%d_h", iter, iter2), Form("mean_multARCutPhiD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multARCutPhiU_p[iter][iter2] = new TH1F(Form("mean_multARCutPhiU_%d_%d_h", iter, iter2), Form("mean_multARCutPhiU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multAEta_p[iter][iter2] = new TH1F(Form("mean_multAEta_%d_%d_h", iter, iter2), Form("mean_multAEta_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multAEtaD_p[iter][iter2] = new TH1F(Form("mean_multAEtaD_%d_%d_h", iter, iter2), Form("mean_multAEtaD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multAEtaU_p[iter][iter2] = new TH1F(Form("mean_multAEtaU_%d_%d_h", iter, iter2), Form("mean_multAEtaU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multAEtaCut_p[iter][iter2] = new TH1F(Form("mean_multAEtaCut_%d_%d_h", iter, iter2), Form("mean_multAEtaCut_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multAEtaCutD_p[iter][iter2] = new TH1F(Form("mean_multAEtaCutD_%d_%d_h", iter, iter2), Form("mean_multAEtaCutD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multAEtaCutU_p[iter][iter2] = new TH1F(Form("mean_multAEtaCutU_%d_%d_h", iter, iter2), Form("mean_multAEtaCutU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multAPhi_p[iter][iter2] = new TH1F(Form("mean_multAPhi_%d_%d_h", iter, iter2), Form("mean_multAPhi_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multAPhiD_p[iter][iter2] = new TH1F(Form("mean_multAPhiD_%d_%d_h", iter, iter2), Form("mean_multAPhiD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multAPhiU_p[iter][iter2] = new TH1F(Form("mean_multAPhiU_%d_%d_h", iter, iter2), Form("mean_multAPhiU_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multAPhiCut_p[iter][iter2] = new TH1F(Form("mean_multAPhiCut_%d_%d_h", iter, iter2), Form("mean_multAPhiCut_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multAPhiCutD_p[iter][iter2] = new TH1F(Form("mean_multAPhiCutD_%d_%d_h", iter, iter2), Form("mean_multAPhiCutD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_multAPhiCutU_p[iter][iter2] = new TH1F(Form("mean_multAPhiCutU_%d_%d_h", iter, iter2), Form("mean_multAPhiCutU_%d_%d_h", iter, iter2), 100, -2000, 2000);
    }
  }

  for(Int_t jEntry = 0; jEntry < (Int_t)anaTree_p->GetEntries(); jEntry++){
    anaTree_p->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(isEventCut(setNum, sType, centLow, centHi, isHighPtTrk)) continue;

    if(TMath::Abs(AlgJtEta_[setNum][0]) > midRapCut2 || TMath::Abs(AlgJtEta_[setNum][1]) > midRapCut2) continue;

    Float_t weight = 1.0;
    if(montecarlo){
      weight = pthatWeight_;
      if(isHI(sType)) weight *= centWeight_[setNum];
    }

    for(Int_t iter = 0; iter < 6; iter++){
      for(Int_t iter2 = 0; iter2 < nBins2-1; iter2++){
	if(strcmp(gR.c_str(), "r") == 0){
          mean_projAR2_p[iter][iter2]->Fill(rAlgImbProjAR_[setCorrNum][iter][iter2], weight);

          if(AlgJtAsymm_[setNum] < 0.22) mean_projAR2D_p[iter][iter2]->Fill(rAlgImbProjAR_[setCorrNum][iter][iter2], weight);
	  else mean_projAR2U_p[iter][iter2]->Fill(rAlgImbProjAR_[setCorrNum][iter][iter2], weight);
	}
        else if(strcmp(gR.c_str(), "g") == 0){
          mean_projAR2_p[iter][iter2]->Fill(gAlgImbProjAR_[setNum][iter][iter2], weight);

          if(AlgJtAsymm_[setNum] < 0.22) mean_projAR2D_p[iter][iter2]->Fill(gAlgImbProjAR_[setNum][iter][iter2], weight);
	  else mean_projAR2U_p[iter][iter2]->Fill(gAlgImbProjAR_[setNum][iter][iter2], weight);
	}
      }


      if(strcmp(gR.c_str(), "r") == 0){
	mean_projAR2_p[iter][nBins2-1]->Fill(rAlgImbProjAR_[setCorrNum][iter][nBins2-1] + rAlgImbProjAR_[setCorrNum][iter][nBins2], weight);
	
	if(AlgJtAsymm_[setNum] < 0.22) mean_projAR2D_p[iter][nBins2-1]->Fill(rAlgImbProjAR_[setCorrNum][iter][nBins2-1] + rAlgImbProjAR_[setCorrNum][iter][nBins2], weight);
	else mean_projAR2U_p[iter][nBins2-1]->Fill(rAlgImbProjAR_[setCorrNum][iter][nBins2-1] + rAlgImbProjAR_[setCorrNum][iter][nBins2], weight);
      }
      else if(strcmp(gR.c_str(), "g") == 0){
	mean_projAR2_p[iter][nBins2-1]->Fill(gAlgImbProjAR_[setNum][iter][nBins2-1] + gAlgImbProjAR_[setNum][iter][nBins2], weight);
	
	if(AlgJtAsymm_[setNum] < 0.22) mean_projAR2D_p[iter][nBins2-1]->Fill(gAlgImbProjAR_[setNum][iter][nBins2-1] + gAlgImbProjAR_[setNum][iter][nBins2], weight);
	else mean_projAR2U_p[iter][nBins2-1]->Fill(gAlgImbProjAR_[setNum][iter][nBins2-1] + gAlgImbProjAR_[setNum][iter][nBins2], weight);
      }
    }


    if(TMath::Abs(AlgJtEta_[setNum][0]) > midRapCut || TMath::Abs(AlgJtEta_[setNum][1]) > midRapCut) continue;

    evtsPass++;

    if(AlgJtAsymm_[setNum] < 0.22) evtsPassD++;
    else evtsPassU++;


    for(Int_t iter = 0; iter < 6; iter++){
      for(Int_t iter2 = 0; iter2 < nBins; iter2++){
	if(strcmp(gR.c_str(), "r") == 0){
	  mean_projAR_p[iter][iter2]->Fill(rAlgImbProjAR_[setCorrNum][iter][iter2], weight);
	  mean_projAEta_p[iter][iter2]->Fill(rAlgImbProjAEta_[setCorrNum][iter][iter2], weight);
	  mean_projAPhi_p[iter][iter2]->Fill(rAlgImbProjAPhi_[setCorrNum][iter][iter2], weight);
	  mean_projARCut_p[iter][iter2]->Fill(rAlgImbProjAR_Cut_[setCorrNum][iter][iter2], weight);
	  mean_projARCutEta_p[iter][iter2]->Fill(rAlgImbProjAR_CutEta_[setCorrNum][iter][iter2], weight);
	  mean_projARCutPhi_p[iter][iter2]->Fill(rAlgImbProjAR_CutPhi_[setCorrNum][iter][iter2], weight);
	  mean_projAEtaCut_p[iter][iter2]->Fill(rAlgImbProjAEta_Cut_[setCorrNum][iter][iter2], weight);
	  mean_projAPhiCut_p[iter][iter2]->Fill(rAlgImbProjAPhi_Cut_[setCorrNum][iter][iter2], weight);

	  mean_multAR_p[iter][iter2]->Fill(rAlgMultAR_[setCorrNum][iter][iter2], weight);
	  mean_multAEta_p[iter][iter2]->Fill(rAlgMultAEta_[setCorrNum][iter][iter2], weight);
	  mean_multAPhi_p[iter][iter2]->Fill(rAlgMultAPhi_[setCorrNum][iter][iter2], weight);
	  mean_multARCut_p[iter][iter2]->Fill(rAlgMultAR_Cut_[setCorrNum][iter][iter2], weight);
	  mean_multARCutEta_p[iter][iter2]->Fill(rAlgMultAR_CutEta_[setCorrNum][iter][iter2], weight);
	  mean_multARCutPhi_p[iter][iter2]->Fill(rAlgMultAR_CutPhi_[setCorrNum][iter][iter2], weight);
	  mean_multAEtaCut_p[iter][iter2]->Fill(rAlgMultAEta_Cut_[setCorrNum][iter][iter2], weight);
	  mean_multAPhiCut_p[iter][iter2]->Fill(rAlgMultAPhi_Cut_[setCorrNum][iter][iter2], weight);

	  if(AlgJtAsymm_[setNum] < 0.22){
	    mean_projARD_p[iter][iter2]->Fill(rAlgImbProjAR_[setCorrNum][iter][iter2], weight);
	    mean_projAEtaD_p[iter][iter2]->Fill(rAlgImbProjAEta_[setCorrNum][iter][iter2], weight);
	    mean_projAPhiD_p[iter][iter2]->Fill(rAlgImbProjAPhi_[setCorrNum][iter][iter2], weight);
	    mean_projARCutD_p[iter][iter2]->Fill(rAlgImbProjAR_Cut_[setCorrNum][iter][iter2], weight);
	    mean_projARCutEtaD_p[iter][iter2]->Fill(rAlgImbProjAR_CutEta_[setCorrNum][iter][iter2], weight);
	    mean_projARCutPhiD_p[iter][iter2]->Fill(rAlgImbProjAR_CutPhi_[setCorrNum][iter][iter2], weight);
	    mean_projAEtaCutD_p[iter][iter2]->Fill(rAlgImbProjAEta_Cut_[setCorrNum][iter][iter2], weight);
	    mean_projAPhiCutD_p[iter][iter2]->Fill(rAlgImbProjAPhi_Cut_[setCorrNum][iter][iter2], weight);

	    mean_multARD_p[iter][iter2]->Fill(rAlgMultAR_[setCorrNum][iter][iter2], weight);
	    mean_multAEtaD_p[iter][iter2]->Fill(rAlgMultAEta_[setCorrNum][iter][iter2], weight);
	    mean_multAPhiD_p[iter][iter2]->Fill(rAlgMultAPhi_[setCorrNum][iter][iter2], weight);
	    mean_multARCutD_p[iter][iter2]->Fill(rAlgMultAR_Cut_[setCorrNum][iter][iter2], weight);
	    mean_multARCutEtaD_p[iter][iter2]->Fill(rAlgMultAR_CutEta_[setCorrNum][iter][iter2], weight);
	    mean_multARCutPhiD_p[iter][iter2]->Fill(rAlgMultAR_CutPhi_[setCorrNum][iter][iter2], weight);
	    mean_multAEtaCutD_p[iter][iter2]->Fill(rAlgMultAEta_Cut_[setCorrNum][iter][iter2], weight);
	    mean_multAPhiCutD_p[iter][iter2]->Fill(rAlgMultAPhi_Cut_[setCorrNum][iter][iter2], weight);
	  }
	  else{
	    mean_projARU_p[iter][iter2]->Fill(rAlgImbProjAR_[setCorrNum][iter][iter2], weight);
	    mean_projAEtaU_p[iter][iter2]->Fill(rAlgImbProjAEta_[setCorrNum][iter][iter2], weight);
	    mean_projAPhiU_p[iter][iter2]->Fill(rAlgImbProjAPhi_[setCorrNum][iter][iter2], weight);
	    mean_projARCutU_p[iter][iter2]->Fill(rAlgImbProjAR_Cut_[setCorrNum][iter][iter2], weight);
	    mean_projARCutEtaU_p[iter][iter2]->Fill(rAlgImbProjAR_CutEta_[setCorrNum][iter][iter2], weight);
	    mean_projARCutPhiU_p[iter][iter2]->Fill(rAlgImbProjAR_CutPhi_[setCorrNum][iter][iter2], weight);
	    mean_projAEtaCutU_p[iter][iter2]->Fill(rAlgImbProjAEta_Cut_[setCorrNum][iter][iter2], weight);
	    mean_projAPhiCutU_p[iter][iter2]->Fill(rAlgImbProjAPhi_Cut_[setCorrNum][iter][iter2], weight);

	    mean_multARU_p[iter][iter2]->Fill(rAlgMultAR_[setCorrNum][iter][iter2], weight);
	    mean_multAEtaU_p[iter][iter2]->Fill(rAlgMultAEta_[setCorrNum][iter][iter2], weight);
	    mean_multAPhiU_p[iter][iter2]->Fill(rAlgMultAPhi_[setCorrNum][iter][iter2], weight);
	    mean_multARCutU_p[iter][iter2]->Fill(rAlgMultAR_Cut_[setCorrNum][iter][iter2], weight);
	    mean_multARCutEtaU_p[iter][iter2]->Fill(rAlgMultAR_CutEta_[setCorrNum][iter][iter2], weight);
	    mean_multARCutPhiU_p[iter][iter2]->Fill(rAlgMultAR_CutPhi_[setCorrNum][iter][iter2], weight);
	    mean_multAEtaCutU_p[iter][iter2]->Fill(rAlgMultAEta_Cut_[setCorrNum][iter][iter2], weight);
	    mean_multAPhiCutU_p[iter][iter2]->Fill(rAlgMultAPhi_Cut_[setCorrNum][iter][iter2], weight);
	  }
	}
	else{
	  mean_projAR_p[iter][iter2]->Fill(gAlgImbProjAR_[setNum][iter][iter2], weight);
	  mean_projAEta_p[iter][iter2]->Fill(gAlgImbProjAEta_[setNum][iter][iter2], weight);
	  mean_projAPhi_p[iter][iter2]->Fill(gAlgImbProjAPhi_[setNum][iter][iter2], weight);
	  mean_projARCut_p[iter][iter2]->Fill(gAlgImbProjAR_Cut_[setNum][iter][iter2], weight);
	  mean_projARCutEta_p[iter][iter2]->Fill(gAlgImbProjAR_CutEta_[setNum][iter][iter2], weight);
	  mean_projARCutPhi_p[iter][iter2]->Fill(gAlgImbProjAR_CutPhi_[setNum][iter][iter2], weight);
	  mean_projAEtaCut_p[iter][iter2]->Fill(gAlgImbProjAEta_Cut_[setNum][iter][iter2], weight);
	  mean_projAPhiCut_p[iter][iter2]->Fill(gAlgImbProjAPhi_Cut_[setNum][iter][iter2], weight);

	  mean_multAR_p[iter][iter2]->Fill(gAlgMultAR_[setNum][iter][iter2], weight);
	  mean_multAEta_p[iter][iter2]->Fill(gAlgMultAEta_[setNum][iter][iter2], weight);
	  mean_multAPhi_p[iter][iter2]->Fill(gAlgMultAPhi_[setNum][iter][iter2], weight);
	  mean_multARCut_p[iter][iter2]->Fill(gAlgMultAR_Cut_[setNum][iter][iter2], weight);
	  mean_multARCutEta_p[iter][iter2]->Fill(gAlgMultAR_CutEta_[setNum][iter][iter2], weight);
	  mean_multARCutPhi_p[iter][iter2]->Fill(gAlgMultAR_CutPhi_[setNum][iter][iter2], weight);
	  mean_multAEtaCut_p[iter][iter2]->Fill(gAlgMultAEta_Cut_[setNum][iter][iter2], weight);
	  mean_multAPhiCut_p[iter][iter2]->Fill(gAlgMultAPhi_Cut_[setNum][iter][iter2], weight);

	  if(AlgJtAsymm_[setNum] < 0.22){
	    mean_projARD_p[iter][iter2]->Fill(gAlgImbProjAR_[setNum][iter][iter2], weight);
	    mean_projAEtaD_p[iter][iter2]->Fill(gAlgImbProjAEta_[setNum][iter][iter2], weight);
	    mean_projAPhiD_p[iter][iter2]->Fill(gAlgImbProjAPhi_[setNum][iter][iter2], weight);
	    mean_projARCutD_p[iter][iter2]->Fill(gAlgImbProjAR_Cut_[setNum][iter][iter2], weight);
	    mean_projARCutEtaD_p[iter][iter2]->Fill(gAlgImbProjAR_CutEta_[setNum][iter][iter2], weight);
	    mean_projARCutPhiD_p[iter][iter2]->Fill(gAlgImbProjAR_CutPhi_[setNum][iter][iter2], weight);
	    mean_projAEtaCutD_p[iter][iter2]->Fill(gAlgImbProjAEta_Cut_[setNum][iter][iter2], weight);
	    mean_projAPhiCutD_p[iter][iter2]->Fill(gAlgImbProjAPhi_Cut_[setNum][iter][iter2], weight);

	    mean_multARD_p[iter][iter2]->Fill(gAlgMultAR_[setNum][iter][iter2], weight);
	    mean_multAEtaD_p[iter][iter2]->Fill(gAlgMultAEta_[setNum][iter][iter2], weight);
	    mean_multAPhiD_p[iter][iter2]->Fill(gAlgMultAPhi_[setNum][iter][iter2], weight);
	    mean_multARCutD_p[iter][iter2]->Fill(gAlgMultAR_Cut_[setNum][iter][iter2], weight);
	    mean_multARCutEtaD_p[iter][iter2]->Fill(gAlgMultAR_CutEta_[setNum][iter][iter2], weight);
	    mean_multARCutPhiD_p[iter][iter2]->Fill(gAlgMultAR_CutPhi_[setNum][iter][iter2], weight);
	    mean_multAEtaCutD_p[iter][iter2]->Fill(gAlgMultAEta_Cut_[setNum][iter][iter2], weight);
	    mean_multAPhiCutD_p[iter][iter2]->Fill(gAlgMultAPhi_Cut_[setNum][iter][iter2], weight);
	  }
	  else{
	    mean_projARU_p[iter][iter2]->Fill(gAlgImbProjAR_[setNum][iter][iter2], weight);
	    mean_projAEtaU_p[iter][iter2]->Fill(gAlgImbProjAEta_[setNum][iter][iter2], weight);
	    mean_projAPhiU_p[iter][iter2]->Fill(gAlgImbProjAPhi_[setNum][iter][iter2], weight);
	    mean_projARCutU_p[iter][iter2]->Fill(gAlgImbProjAR_Cut_[setNum][iter][iter2], weight);
	    mean_projARCutEtaU_p[iter][iter2]->Fill(gAlgImbProjAR_CutEta_[setNum][iter][iter2], weight);
	    mean_projARCutPhiU_p[iter][iter2]->Fill(gAlgImbProjAR_CutPhi_[setNum][iter][iter2], weight);
	    mean_projAEtaCutU_p[iter][iter2]->Fill(gAlgImbProjAEta_Cut_[setNum][iter][iter2], weight);
	    mean_projAPhiCutU_p[iter][iter2]->Fill(gAlgImbProjAPhi_Cut_[setNum][iter][iter2], weight);

	    mean_multARU_p[iter][iter2]->Fill(gAlgMultAR_[setNum][iter][iter2], weight);
	    mean_multAEtaU_p[iter][iter2]->Fill(gAlgMultAEta_[setNum][iter][iter2], weight);
	    mean_multAPhiU_p[iter][iter2]->Fill(gAlgMultAPhi_[setNum][iter][iter2], weight);
	    mean_multARCutU_p[iter][iter2]->Fill(gAlgMultAR_Cut_[setNum][iter][iter2], weight);
	    mean_multARCutEtaU_p[iter][iter2]->Fill(gAlgMultAR_CutEta_[setNum][iter][iter2], weight);
	    mean_multARCutPhiU_p[iter][iter2]->Fill(gAlgMultAR_CutPhi_[setNum][iter][iter2], weight);
	    mean_multAEtaCutU_p[iter][iter2]->Fill(gAlgMultAEta_Cut_[setNum][iter][iter2], weight);
	    mean_multAPhiCutU_p[iter][iter2]->Fill(gAlgMultAPhi_Cut_[setNum][iter][iter2], weight);
	  }
	}
	
      }
    }
  }


  for(Int_t jEntry = 0; jEntry < (Int_t)anaTree_p->GetEntries(); jEntry++){
    anaTree_p->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(isEventCut(setNum, sType, centLow, centHi, isHighPtTrk)) continue;

    if(TMath::Abs(AlgJtEta_[setNum][0]) > 1.6 || TMath::Abs(AlgJtEta_[setNum][1]) > 1.6) continue;

    Float_t weight = 1.0;
    if(montecarlo){
      weight = pthatWeight_;
      if(isHI(sType)) weight *= centWeight_[setNum];
    }


    for(Int_t iter = 0; iter < 6; iter++){
      for(Int_t iter2 = 0; iter2 < nBins; iter2++){
        if(strcmp(gR.c_str(), "r") == 0){
          mean_projARMID_p[iter][iter2]->Fill(rAlgImbProjAR_[setCorrNum][iter][iter2], weight);

	  if(AlgJtAsymm_[setNum] < 0.22) mean_projARMIDD_p[iter][iter2]->Fill(rAlgImbProjAR_[setCorrNum][iter][iter2], weight);
	  else mean_projARMIDU_p[iter][iter2]->Fill(rAlgImbProjAR_[setCorrNum][iter][iter2], weight);
	}
	else if(strcmp(gR.c_str(), "g") == 0){
          mean_projARMID_p[iter][iter2]->Fill(gAlgImbProjAR_[setNum][iter][iter2], weight);

	  if(AlgJtAsymm_[setNum] < 0.22) mean_projARMIDD_p[iter][iter2]->Fill(gAlgImbProjAR_[setNum][iter][iter2], weight);
	  else mean_projARMIDU_p[iter][iter2]->Fill(gAlgImbProjAR_[setNum][iter][iter2], weight);
	}
      }
    }

    if((AlgJtEta_[setNum][0] <  midRapCut && AlgJtEta_[setNum][1] < midRapCut) || (AlgJtEta_[setNum][0] > -midRapCut && AlgJtEta_[setNum][1] > -midRapCut)){
      for(Int_t iter = 0; iter < 6; iter++){
	for(Int_t iter2 = 0; iter2 < nBins; iter2++){
	  if(strcmp(gR.c_str(), "r") == 0){
	    mean_projARFOR_p[iter][iter2]->Fill(rAlgImbProjAR_FOR_[setCorrNum][iter][iter2], weight);

	    if(AlgJtAsymm_[setNum] < 0.22) mean_projARFORD_p[iter][iter2]->Fill(rAlgImbProjAR_FOR_[setCorrNum][iter][iter2], weight);
	    else mean_projARFORU_p[iter][iter2]->Fill(rAlgImbProjAR_FOR_[setCorrNum][iter][iter2], weight);

	    if(TMath::Abs(AlgJtEta_[setNum][0]) < midRapCut && TMath::Abs(AlgJtEta_[setNum][1]) < midRapCut){
	      mean_projARFORMID_p[iter][iter2]->Fill(rAlgImbProjAR_FOR_[setCorrNum][iter][iter2], weight);
	      
	      if(AlgJtAsymm_[setNum] < 0.22) mean_projARFORMIDD_p[iter][iter2]->Fill(rAlgImbProjAR_FOR_[setCorrNum][iter][iter2], weight);
	      else mean_projARFORMIDU_p[iter][iter2]->Fill(rAlgImbProjAR_FOR_[setCorrNum][iter][iter2], weight);
	    }
	    else{
	      mean_projARFORFOR_p[iter][iter2]->Fill(rAlgImbProjAR_FOR_[setCorrNum][iter][iter2], weight);
	      
	      if(AlgJtAsymm_[setNum] < 0.22) mean_projARFORFORD_p[iter][iter2]->Fill(rAlgImbProjAR_FOR_[setCorrNum][iter][iter2], weight);
	      else mean_projARFORFORU_p[iter][iter2]->Fill(rAlgImbProjAR_FOR_[setCorrNum][iter][iter2], weight);
	    }
	  }
	  else if(strcmp(gR.c_str(), "g") == 0){
	    mean_projARFOR_p[iter][iter2]->Fill(gAlgImbProjAR_FOR_[setNum][iter][iter2], weight);
	    
	    if(AlgJtAsymm_[setNum] < 0.22) mean_projARFORD_p[iter][iter2]->Fill(gAlgImbProjAR_FOR_[setNum][iter][iter2], weight);
	    else mean_projARFORU_p[iter][iter2]->Fill(gAlgImbProjAR_FOR_[setNum][iter][iter2], weight);


	    if(TMath::Abs(AlgJtEta_[setNum][0]) < midRapCut && TMath::Abs(AlgJtEta_[setNum][1]) < midRapCut){
	      mean_projARFORMID_p[iter][iter2]->Fill(gAlgImbProjAR_FOR_[setNum][iter][iter2], weight);
	      
	      if(AlgJtAsymm_[setNum] < 0.22) mean_projARFORMIDD_p[iter][iter2]->Fill(gAlgImbProjAR_FOR_[setNum][iter][iter2], weight);
	      else mean_projARFORMIDU_p[iter][iter2]->Fill(gAlgImbProjAR_FOR_[setNum][iter][iter2], weight);
	    }
	    else{
	      mean_projARFORFOR_p[iter][iter2]->Fill(gAlgImbProjAR_FOR_[setNum][iter][iter2], weight);
	      
	      if(AlgJtAsymm_[setNum] < 0.22) mean_projARFORFORD_p[iter][iter2]->Fill(gAlgImbProjAR_FOR_[setNum][iter][iter2], weight);
	      else mean_projARFORFORU_p[iter][iter2]->Fill(gAlgImbProjAR_FOR_[setNum][iter][iter2], weight);
	    }
	  }
	}
      }
    }
  }


  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      if(mean_projAR_p[iter][iter2]->GetEntries() != 0){
	imbARHist_p[iter]->SetBinContent(iter2+1, mean_projAR_p[iter][iter2]->GetMean());
	imbARHist_p[iter]->SetBinError(iter2+1, mean_projAR_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARD_p[iter][iter2]->GetEntries() != 0){
	imbARDHist_p[iter]->SetBinContent(iter2+1, mean_projARD_p[iter][iter2]->GetMean());
	imbARDHist_p[iter]->SetBinError(iter2+1, mean_projARD_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARU_p[iter][iter2]->GetEntries() != 0){
	imbARUHist_p[iter]->SetBinContent(iter2+1, mean_projARU_p[iter][iter2]->GetMean());
	imbARUHist_p[iter]->SetBinError(iter2+1, mean_projARU_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARCut_p[iter][iter2]->GetEntries() != 0){
	imbARCutHist_p[iter]->SetBinContent(iter2+1, mean_projARCut_p[iter][iter2]->GetMean());
	imbARCutHist_p[iter]->SetBinError(iter2+1, mean_projARCut_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARCutD_p[iter][iter2]->GetEntries() != 0){
	imbARCutDHist_p[iter]->SetBinContent(iter2+1, mean_projARCutD_p[iter][iter2]->GetMean());
	imbARCutDHist_p[iter]->SetBinError(iter2+1, mean_projARCutD_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARCutU_p[iter][iter2]->GetEntries() != 0){
	imbARCutUHist_p[iter]->SetBinContent(iter2+1, mean_projARCutU_p[iter][iter2]->GetMean());
	imbARCutUHist_p[iter]->SetBinError(iter2+1, mean_projARCutU_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARCutEta_p[iter][iter2]->GetEntries() != 0){
	imbARCutEtaHist_p[iter]->SetBinContent(iter2+1, mean_projARCutEta_p[iter][iter2]->GetMean());
	imbARCutEtaHist_p[iter]->SetBinError(iter2+1, mean_projARCutEta_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARCutEtaD_p[iter][iter2]->GetEntries() != 0){
	imbARCutEtaDHist_p[iter]->SetBinContent(iter2+1, mean_projARCutEtaD_p[iter][iter2]->GetMean());
	imbARCutEtaDHist_p[iter]->SetBinError(iter2+1, mean_projARCutEtaD_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARCutEtaU_p[iter][iter2]->GetEntries() != 0){
	imbARCutEtaUHist_p[iter]->SetBinContent(iter2+1, mean_projARCutEtaU_p[iter][iter2]->GetMean());
	imbARCutEtaUHist_p[iter]->SetBinError(iter2+1, mean_projARCutEtaU_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARCutPhi_p[iter][iter2]->GetEntries() != 0){
	imbARCutPhiHist_p[iter]->SetBinContent(iter2+1, mean_projARCutPhi_p[iter][iter2]->GetMean());
	imbARCutPhiHist_p[iter]->SetBinError(iter2+1, mean_projARCutPhi_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARCutPhiD_p[iter][iter2]->GetEntries() != 0){
	imbARCutPhiDHist_p[iter]->SetBinContent(iter2+1, mean_projARCutPhiD_p[iter][iter2]->GetMean());
	imbARCutPhiDHist_p[iter]->SetBinError(iter2+1, mean_projARCutPhiD_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARCutPhiU_p[iter][iter2]->GetEntries() != 0){
	imbARCutPhiUHist_p[iter]->SetBinContent(iter2+1, mean_projARCutPhiU_p[iter][iter2]->GetMean());
	imbARCutPhiUHist_p[iter]->SetBinError(iter2+1, mean_projARCutPhiU_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAEta_p[iter][iter2]->GetEntries() != 0){
	imbAEtaHist_p[iter]->SetBinContent(iter2+1, mean_projAEta_p[iter][iter2]->GetMean());
	imbAEtaHist_p[iter]->SetBinError(iter2+1, mean_projAEta_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAEtaD_p[iter][iter2]->GetEntries() != 0){
	imbAEtaDHist_p[iter]->SetBinContent(iter2+1, mean_projAEtaD_p[iter][iter2]->GetMean());
	imbAEtaDHist_p[iter]->SetBinError(iter2+1, mean_projAEtaD_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAEtaU_p[iter][iter2]->GetEntries() != 0){
	imbAEtaUHist_p[iter]->SetBinContent(iter2+1, mean_projAEtaU_p[iter][iter2]->GetMean());
	imbAEtaUHist_p[iter]->SetBinError(iter2+1, mean_projAEtaU_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAEtaCut_p[iter][iter2]->GetEntries() != 0){
	imbAEtaCutHist_p[iter]->SetBinContent(iter2+1, mean_projAEtaCut_p[iter][iter2]->GetMean());
	imbAEtaCutHist_p[iter]->SetBinError(iter2+1, mean_projAEtaCut_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAEtaCutD_p[iter][iter2]->GetEntries() != 0){
	imbAEtaCutDHist_p[iter]->SetBinContent(iter2+1, mean_projAEtaCutD_p[iter][iter2]->GetMean());
	imbAEtaCutDHist_p[iter]->SetBinError(iter2+1, mean_projAEtaCutD_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAEtaCutU_p[iter][iter2]->GetEntries() != 0){
	imbAEtaCutUHist_p[iter]->SetBinContent(iter2+1, mean_projAEtaCutU_p[iter][iter2]->GetMean());
	imbAEtaCutUHist_p[iter]->SetBinError(iter2+1, mean_projAEtaCutU_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAPhi_p[iter][iter2]->GetEntries() != 0){
	imbAPhiHist_p[iter]->SetBinContent(iter2+1, mean_projAPhi_p[iter][iter2]->GetMean());
	imbAPhiHist_p[iter]->SetBinError(iter2+1, mean_projAPhi_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAPhiD_p[iter][iter2]->GetEntries() != 0){
	imbAPhiDHist_p[iter]->SetBinContent(iter2+1, mean_projAPhiD_p[iter][iter2]->GetMean());
	imbAPhiDHist_p[iter]->SetBinError(iter2+1, mean_projAPhiD_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAPhiU_p[iter][iter2]->GetEntries() != 0){
	imbAPhiUHist_p[iter]->SetBinContent(iter2+1, mean_projAPhiU_p[iter][iter2]->GetMean());
	imbAPhiUHist_p[iter]->SetBinError(iter2+1, mean_projAPhiU_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAPhiCut_p[iter][iter2]->GetEntries() != 0){
	imbAPhiCutHist_p[iter]->SetBinContent(iter2+1, mean_projAPhiCut_p[iter][iter2]->GetMean());
	imbAPhiCutHist_p[iter]->SetBinError(iter2+1, mean_projAPhiCut_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAPhiCutD_p[iter][iter2]->GetEntries() != 0){
	imbAPhiCutDHist_p[iter]->SetBinContent(iter2+1, mean_projAPhiCutD_p[iter][iter2]->GetMean());
	imbAPhiCutDHist_p[iter]->SetBinError(iter2+1, mean_projAPhiCutD_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAPhiCutU_p[iter][iter2]->GetEntries() != 0){
	imbAPhiCutUHist_p[iter]->SetBinContent(iter2+1, mean_projAPhiCutU_p[iter][iter2]->GetMean());
	imbAPhiCutUHist_p[iter]->SetBinError(iter2+1, mean_projAPhiCutU_p[iter][iter2]->GetMeanError());
      }



      if(mean_projARFOR_p[iter][iter2]->GetEntries() != 0){
	imbARFORHist_p[iter]->SetBinContent(iter2+1, mean_projARFOR_p[iter][iter2]->GetMean());
	imbARFORHist_p[iter]->SetBinError(iter2+1, mean_projARFOR_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARFORD_p[iter][iter2]->GetEntries() != 0){
	imbARFORDHist_p[iter]->SetBinContent(iter2+1, mean_projARFORD_p[iter][iter2]->GetMean());
	imbARFORDHist_p[iter]->SetBinError(iter2+1, mean_projARFORD_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARFORU_p[iter][iter2]->GetEntries() != 0){
	imbARFORUHist_p[iter]->SetBinContent(iter2+1, mean_projARFORU_p[iter][iter2]->GetMean());
	imbARFORUHist_p[iter]->SetBinError(iter2+1, mean_projARFORU_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARFORMID_p[iter][iter2]->GetEntries() != 0){
	imbARFORMIDHist_p[iter]->SetBinContent(iter2+1, mean_projARFORMID_p[iter][iter2]->GetMean());
	imbARFORMIDHist_p[iter]->SetBinError(iter2+1, mean_projARFORMID_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARFORMIDD_p[iter][iter2]->GetEntries() != 0){
	imbARFORMIDDHist_p[iter]->SetBinContent(iter2+1, mean_projARFORMIDD_p[iter][iter2]->GetMean());
	imbARFORMIDDHist_p[iter]->SetBinError(iter2+1, mean_projARFORMIDD_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARFORMIDU_p[iter][iter2]->GetEntries() != 0){
	imbARFORMIDUHist_p[iter]->SetBinContent(iter2+1, mean_projARFORMIDU_p[iter][iter2]->GetMean());
	imbARFORMIDUHist_p[iter]->SetBinError(iter2+1, mean_projARFORMIDU_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARFORFOR_p[iter][iter2]->GetEntries() != 0){
	imbARFORFORHist_p[iter]->SetBinContent(iter2+1, mean_projARFORFOR_p[iter][iter2]->GetMean());
	imbARFORFORHist_p[iter]->SetBinError(iter2+1, mean_projARFORFOR_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARFORFORD_p[iter][iter2]->GetEntries() != 0){
	imbARFORFORDHist_p[iter]->SetBinContent(iter2+1, mean_projARFORFORD_p[iter][iter2]->GetMean());
	imbARFORFORDHist_p[iter]->SetBinError(iter2+1, mean_projARFORFORD_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARFORFORU_p[iter][iter2]->GetEntries() != 0){
	imbARFORFORUHist_p[iter]->SetBinContent(iter2+1, mean_projARFORFORU_p[iter][iter2]->GetMean());
	imbARFORFORUHist_p[iter]->SetBinError(iter2+1, mean_projARFORFORU_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARMID_p[iter][iter2]->GetEntries() != 0){
	imbARMIDHist_p[iter]->SetBinContent(iter2+1, mean_projARMID_p[iter][iter2]->GetMean());
	imbARMIDHist_p[iter]->SetBinError(iter2+1, mean_projARMID_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARMIDD_p[iter][iter2]->GetEntries() != 0){
	imbARMIDDHist_p[iter]->SetBinContent(iter2+1, mean_projARMIDD_p[iter][iter2]->GetMean());
	imbARMIDDHist_p[iter]->SetBinError(iter2+1, mean_projARMIDD_p[iter][iter2]->GetMeanError());
      }
      if(mean_projARMIDU_p[iter][iter2]->GetEntries() != 0){
	imbARMIDUHist_p[iter]->SetBinContent(iter2+1, mean_projARMIDU_p[iter][iter2]->GetMean());
	imbARMIDUHist_p[iter]->SetBinError(iter2+1, mean_projARMIDU_p[iter][iter2]->GetMeanError());
      }


      //MULTIPLICITY PROJECTIONS

      if(mean_multAR_p[iter][iter2]->GetEntries() != 0){
	multARHist_p[iter]->SetBinContent(iter2+1, mean_multAR_p[iter][iter2]->GetMean());
	multARHist_p[iter]->SetBinError(iter2+1, mean_multAR_p[iter][iter2]->GetMeanError());
      }
      if(mean_multARD_p[iter][iter2]->GetEntries() != 0){
	multARDHist_p[iter]->SetBinContent(iter2+1, mean_multARD_p[iter][iter2]->GetMean());
	multARDHist_p[iter]->SetBinError(iter2+1, mean_multARD_p[iter][iter2]->GetMeanError());
      }
      if(mean_multARU_p[iter][iter2]->GetEntries() != 0){
	multARUHist_p[iter]->SetBinContent(iter2+1, mean_multARU_p[iter][iter2]->GetMean());
	multARUHist_p[iter]->SetBinError(iter2+1, mean_multARU_p[iter][iter2]->GetMeanError());
      }
      if(mean_multARCut_p[iter][iter2]->GetEntries() != 0){
	multARCutHist_p[iter]->SetBinContent(iter2+1, mean_multARCut_p[iter][iter2]->GetMean());
	multARCutHist_p[iter]->SetBinError(iter2+1, mean_multARCut_p[iter][iter2]->GetMeanError());
      }
      if(mean_multARCutD_p[iter][iter2]->GetEntries() != 0){
	multARCutDHist_p[iter]->SetBinContent(iter2+1, mean_multARCutD_p[iter][iter2]->GetMean());
	multARCutDHist_p[iter]->SetBinError(iter2+1, mean_multARCutD_p[iter][iter2]->GetMeanError());
      }
      if(mean_multARCutU_p[iter][iter2]->GetEntries() != 0){
	multARCutUHist_p[iter]->SetBinContent(iter2+1, mean_multARCutU_p[iter][iter2]->GetMean());
	multARCutUHist_p[iter]->SetBinError(iter2+1, mean_multARCutU_p[iter][iter2]->GetMeanError());
      }
      if(mean_multARCutEta_p[iter][iter2]->GetEntries() != 0){
	multARCutEtaHist_p[iter]->SetBinContent(iter2+1, mean_multARCutEta_p[iter][iter2]->GetMean());
	multARCutEtaHist_p[iter]->SetBinError(iter2+1, mean_multARCutEta_p[iter][iter2]->GetMeanError());
      }
      if(mean_multARCutEtaD_p[iter][iter2]->GetEntries() != 0){
	multARCutEtaDHist_p[iter]->SetBinContent(iter2+1, mean_multARCutEtaD_p[iter][iter2]->GetMean());
	multARCutEtaDHist_p[iter]->SetBinError(iter2+1, mean_multARCutEtaD_p[iter][iter2]->GetMeanError());
      }
      if(mean_multARCutEtaU_p[iter][iter2]->GetEntries() != 0){
	multARCutEtaUHist_p[iter]->SetBinContent(iter2+1, mean_multARCutEtaU_p[iter][iter2]->GetMean());
	multARCutEtaUHist_p[iter]->SetBinError(iter2+1, mean_multARCutEtaU_p[iter][iter2]->GetMeanError());
      }
      if(mean_multARCutPhi_p[iter][iter2]->GetEntries() != 0){
	multARCutPhiHist_p[iter]->SetBinContent(iter2+1, mean_multARCutPhi_p[iter][iter2]->GetMean());
	multARCutPhiHist_p[iter]->SetBinError(iter2+1, mean_multARCutPhi_p[iter][iter2]->GetMeanError());
      }
      if(mean_multARCutPhiD_p[iter][iter2]->GetEntries() != 0){
	multARCutPhiDHist_p[iter]->SetBinContent(iter2+1, mean_multARCutPhiD_p[iter][iter2]->GetMean());
	multARCutPhiDHist_p[iter]->SetBinError(iter2+1, mean_multARCutPhiD_p[iter][iter2]->GetMeanError());
      }
      if(mean_multARCutPhiU_p[iter][iter2]->GetEntries() != 0){
	multARCutPhiUHist_p[iter]->SetBinContent(iter2+1, mean_multARCutPhiU_p[iter][iter2]->GetMean());
	multARCutPhiUHist_p[iter]->SetBinError(iter2+1, mean_multARCutPhiU_p[iter][iter2]->GetMeanError());
      }
      if(mean_multAEta_p[iter][iter2]->GetEntries() != 0){
	multAEtaHist_p[iter]->SetBinContent(iter2+1, mean_multAEta_p[iter][iter2]->GetMean());
	multAEtaHist_p[iter]->SetBinError(iter2+1, mean_multAEta_p[iter][iter2]->GetMeanError());
      }
      if(mean_multAEtaD_p[iter][iter2]->GetEntries() != 0){
	multAEtaDHist_p[iter]->SetBinContent(iter2+1, mean_multAEtaD_p[iter][iter2]->GetMean());
	multAEtaDHist_p[iter]->SetBinError(iter2+1, mean_multAEtaD_p[iter][iter2]->GetMeanError());
      }
      if(mean_multAEtaU_p[iter][iter2]->GetEntries() != 0){
	multAEtaUHist_p[iter]->SetBinContent(iter2+1, mean_multAEtaU_p[iter][iter2]->GetMean());
	multAEtaUHist_p[iter]->SetBinError(iter2+1, mean_multAEtaU_p[iter][iter2]->GetMeanError());
      }
      if(mean_multAEtaCut_p[iter][iter2]->GetEntries() != 0){
	multAEtaCutHist_p[iter]->SetBinContent(iter2+1, mean_multAEtaCut_p[iter][iter2]->GetMean());
	multAEtaCutHist_p[iter]->SetBinError(iter2+1, mean_multAEtaCut_p[iter][iter2]->GetMeanError());
      }
      if(mean_multAEtaCutD_p[iter][iter2]->GetEntries() != 0){
	multAEtaCutDHist_p[iter]->SetBinContent(iter2+1, mean_multAEtaCutD_p[iter][iter2]->GetMean());
	multAEtaCutDHist_p[iter]->SetBinError(iter2+1, mean_multAEtaCutD_p[iter][iter2]->GetMeanError());
      }
      if(mean_multAEtaCutU_p[iter][iter2]->GetEntries() != 0){
	multAEtaCutUHist_p[iter]->SetBinContent(iter2+1, mean_multAEtaCutU_p[iter][iter2]->GetMean());
	multAEtaCutUHist_p[iter]->SetBinError(iter2+1, mean_multAEtaCutU_p[iter][iter2]->GetMeanError());
      }
      if(mean_multAPhi_p[iter][iter2]->GetEntries() != 0){
	multAPhiHist_p[iter]->SetBinContent(iter2+1, mean_multAPhi_p[iter][iter2]->GetMean());
	multAPhiHist_p[iter]->SetBinError(iter2+1, mean_multAPhi_p[iter][iter2]->GetMeanError());
      }
      if(mean_multAPhiD_p[iter][iter2]->GetEntries() != 0){
	multAPhiDHist_p[iter]->SetBinContent(iter2+1, mean_multAPhiD_p[iter][iter2]->GetMean());
	multAPhiDHist_p[iter]->SetBinError(iter2+1, mean_multAPhiD_p[iter][iter2]->GetMeanError());
      }
      if(mean_multAPhiU_p[iter][iter2]->GetEntries() != 0){
	multAPhiUHist_p[iter]->SetBinContent(iter2+1, mean_multAPhiU_p[iter][iter2]->GetMean());
	multAPhiUHist_p[iter]->SetBinError(iter2+1, mean_multAPhiU_p[iter][iter2]->GetMeanError());
      }
      if(mean_multAPhiCut_p[iter][iter2]->GetEntries() != 0){
	multAPhiCutHist_p[iter]->SetBinContent(iter2+1, mean_multAPhiCut_p[iter][iter2]->GetMean());
	multAPhiCutHist_p[iter]->SetBinError(iter2+1, mean_multAPhiCut_p[iter][iter2]->GetMeanError());
      }
      if(mean_multAPhiCutD_p[iter][iter2]->GetEntries() != 0){
	multAPhiCutDHist_p[iter]->SetBinContent(iter2+1, mean_multAPhiCutD_p[iter][iter2]->GetMean());
	multAPhiCutDHist_p[iter]->SetBinError(iter2+1, mean_multAPhiCutD_p[iter][iter2]->GetMeanError());
      }
      if(mean_multAPhiCutU_p[iter][iter2]->GetEntries() != 0){
	multAPhiCutUHist_p[iter]->SetBinContent(iter2+1, mean_multAPhiCutU_p[iter][iter2]->GetMean());
	multAPhiCutUHist_p[iter]->SetBinError(iter2+1, mean_multAPhiCutU_p[iter][iter2]->GetMeanError());
      }
    }
  }

  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins2; iter2++){
      if(mean_projAR2_p[iter][iter2]->GetEntries() != 0){
	imbAR2Hist_p[iter]->SetBinContent(iter2+1, mean_projAR2_p[iter][iter2]->GetMean());
	imbAR2Hist_p[iter]->SetBinError(iter2+1, mean_projAR2_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAR2D_p[iter][iter2]->GetEntries() != 0){
	imbAR2DHist_p[iter]->SetBinContent(iter2+1, mean_projAR2D_p[iter][iter2]->GetMean());
	imbAR2DHist_p[iter]->SetBinError(iter2+1, mean_projAR2D_p[iter][iter2]->GetMeanError());
      }
      if(mean_projAR2U_p[iter][iter2]->GetEntries() != 0){
	imbAR2UHist_p[iter]->SetBinContent(iter2+1, mean_projAR2U_p[iter][iter2]->GetMean());
	imbAR2UHist_p[iter]->SetBinError(iter2+1, mean_projAR2U_p[iter][iter2]->GetMeanError());
      }
    }
  }


  outFile_p = new TFile(Form("%s.root", outName.c_str()), "UPDATE");
  std::cout << outName << std::endl;

  for(Int_t iter = 0; iter < 6; iter++){
    imbARHist_p[iter]->Write("", TObject::kOverwrite);
    imbARDHist_p[iter]->Write("", TObject::kOverwrite);
    imbARUHist_p[iter]->Write("", TObject::kOverwrite);

    imbAR2Hist_p[iter]->Write("", TObject::kOverwrite);
    imbAR2DHist_p[iter]->Write("", TObject::kOverwrite);
    imbAR2UHist_p[iter]->Write("", TObject::kOverwrite);

    imbARCutHist_p[iter]->Write("", TObject::kOverwrite);
    imbARCutDHist_p[iter]->Write("", TObject::kOverwrite);
    imbARCutUHist_p[iter]->Write("", TObject::kOverwrite);
    imbARCutEtaHist_p[iter]->Write("", TObject::kOverwrite);
    imbARCutEtaDHist_p[iter]->Write("", TObject::kOverwrite);
    imbARCutEtaUHist_p[iter]->Write("", TObject::kOverwrite);
    imbARCutPhiHist_p[iter]->Write("", TObject::kOverwrite);
    imbARCutPhiDHist_p[iter]->Write("", TObject::kOverwrite);
    imbARCutPhiUHist_p[iter]->Write("", TObject::kOverwrite);
    imbAEtaHist_p[iter]->Write("", TObject::kOverwrite);
    imbAEtaDHist_p[iter]->Write("", TObject::kOverwrite);
    imbAEtaUHist_p[iter]->Write("", TObject::kOverwrite);
    imbAEtaCutHist_p[iter]->Write("", TObject::kOverwrite);
    imbAEtaCutDHist_p[iter]->Write("", TObject::kOverwrite);
    imbAEtaCutUHist_p[iter]->Write("", TObject::kOverwrite);
    imbAPhiHist_p[iter]->Write("", TObject::kOverwrite);
    imbAPhiDHist_p[iter]->Write("", TObject::kOverwrite);
    imbAPhiUHist_p[iter]->Write("", TObject::kOverwrite);
    imbAPhiCutHist_p[iter]->Write("", TObject::kOverwrite);
    imbAPhiCutDHist_p[iter]->Write("", TObject::kOverwrite);
    imbAPhiCutUHist_p[iter]->Write("", TObject::kOverwrite);


    imbARFORHist_p[iter]->Write("", TObject::kOverwrite);
    imbARFORDHist_p[iter]->Write("", TObject::kOverwrite);
    imbARFORUHist_p[iter]->Write("", TObject::kOverwrite);
    imbARFORMIDHist_p[iter]->Write("", TObject::kOverwrite);
    imbARFORMIDDHist_p[iter]->Write("", TObject::kOverwrite);
    imbARFORMIDUHist_p[iter]->Write("", TObject::kOverwrite);
    imbARFORFORHist_p[iter]->Write("", TObject::kOverwrite);
    imbARFORFORDHist_p[iter]->Write("", TObject::kOverwrite);
    imbARFORFORUHist_p[iter]->Write("", TObject::kOverwrite);
    imbARMIDHist_p[iter]->Write("", TObject::kOverwrite);
    imbARMIDDHist_p[iter]->Write("", TObject::kOverwrite);
    imbARMIDUHist_p[iter]->Write("", TObject::kOverwrite);


    multARHist_p[iter]->Write("", TObject::kOverwrite);
    multARDHist_p[iter]->Write("", TObject::kOverwrite);
    multARUHist_p[iter]->Write("", TObject::kOverwrite);
    multARCutHist_p[iter]->Write("", TObject::kOverwrite);
    multARCutDHist_p[iter]->Write("", TObject::kOverwrite);
    multARCutUHist_p[iter]->Write("", TObject::kOverwrite);
    multARCutEtaHist_p[iter]->Write("", TObject::kOverwrite);
    multARCutEtaDHist_p[iter]->Write("", TObject::kOverwrite);
    multARCutEtaUHist_p[iter]->Write("", TObject::kOverwrite);
    multARCutPhiHist_p[iter]->Write("", TObject::kOverwrite);
    multARCutPhiDHist_p[iter]->Write("", TObject::kOverwrite);
    multARCutPhiUHist_p[iter]->Write("", TObject::kOverwrite);
    multAEtaHist_p[iter]->Write("", TObject::kOverwrite);
    multAEtaDHist_p[iter]->Write("", TObject::kOverwrite);
    multAEtaUHist_p[iter]->Write("", TObject::kOverwrite);
    multAEtaCutHist_p[iter]->Write("", TObject::kOverwrite);
    multAEtaCutDHist_p[iter]->Write("", TObject::kOverwrite);
    multAEtaCutUHist_p[iter]->Write("", TObject::kOverwrite);
    multAPhiHist_p[iter]->Write("", TObject::kOverwrite);
    multAPhiDHist_p[iter]->Write("", TObject::kOverwrite);
    multAPhiUHist_p[iter]->Write("", TObject::kOverwrite);
    multAPhiCutHist_p[iter]->Write("", TObject::kOverwrite);
    multAPhiCutDHist_p[iter]->Write("", TObject::kOverwrite);
    multAPhiCutUHist_p[iter]->Write("", TObject::kOverwrite);
  } 

  outFile_p->Close();
  delete outFile_p;
  outFile_p = 0;

  if(!montecarlo){
    outFile_p = new TFile("ajIncl_DelR_allPtBins.root", "UPDATE");
    for(Int_t iter = 0; iter < 6; iter++){
      imbARHist_p[iter]->Write("", TObject::kOverwrite);
    }
    outFile_p->Close();
    delete outFile_p;
    outFile_p = 0;
  }

  std::ofstream txtFile;
  txtFile.open(Form("%s.txt", outName.c_str()), std::ofstream::out | std::ofstream::app);
  if(centLow == 0) txtFile << Form("ProjAR (D) (U) %s,\n", algType[setNum].c_str());
  txtFile << Form("   %d-%d: %d (%d) (%d).\n", centLow, centHi, evtsPass, evtsPassD, evtsPassU);
  txtFile.close();

  for(Int_t iter = 0; iter < 6; iter++){
    CleanHist(imbARHist_p[iter]);
    CleanHist(imbARDHist_p[iter]);
    CleanHist(imbARUHist_p[iter]);

    CleanHist(imbAR2Hist_p[iter]);
    CleanHist(imbAR2DHist_p[iter]);
    CleanHist(imbAR2UHist_p[iter]);

    CleanHist(imbARCutHist_p[iter]);
    CleanHist(imbARCutDHist_p[iter]);
    CleanHist(imbARCutUHist_p[iter]);
    CleanHist(imbARCutEtaHist_p[iter]);
    CleanHist(imbARCutEtaDHist_p[iter]);
    CleanHist(imbARCutEtaUHist_p[iter]);
    CleanHist(imbARCutPhiHist_p[iter]);
    CleanHist(imbARCutPhiDHist_p[iter]);
    CleanHist(imbARCutPhiUHist_p[iter]);
    CleanHist(imbAEtaHist_p[iter]);
    CleanHist(imbAEtaDHist_p[iter]);
    CleanHist(imbAEtaUHist_p[iter]);
    CleanHist(imbAEtaCutHist_p[iter]);
    CleanHist(imbAEtaCutDHist_p[iter]);
    CleanHist(imbAEtaCutUHist_p[iter]);
    CleanHist(imbAPhiHist_p[iter]);
    CleanHist(imbAPhiDHist_p[iter]);
    CleanHist(imbAPhiUHist_p[iter]);
    CleanHist(imbAPhiCutHist_p[iter]);
    CleanHist(imbAPhiCutDHist_p[iter]);
    CleanHist(imbAPhiCutUHist_p[iter]);

    CleanHist(imbARFORHist_p[iter]);
    CleanHist(imbARFORDHist_p[iter]);
    CleanHist(imbARFORUHist_p[iter]);
    CleanHist(imbARFORMIDHist_p[iter]);
    CleanHist(imbARFORMIDDHist_p[iter]);
    CleanHist(imbARFORMIDUHist_p[iter]);
    CleanHist(imbARFORFORHist_p[iter]);
    CleanHist(imbARFORFORDHist_p[iter]);
    CleanHist(imbARFORFORUHist_p[iter]);
    CleanHist(imbARMIDHist_p[iter]);
    CleanHist(imbARMIDDHist_p[iter]);
    CleanHist(imbARMIDUHist_p[iter]);

    CleanHist(multARHist_p[iter]);
    CleanHist(multARDHist_p[iter]);
    CleanHist(multARUHist_p[iter]);
    CleanHist(multARCutHist_p[iter]);
    CleanHist(multARCutDHist_p[iter]);
    CleanHist(multARCutUHist_p[iter]);
    CleanHist(multARCutEtaHist_p[iter]);
    CleanHist(multARCutEtaDHist_p[iter]);
    CleanHist(multARCutEtaUHist_p[iter]);
    CleanHist(multARCutPhiHist_p[iter]);
    CleanHist(multARCutPhiDHist_p[iter]);
    CleanHist(multARCutPhiUHist_p[iter]);
    CleanHist(multAEtaHist_p[iter]);
    CleanHist(multAEtaDHist_p[iter]);
    CleanHist(multAEtaUHist_p[iter]);
    CleanHist(multAEtaCutHist_p[iter]);
    CleanHist(multAEtaCutDHist_p[iter]);
    CleanHist(multAEtaCutUHist_p[iter]);
    CleanHist(multAPhiHist_p[iter]);
    CleanHist(multAPhiDHist_p[iter]);
    CleanHist(multAPhiUHist_p[iter]);
    CleanHist(multAPhiCutHist_p[iter]);
    CleanHist(multAPhiCutDHist_p[iter]);
    CleanHist(multAPhiCutUHist_p[iter]);


    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      CleanHist(mean_projAR_p[iter][iter2]);
      CleanHist(mean_projARD_p[iter][iter2]);
      CleanHist(mean_projARU_p[iter][iter2]);

      if(iter2 < nBins2){
	CleanHist(mean_projAR2_p[iter][iter2]);
	CleanHist(mean_projAR2D_p[iter][iter2]);
	CleanHist(mean_projAR2U_p[iter][iter2]);
      }

      CleanHist(mean_projARCut_p[iter][iter2]);
      CleanHist(mean_projARCutD_p[iter][iter2]);
      CleanHist(mean_projARCutU_p[iter][iter2]);
      CleanHist(mean_projARCutEta_p[iter][iter2]);
      CleanHist(mean_projARCutEtaD_p[iter][iter2]);
      CleanHist(mean_projARCutEtaU_p[iter][iter2]);
      CleanHist(mean_projARCutPhi_p[iter][iter2]);
      CleanHist(mean_projARCutPhiD_p[iter][iter2]);
      CleanHist(mean_projARCutPhiU_p[iter][iter2]);
      CleanHist(mean_projAEta_p[iter][iter2]);
      CleanHist(mean_projAEtaD_p[iter][iter2]);
      CleanHist(mean_projAEtaU_p[iter][iter2]);
      CleanHist(mean_projAEtaCut_p[iter][iter2]);
      CleanHist(mean_projAEtaCutD_p[iter][iter2]);
      CleanHist(mean_projAEtaCutU_p[iter][iter2]);
      CleanHist(mean_projAPhi_p[iter][iter2]);
      CleanHist(mean_projAPhiD_p[iter][iter2]);
      CleanHist(mean_projAPhiU_p[iter][iter2]);
      CleanHist(mean_projAPhiCut_p[iter][iter2]);
      CleanHist(mean_projAPhiCutD_p[iter][iter2]);
      CleanHist(mean_projAPhiCutU_p[iter][iter2]);

      CleanHist(mean_projARFOR_p[iter][iter2]);
      CleanHist(mean_projARFORD_p[iter][iter2]);
      CleanHist(mean_projARFORU_p[iter][iter2]);
      CleanHist(mean_projARFORMID_p[iter][iter2]);
      CleanHist(mean_projARFORMIDD_p[iter][iter2]);
      CleanHist(mean_projARFORMIDU_p[iter][iter2]);
      CleanHist(mean_projARFORFOR_p[iter][iter2]);
      CleanHist(mean_projARFORFORD_p[iter][iter2]);
      CleanHist(mean_projARFORFORU_p[iter][iter2]);
      CleanHist(mean_projARMID_p[iter][iter2]);
      CleanHist(mean_projARMIDD_p[iter][iter2]);
      CleanHist(mean_projARMIDU_p[iter][iter2]);

      CleanHist(mean_multAR_p[iter][iter2]);
      CleanHist(mean_multARD_p[iter][iter2]);
      CleanHist(mean_multARU_p[iter][iter2]);
      CleanHist(mean_multARCut_p[iter][iter2]);
      CleanHist(mean_multARCutD_p[iter][iter2]);
      CleanHist(mean_multARCutU_p[iter][iter2]);
      CleanHist(mean_multARCutEta_p[iter][iter2]);
      CleanHist(mean_multARCutEtaD_p[iter][iter2]);
      CleanHist(mean_multARCutEtaU_p[iter][iter2]);
      CleanHist(mean_multARCutPhi_p[iter][iter2]);
      CleanHist(mean_multARCutPhiD_p[iter][iter2]);
      CleanHist(mean_multARCutPhiU_p[iter][iter2]);
      CleanHist(mean_multAEta_p[iter][iter2]);
      CleanHist(mean_multAEtaD_p[iter][iter2]);
      CleanHist(mean_multAEtaU_p[iter][iter2]);
      CleanHist(mean_multAEtaCut_p[iter][iter2]);
      CleanHist(mean_multAEtaCutD_p[iter][iter2]);
      CleanHist(mean_multAEtaCutU_p[iter][iter2]);
      CleanHist(mean_multAPhi_p[iter][iter2]);
      CleanHist(mean_multAPhiD_p[iter][iter2]);
      CleanHist(mean_multAPhiU_p[iter][iter2]);
      CleanHist(mean_multAPhiCut_p[iter][iter2]);
      CleanHist(mean_multAPhiCutD_p[iter][iter2]);
      CleanHist(mean_multAPhiCutU_p[iter][iter2]);
    }
  }
  return;
}


int makeDiJetHists(const std::string inName, sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  TH1::SetDefaultSumw2();
  Bool_t montecarlo = isMonteCarlo(sType);

  setFileTag(inName);

  inFile_p = new TFile(inName.c_str(), "READ");
  GetDiJetAnaSkim(inFile_p, sType);

  std::cout << "AnaSkim Loaded" << std::endl;

  if(montecarlo)
    trackTreeAna_p->AddFriend(genTreeAna_p);

  jetTreeAna_p->AddFriend(trackTreeAna_p);

  std::string outName = inName;
  const std::string cutString[2] = {"AnaSkim", ".root"};
  const std::string repString[2] = {"Hist", ""};

  std::cout << "Replace string" << std::endl;

  for(Int_t iter = 0; iter < 2; iter++){
    std::size_t strIndex = outName.find(cutString[iter]);
    if(!(strIndex == std::string::npos)){
      outName.replace(strIndex, cutString[iter].length(), repString[iter]);
    }
  }

  const std::string Corr[2] = {"", "Corr"};
  const std::string Tight[2] = {"", "Tight"};

  for(Int_t corrIter = 1; corrIter < 2; corrIter++){
    for(Int_t tightIter = 0; tightIter < 1; tightIter++){
      for(Int_t algIter = 4; algIter < 8; algIter++){

	if(algIter == 4  || algIter == 7){	
	  makeImbAHist(jetTreeAna_p, outName, algIter, 0, 19, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	  
	  if(isHI(sType)){
	    makeImbAHist(jetTreeAna_p, outName, algIter, 20, 59, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	    makeImbAHist(jetTreeAna_p, outName, algIter, 60, 99, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	    makeImbAHist(jetTreeAna_p, outName, algIter, 100, 199, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	  }	
	}
	
	/*	
	makeImbACNCHist(jetTreeAna_p, outName, algIter, 0, 59, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	if(isHI(sType))	makeImbACNCHist(jetTreeAna_p, outName, algIter, 60, 199, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);      
	*/
	
	std::cout << algIter << std::endl;

	if(algIter == 4 || algIter == 7){
	  makeImbARHist(jetTreeAna_p, outName, "r", algIter, 0, 59, Corr[corrIter], sType, isHighPtTrk);
	  if(isHI(sType)) makeImbARHist(jetTreeAna_p, outName, "r", algIter, 60, 199, Corr[corrIter], sType, isHighPtTrk);
	}
	
	if(algIter == 7 && montecarlo){
	  makeImbARHist(jetTreeAna_p, outName, "g", algIter, 0, 59, Corr[corrIter], sType, isHighPtTrk);
	  if(isHI(sType)) makeImbARHist(jetTreeAna_p, outName, "g", algIter, 60, 199, Corr[corrIter], sType, isHighPtTrk);
	}
      }
    }
  }

  return 0;
}
