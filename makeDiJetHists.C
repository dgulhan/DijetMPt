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

const Int_t evtCutPos[7] = {2, 2, 2, 6, 6, 6, 6};
const Bool_t inVenn = false;
const Bool_t outVenn = false;

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
}


void BookHist(TH1F* inHist_p[6], const std::string CNCR, const Int_t nBins, Float_t xArr[], Float_t upBound, const Float_t niceNum[4])
{
  for(Int_t iter = 0; iter < 6; iter++){
    inHist_p[iter] = new TH1F(Form("rImbA%sHist_%d_p", CNCR.c_str(), iter), Form("rImbA%sHist_%d_p", CNCR.c_str(), iter), nBins, xArr);

    inHist_p[iter]->GetXaxis()->SetLimits(0.00, upBound);
    niceTH1(inHist_p[iter], niceNum[0], niceNum[1], niceNum[2], niceNum[3]);
      
  }
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

  std::string title[6];
  TH1F* rImbAHist_p[6];
  TH1F* mean_rProjA_p[6][nBins];
  
  InitHist(rImbAHist_p);
  BookHist(rImbAHist_p, "", nBins, xArr, 0.50, niceNumCNC);

  for(Int_t iter = 0; iter < 6; iter++){
    title[iter] = Form("r%sImbProjA%s%s%s_%s_%s_h", algType[setNum].c_str(), FPT[iter].c_str(), Corr.c_str(), Tight.c_str(), centString.c_str(), fileTag.c_str());
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
    rImbAHist_p[iter]->Write(title[iter].c_str());
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

  std::string titleC[6];
  std::string titleNC[6];
  std::string titleC0[6];
  std::string titleC1[6];
  std::string titleC2[6];
  std::string titleC3[6];

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

  BookHist(rImbACHist_p, "C", nBins, xArr, 0.50, niceNumCNC);
  BookHist(rImbANCHist_p, "NC", nBins, xArr, 0.50, niceNumCNC);
  BookHist(rImbAC0Hist_p, "C0", nBins, xArr, 0.50, niceNumCNC);
  BookHist(rImbAC1Hist_p, "C1", nBins, xArr, 0.50, niceNumCNC);
  BookHist(rImbAC2Hist_p, "C2", nBins, xArr, 0.50, niceNumCNC);
  BookHist(rImbAC3Hist_p, "C3", nBins, xArr, 0.50, niceNumCNC);

  for(Int_t iter = 0; iter < 6; iter++){
    titleC[iter] = Form("r%sImbProjAC%s%s%s_%s_%s_h", algType[setNum].c_str(), FPT[iter].c_str(), Corr.c_str(), Tight.c_str(), centString.c_str(), fileTag.c_str());
    titleNC[iter] = Form("r%sImbProjANC%s%s%s_%s_%s_h", algType[setNum].c_str(), FPT[iter].c_str(), Corr.c_str(), Tight.c_str(), centString.c_str(), fileTag.c_str());
    titleC0[iter] = Form("r%sImbProjAC0%s%s%s_%s_%s_h", algType[setNum].c_str(), FPT[iter].c_str(), Corr.c_str(), Tight.c_str(), centString.c_str(), fileTag.c_str());
    titleC1[iter] = Form("r%sImbProjAC1%s%s%s_%s_%s_h", algType[setNum].c_str(), FPT[iter].c_str(), Corr.c_str(), Tight.c_str(), centString.c_str(), fileTag.c_str());
    titleC2[iter] = Form("r%sImbProjAC2%s%s%s_%s_%s_h", algType[setNum].c_str(), FPT[iter].c_str(), Corr.c_str(), Tight.c_str(), centString.c_str(), fileTag.c_str());
    titleC3[iter] = Form("r%sImbProjAC3%s%s%s_%s_%s_h", algType[setNum].c_str(), FPT[iter].c_str(), Corr.c_str(), Tight.c_str(), centString.c_str(), fileTag.c_str());


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

	  if(TMath::Abs(AlgJtEta_[setNum][0]) < 0.5 && TMath::Abs(AlgJtEta_[setNum][1]) < 0.5){
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
    rImbACHist_p[iter]->Write(titleC[iter].c_str());
    rImbANCHist_p[iter]->Write(titleNC[iter].c_str());

    rImbAC0Hist_p[iter]->Write(titleC0[iter].c_str());
    rImbAC1Hist_p[iter]->Write(titleC1[iter].c_str());
    rImbAC2Hist_p[iter]->Write(titleC2[iter].c_str());
    rImbAC3Hist_p[iter]->Write(titleC3[iter].c_str());
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



void makeImbARHist(TTree* anaTree_p, const std::string outName, Int_t setNum, Int_t centLow, Int_t centHi, const std::string Corr = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);                                                    

  Int_t setCorrNum = setNum;
  if(!strcmp("Corr", Corr.c_str()))
    setCorrNum = setNum + 8;

  Int_t evtsPass = 0;
  Int_t evtsPassU = 0;
  Int_t evtsPassD = 0;

  const Int_t nBins = 10;
  Float_t rBins[11] = {0.0001, 0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 1.999};

  const std::string centString = getCentString(sType, centLow, centHi);

  std::string titleR[6];
  TH1F* rImbARHist_p[6];
  TH1F* mean_rProjAR_p[6][nBins];
  std::string titleRD[6];
  TH1F* rImbARDHist_p[6];
  TH1F* mean_rProjARD_p[6][nBins];
  std::string titleRU[6];
  TH1F* rImbARUHist_p[6];
  TH1F* mean_rProjARU_p[6][nBins];


  InitHist(rImbARHist_p);
  BookHist(rImbARHist_p, "R", nBins, rBins, 2.00, niceNumR);
  InitHist(rImbARDHist_p);
  BookHist(rImbARDHist_p, "RD", nBins, rBins, 2.00, niceNumR);
  InitHist(rImbARUHist_p);
  BookHist(rImbARUHist_p, "RU", nBins, rBins, 2.00, niceNumR);


  for(Int_t iter = 0; iter < 6; iter++){
    titleR[iter] = Form("r%sImbProjAR%s%s_%s_%s_h", algType[setNum].c_str(), FPT[iter].c_str(), Corr.c_str(), centString.c_str(), fileTag.c_str());
    titleRD[iter] = Form("r%sImbProjARD%s%s_%s_%s_h", algType[setNum].c_str(), FPT[iter].c_str(), Corr.c_str(), centString.c_str(), fileTag.c_str());
    titleRU[iter] = Form("r%sImbProjARU%s%s_%s_%s_h", algType[setNum].c_str(), FPT[iter].c_str(), Corr.c_str(), centString.c_str(), fileTag.c_str());

    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      mean_rProjAR_p[iter][iter2] = new TH1F(Form("mean_rProjAR_%d_%d_h", iter, iter2), Form("mean_rProjAR_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_rProjARD_p[iter][iter2] = new TH1F(Form("mean_rProjARD_%d_%d_h", iter, iter2), Form("mean_rProjARD_%d_%d_h", iter, iter2), 100, -2000, 2000);
      mean_rProjARU_p[iter][iter2] = new TH1F(Form("mean_rProjARU_%d_%d_h", iter, iter2), Form("mean_rProjARU_%d_%d_h", iter, iter2), 100, -2000, 2000);
    }
  }

  for(Int_t jEntry = 0; jEntry < (Int_t)anaTree_p->GetEntries(); jEntry++){
    anaTree_p->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(isEventCut(setNum, sType, centLow, centHi, isHighPtTrk)) continue;

    if(TMath::Abs(AlgJtEta_[setNum][0]) > 0.5 || TMath::Abs(AlgJtEta_[setNum][1]) > 0.5) continue;

    evtsPass++;

    if(AlgJtAsymm_[setNum] < 0.22) evtsPassD++;
    else evtsPassU++;

    for(Int_t iter = 0; iter < 6; iter++){
      for(Int_t iter2 = 0; iter2 < nBins; iter2++){

	Float_t weight = 1.0;
	if(montecarlo){
	  weight = pthatWeight_;
	  if(isHI(sType)) weight *= centWeight_[setNum];
	}

	mean_rProjAR_p[iter][iter2]->Fill(rAlgImbProjAR_[setCorrNum][iter][iter2], weight);

	if(AlgJtAsymm_[setNum] < 0.22) mean_rProjARD_p[iter][iter2]->Fill(rAlgImbProjAR_[setCorrNum][iter][iter2], weight);
	else mean_rProjARU_p[iter][iter2]->Fill(rAlgImbProjAR_[setCorrNum][iter][iter2], weight);
	
      }
    }
  }

  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      if(mean_rProjAR_p[iter][iter2]->GetEntries() != 0){
	rImbARHist_p[iter]->SetBinContent(iter2+1, mean_rProjAR_p[iter][iter2]->GetMean());
	rImbARHist_p[iter]->SetBinError(iter2+1, mean_rProjAR_p[iter][iter2]->GetMeanError());
      }

      if(mean_rProjARD_p[iter][iter2]->GetEntries() != 0){
	rImbARDHist_p[iter]->SetBinContent(iter2+1, mean_rProjARD_p[iter][iter2]->GetMean());
	rImbARDHist_p[iter]->SetBinError(iter2+1, mean_rProjARD_p[iter][iter2]->GetMeanError());
      }

      if(mean_rProjARU_p[iter][iter2]->GetEntries() != 0){
	rImbARUHist_p[iter]->SetBinContent(iter2+1, mean_rProjARU_p[iter][iter2]->GetMean());
	rImbARUHist_p[iter]->SetBinError(iter2+1, mean_rProjARU_p[iter][iter2]->GetMeanError());
      }
    }
  }

  outFile_p = new TFile(Form("%s.root", outName.c_str()), "UPDATE");
  std::cout << outName << std::endl;

  for(Int_t iter = 0; iter < 6; iter++){
    rImbARHist_p[iter]->Write(titleR[iter].c_str());
    rImbARDHist_p[iter]->Write(titleRD[iter].c_str());
    rImbARUHist_p[iter]->Write(titleRU[iter].c_str());
  } 

  outFile_p->Close();
  delete outFile_p;
  outFile_p = 0;

  std::ofstream txtFile;
  txtFile.open(Form("%s.txt", outName.c_str()), std::ofstream::out | std::ofstream::app);
  if(centLow == 0) txtFile << Form("ProjAR (D) (U) %s,\n", algType[setNum].c_str());
  txtFile << Form("   %d-%d: %d (%d) (%d).\n", centLow, centHi, evtsPass, evtsPassD, evtsPassU);
  txtFile.close();

  for(Int_t iter = 0; iter < 6; iter++){
    delete rImbARHist_p[iter];
    rImbARHist_p[iter] = 0;
    delete rImbARDHist_p[iter];
    rImbARDHist_p[iter] = 0;
    delete rImbARUHist_p[iter];
    rImbARUHist_p[iter] = 0;

    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      delete mean_rProjAR_p[iter][iter2];
      mean_rProjAR_p[iter][iter2] = 0;

      delete mean_rProjARD_p[iter][iter2];
      mean_rProjARD_p[iter][iter2] = 0;

      delete mean_rProjARU_p[iter][iter2];
      mean_rProjARU_p[iter][iter2] = 0;
    }
  }
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
      for(Int_t algIter = 0; algIter < 7; algIter++){

	makeImbAHist(jetTreeAna_p, outName, algIter, 0, 19, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);

	if(isHI(sType)){
	  makeImbAHist(jetTreeAna_p, outName, algIter, 20, 59, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	  makeImbAHist(jetTreeAna_p, outName, algIter, 60, 99, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	  makeImbAHist(jetTreeAna_p, outName, algIter, 100, 199, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	}	

	
	makeImbACNCHist(jetTreeAna_p, outName, algIter, 0, 59, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	if(isHI(sType))	makeImbACNCHist(jetTreeAna_p, outName, algIter, 60, 199, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);      
	
	
	makeImbARHist(jetTreeAna_p, outName, algIter, 0, 59, Corr[corrIter], sType, isHighPtTrk);
	if(isHI(sType)) makeImbARHist(jetTreeAna_p, outName, algIter, 60, 199, Corr[corrIter], sType, isHighPtTrk);
      }
    }
  }

  return 0;
}
