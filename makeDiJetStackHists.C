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

#include <fstream>

TFile* inFile_p = 0;
TFile* outFile_p = 0;

enum sampleType{
  kHIDATA, //0
  kHIMC,   //1                                                                                       
  kPPDATA, //2                                                                                      
  kPPMC,   //3                                                                                      
  kPADATA, //4                                                                                      
  kPAMC    //5                                                                                           
};

const std::string algType[10] = {"Pu3Calo", "Pu4Calo", "Pu5Calo", "Vs2Calo", "Vs3Calo", "Vs4Calo", "Vs5Calo", "T", "PuPF", "VsPF"};

Bool_t isMonteCarlo(sampleType sType = kHIDATA){
  if(sType == kHIMC || sType == kPPMC || sType == kPAMC) return true;
  else return false;
}


Bool_t isHI(sampleType sType = kHIDATA){
  if(sType == kHIDATA || sType == kHIMC) return true;
  else return false;
}

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

void makeImbAHist(TFile* inF_p, const std::string outName, const std::string gR, Int_t setNum, Int_t centLow, Int_t centHi, const std::string Corr = "", const std::string Tight = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);

  Int_t setCorrNum = setNum;
  if(!strcmp("Corr", Corr.c_str()))
    setCorrNum = setNum + 8;

  const Int_t nBins = retBinNumber(Tight);
  Float_t xArr[nBins+1];
  Float_t xArrCut[nBins+1];

  if(nBins == 4) getBinArr(nBins, xArr, xArrCut, loose);
  else if(nBins == 8) getBinArr(nBins, xArr, xArrCut, tight);

  const std::string centString = getCentString(sType, centLow, centHi);

  TH1F* rImbAHist_p[6];
  TH1F* mean_rProjA_p[6][nBins];
  
  InitHist(rImbAHist_p);
  BookHist(rImbAHist_p, gR, algType[setNum], "Proj", "", Corr, centString, nBins, xArr, 0.50, niceNumCNC);


  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      mean_rProjA_p[iter][iter2] = (TH1F*)inF_p->Get(Form("mean_%s%sProjA%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2));
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
  
  for(Int_t iter = 0; iter < 6; iter++){
    delete rImbAHist_p[iter];
    rImbAHist_p[iter] = 0;

    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      delete mean_rProjA_p[iter][iter2];
      mean_rProjA_p[iter][iter2] = 0;
    }
  }

}

/*
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
*/


void makeImbARHist(TFile* inF_p, const std::string outName, const std::string gR, Int_t setNum, const std::string rStr, Int_t centLow, Int_t centHi, const std::string Corr = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  const Int_t nBins = 10;
  Float_t rBins[11] = {0.0001, 0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 1.9999};
  const std::string centString = getCentString(sType, centLow, centHi);

  const Int_t nBins2 = 9;
  Float_t rBins2[10] = {0.0001, 0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.7999};

  TH1F* imbARHist_p[6];
  TH1F* mean_projAR_p[6][nBins];
  InitHist(imbARHist_p);
  BookHist(imbARHist_p, gR, algType[setNum], "Proj", rStr.c_str(), Corr, centString, nBins, rBins, 2.00, niceNumR);

  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      mean_projAR_p[iter][iter2] = (TH1F*)inF_p->Get(Form("mean_%s%sProj%s%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), rStr.c_str(), FPT[iter].c_str(), centString.c_str(), iter2));
    }
  }

  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      if(mean_projAR_p[iter][iter2]->GetEntries() != 0){
	imbARHist_p[iter]->SetBinContent(iter2+1, mean_projAR_p[iter][iter2]->GetMean());
	imbARHist_p[iter]->SetBinError(iter2+1, mean_projAR_p[iter][iter2]->GetMeanError());
      }
    }
  }
  outFile_p = new TFile(Form("%s.root", outName.c_str()), "UPDATE");
  std::cout << outName << std::endl;

  for(Int_t iter = 0; iter < 6; iter++){
    imbARHist_p[iter]->Write("", TObject::kOverwrite);
  } 

  outFile_p->Close();
  delete outFile_p;
  outFile_p = 0;

  for(Int_t iter = 0; iter < 6; iter++){
    CleanHist(imbARHist_p[iter]);
  }
  return;
}

int makeDiJetStackHists(const std::string inName, sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  TH1::SetDefaultSumw2();
  Bool_t montecarlo = isMonteCarlo(sType);

  setFileTag(inName);

  inFile_p = new TFile(inName.c_str(), "READ");
 
  std::string outName = inName;
  const std::string cutString[2] = {"MeanHists", ".root"};
  const std::string repString[2] = {"StackHists", ""};

  std::cout << "Replace string" << std::endl;

  for(Int_t iter = 0; iter < 2; iter++){
    std::size_t strIndex = outName.find(cutString[iter]);
    if(!(strIndex == std::string::npos)){
      outName.replace(strIndex, cutString[iter].length(), repString[iter]);
    }
  }

  const std::string Corr[2] = {"", "Corr"};
  const std::string Tight[2] = {"", "Tight"};

  const char* rStr[39] = {"R", "RD", "RU", "Eta", "EtaD", "EtaU", "Phi", "PhiD", "PhiU", "RCut", "RCutD", "RCutU", "RCutEta", "RCutEtaD", "RCutEtaU", "RCutPhi", "RCutPhiD", "RCutPhiU", "EtaCut", "EtaCutD", "EtaCutU", "PhiCut", "PhiCutD", "PhiCutU", "RFOR", "RFORD", "RFORU", "RFORMID", "RFORMIDD", "RFORMIDU", "RFORFOR", "RFORFORD", "RFORFORU", "RMID", "RMIDD", "RMIDU", "R2", "R2U", "R2D"};


  for(Int_t corrIter = 1; corrIter < 2; corrIter++){
    for(Int_t tightIter = 0; tightIter < 1; tightIter++){
      for(Int_t algIter = 4; algIter < 8; algIter++){

	if(algIter == 4){	
	  makeImbAHist(inFile_p, outName, "r", algIter, 0, 19, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	  
	  if(isHI(sType)){
	    makeImbAHist(inFile_p, outName, "r", algIter, 20, 59, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	    makeImbAHist(inFile_p, outName, "r", algIter, 60, 99, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	    makeImbAHist(inFile_p, outName, "r", algIter, 100, 199, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	  }	
	}

	if(algIter == 4 || algIter == 7){	
	  makeImbAHist(inFile_p, outName, "g", algIter, 0, 19, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	  
	  if(isHI(sType)){
	    makeImbAHist(inFile_p, outName, "g", algIter, 20, 59, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	    makeImbAHist(inFile_p, outName, "g", algIter, 60, 99, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	    makeImbAHist(inFile_p, outName, "g", algIter, 100, 199, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	  }	
	}
	
	/*	
	makeImbACNCHist(inFile_p, outName, algIter, 0, 59, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	if(isHI(sType))	makeImbACNCHist(inFile_p, outName, algIter, 60, 199, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);      
	*/
	
	std::cout << algIter << std::endl;


	for(Int_t rIter = 0; rIter < 36; rIter++){

	  if(algIter == 4){
	    makeImbARHist(inFile_p, outName, "r", algIter, rStr[rIter], 0, 59, Corr[corrIter], sType, isHighPtTrk);
	    if(isHI(sType)) makeImbARHist(inFile_p, outName, "r", algIter, rStr[rIter], 60, 199, Corr[corrIter], sType, isHighPtTrk);
	  }
	
	  if((algIter == 4 || algIter == 7) && montecarlo){
	    makeImbARHist(inFile_p, outName, "g", algIter, rStr[rIter], 0, 59, Corr[corrIter], sType, isHighPtTrk);
	    if(isHI(sType)) makeImbARHist(inFile_p, outName, "g", algIter, rStr[rIter], 60, 199, Corr[corrIter], sType, isHighPtTrk);
	  }
	}


      }
    }
  }

  return 0;
}
