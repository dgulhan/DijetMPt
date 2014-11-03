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

void makeImbAMeanHist(TTree* anaTree_p, const std::string outName, const std::string gR, Int_t setNum, Int_t centLow, Int_t centHi, const std::string Corr = "", const std::string Tight = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
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

  TH1F* mean_rProjA_p[6][nBins];
  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      mean_rProjA_p[iter][iter2] = new TH1F(Form("mean_%s%sProjA%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_rProjA%s_%s_%d_h", FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
    }
  }    

  for(Int_t jEntry = 0; jEntry < (Int_t)anaTree_p->GetEntries(); jEntry++){
    anaTree_p->GetEntry(jEntry);

    if(jEntry%1000 == 0) std::cout << jEntry << std::endl;

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
  std::cout << outName << std::endl;
  outFile_p = new TFile(Form("%s.root", outName.c_str()), "UPDATE");

  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      mean_rProjA_p[iter][iter2]->Write("", TObject::kOverwrite);
    }
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
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      CleanHist(mean_rProjA_p[iter][iter2]);
    }
  }
  return;
}


void makeImbACNCMeanHist(TTree* anaTree_p, const std::string outName, Int_t setNum, Int_t centLow, Int_t centHi, const std::string Corr = "", const std::string Tight = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
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

  TH1F* mean_rProjAC_p[6][nBins];
  TH1F* mean_rProjANC_p[6][nBins];
  TH1F* mean_rProjAC0_p[6][nBins];
  TH1F* mean_rProjAC1_p[6][nBins];
  TH1F* mean_rProjAC2_p[6][nBins];
  TH1F* mean_rProjAC3_p[6][nBins];

  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      mean_rProjAC_p[iter][iter2] = new TH1F(Form("mean_rProjAC%s_%d_h", FPT[iter].c_str(), iter2), Form("mean_rProjAC%s_%d_h", FPT[iter].c_str(), iter2), 4000, -2000, 2000);
      mean_rProjANC_p[iter][iter2] = new TH1F(Form("mean_rProjANC%s_%d_h", FPT[iter].c_str(), iter2), Form("mean_rProjANC%s_%d_h", FPT[iter].c_str(), iter2), 4000, -2000, 2000);
      mean_rProjAC0_p[iter][iter2] = new TH1F(Form("mean_rProjAC0%s_%d_h", FPT[iter].c_str(), iter2), Form("mean_rProjAC0%s_%d_h", FPT[iter].c_str(), iter2), 4000, -2000, 2000);
      mean_rProjAC1_p[iter][iter2] = new TH1F(Form("mean_rProjAC1%s_%d_h", FPT[iter].c_str(), iter2), Form("mean_rProjAC1%s_%d_h", FPT[iter].c_str(), iter2), 4000, -2000, 2000);
      mean_rProjAC2_p[iter][iter2] = new TH1F(Form("mean_rProjAC2%s_%d_h", FPT[iter].c_str(), iter2), Form("mean_rProjAC2%s_%d_h", FPT[iter].c_str(), iter2), 4000, -2000, 2000);
      mean_rProjAC3_p[iter][iter2] = new TH1F(Form("mean_rProjAC3%s_%d_h", FPT[iter].c_str(), iter2), Form("mean_rProjAC3%s_%d_h", FPT[iter].c_str(), iter2), 4000, -2000, 2000);
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
  
  outFile_p = new TFile(Form("%s.root", outName.c_str()), "UPDATE");
  std::cout << outName << std::endl;

  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      mean_rProjAC_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_rProjANC_p[iter][iter2]->Write("", TObject::kOverwrite);

      mean_rProjAC0_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_rProjAC1_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_rProjAC2_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_rProjAC3_p[iter][iter2]->Write("", TObject::kOverwrite);
    }
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
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      CleanHist(mean_rProjAC_p[iter][iter2]);
      CleanHist(mean_rProjANC_p[iter][iter2]);

      CleanHist(mean_rProjAC0_p[iter][iter2]);
      CleanHist(mean_rProjAC1_p[iter][iter2]);
      CleanHist(mean_rProjAC2_p[iter][iter2]);
      CleanHist(mean_rProjAC3_p[iter][iter2]);
    }
  }

  return;
}



void makeImbARMeanHist(TTree* anaTree_p, const std::string outName, const std::string gR, Int_t setNum, Int_t centLow, Int_t centHi, const std::string Corr = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);                                                    

  Int_t setCorrNum = setNum;
  if(!strcmp("Corr", Corr.c_str())) setCorrNum = setNum + 8;

  Int_t evtsPass = 0;
  Int_t evtsPassU = 0;
  Int_t evtsPassD = 0;

  const Int_t nBins = 10;
  const std::string centString = getCentString(sType, centLow, centHi);

  const Int_t nBins2 = 9;

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

  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      mean_projAR_p[iter][iter2] = new TH1F(Form("mean_%s%sProjR%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjR%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARD_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARU_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);

      if(iter2 < nBins2){
	mean_projAR2_p[iter][iter2] = new TH1F(Form("mean_%s%sProjR2%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjR2%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
	mean_projAR2D_p[iter][iter2] = new TH1F(Form("mean_%s%sProjR2D%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjR2D%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
	mean_projAR2U_p[iter][iter2] = new TH1F(Form("mean_%s%sProjR2U%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjR2U%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      }

      mean_projARCut_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRCut%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRCut%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARCutD_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRCutD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRCutD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARCutU_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRCutU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRCutU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARCutEta_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRCutEta%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRCutEta%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARCutEtaD_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRCutEtaD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRCutEtaD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARCutEtaU_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRCutEtaU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRCutEtaU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARCutPhi_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRCutPhi%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRCutPhi%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARCutPhiD_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRCutPhiD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRCutPhiD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARCutPhiU_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRCutPhiU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRCutPhiU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projAEta_p[iter][iter2] = new TH1F(Form("mean_%s%sProjEta%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjEta%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projAEtaD_p[iter][iter2] = new TH1F(Form("mean_%s%sProjEtaD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjEtaD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projAEtaU_p[iter][iter2] = new TH1F(Form("mean_%s%sProjEtaU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjEtaU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projAEtaCut_p[iter][iter2] = new TH1F(Form("mean_%s%sProjEtaCut%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjEtaCut%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projAEtaCutD_p[iter][iter2] = new TH1F(Form("mean_%s%sProjEtaCutD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjEtaCutD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projAEtaCutU_p[iter][iter2] = new TH1F(Form("mean_%s%sProjEtaCutU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjEtaCutU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projAPhi_p[iter][iter2] = new TH1F(Form("mean_%s%sProjPhi%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjPhi%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projAPhiD_p[iter][iter2] = new TH1F(Form("mean_%s%sProjPhiD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjPhiD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projAPhiU_p[iter][iter2] = new TH1F(Form("mean_%s%sProjPhiU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjPhiU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projAPhiCut_p[iter][iter2] = new TH1F(Form("mean_%s%sProjPhiCut%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjPhiCut%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projAPhiCutD_p[iter][iter2] = new TH1F(Form("mean_%s%sProjPhiCutD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjPhiCutD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projAPhiCutU_p[iter][iter2] = new TH1F(Form("mean_%s%sProjPhiCutU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjPhiCutU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);


      mean_projARFOR_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRFOR%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRFOR%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARFORD_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRFORD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRFORD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARFORU_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRFORU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRFORU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARFORMID_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRFORMID%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRFORMID%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARFORMIDD_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRFORMIDD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRFORMIDD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARFORMIDU_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRFORMIDU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRFORMIDU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARFORFOR_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRFORFOR%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRFORFOR%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARFORFORD_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRFORFORD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRFORFORD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARFORFORU_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRFORFORU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRFORFORU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARMID_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRMID%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRMID%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARMIDD_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRMIDD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRMIDD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_projARMIDU_p[iter][iter2] = new TH1F(Form("mean_%s%sProjRMIDU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sProjRMIDU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);


      mean_multAR_p[iter][iter2] = new TH1F(Form("mean_%s%sMultR%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultR%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multARD_p[iter][iter2] = new TH1F(Form("mean_%s%sMultRD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultRD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multARU_p[iter][iter2] = new TH1F(Form("mean_%s%sMultRU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultRU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multARCut_p[iter][iter2] = new TH1F(Form("mean_%s%sMultRCut%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultRCut%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multARCutD_p[iter][iter2] = new TH1F(Form("mean_%s%sMultRCutD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultRCutD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multARCutU_p[iter][iter2] = new TH1F(Form("mean_%s%sMultRCutU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultRCutU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multARCutEta_p[iter][iter2] = new TH1F(Form("mean_%s%sMultRCutEta%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultRCutEta%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multARCutEtaD_p[iter][iter2] = new TH1F(Form("mean_%s%sMultRCutEtaD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultRCutEtaD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multARCutEtaU_p[iter][iter2] = new TH1F(Form("mean_%s%sMultRCutEtaU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultRCutEtaU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multARCutPhi_p[iter][iter2] = new TH1F(Form("mean_%s%sMultRCutPhi%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultRCutPhi%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multARCutPhiD_p[iter][iter2] = new TH1F(Form("mean_%s%sMultRCutPhiD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultRCutPhiD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multARCutPhiU_p[iter][iter2] = new TH1F(Form("mean_%s%sMultRCutPhiU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultRCutPhiU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multAEta_p[iter][iter2] = new TH1F(Form("mean_%s%sMultEta%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultEta%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multAEtaD_p[iter][iter2] = new TH1F(Form("mean_%s%sMultEtaD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultEtaD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multAEtaU_p[iter][iter2] = new TH1F(Form("mean_%s%sMultEtaU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultEtaU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multAEtaCut_p[iter][iter2] = new TH1F(Form("mean_%s%sMultEtaCut%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultEtaCut%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multAEtaCutD_p[iter][iter2] = new TH1F(Form("mean_%s%sMultEtaCutD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultEtaCutD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multAEtaCutU_p[iter][iter2] = new TH1F(Form("mean_%s%sMultEtaCutU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultEtaCutU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multAPhi_p[iter][iter2] = new TH1F(Form("mean_%s%sMultPhi%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultPhi%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multAPhiD_p[iter][iter2] = new TH1F(Form("mean_%s%sMultPhiD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultPhiD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multAPhiU_p[iter][iter2] = new TH1F(Form("mean_%s%sMultPhiU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultPhiU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multAPhiCut_p[iter][iter2] = new TH1F(Form("mean_%s%sMultPhiCut%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultPhiCut%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multAPhiCutD_p[iter][iter2] = new TH1F(Form("mean_%s%sMultPhiCutD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultPhiCutD%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
      mean_multAPhiCutU_p[iter][iter2] = new TH1F(Form("mean_%s%sMultPhiCutU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), Form("mean_%s%sMultPhiCutU%s_%s_%d_h", gR.c_str(), algType[setNum].c_str(), FPT[iter].c_str(), centString.c_str(), iter2), 4000, -2000, 2000);
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

  outFile_p = new TFile(Form("%s.root", outName.c_str()), "UPDATE");
  std::cout << outName << std::endl;

  for(Int_t iter = 0; iter < 6; iter++){
    for(Int_t iter2 = 0; iter2 < nBins; iter2++){
      mean_projAR_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARU_p[iter][iter2]->Write("", TObject::kOverwrite);
      
      if(iter2 < nBins2){
	mean_projAR2_p[iter][iter2]->Write("", TObject::kOverwrite);
	mean_projAR2D_p[iter][iter2]->Write("", TObject::kOverwrite);
	mean_projAR2U_p[iter][iter2]->Write("", TObject::kOverwrite);
      }
      
      mean_projARCut_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARCutD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARCutU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARCutEta_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARCutEtaD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARCutEtaU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARCutPhi_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARCutPhiD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARCutPhiU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projAEta_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projAEtaD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projAEtaU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projAEtaCut_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projAEtaCutD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projAEtaCutU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projAPhi_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projAPhiD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projAPhiU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projAPhiCut_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projAPhiCutD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projAPhiCutU_p[iter][iter2]->Write("", TObject::kOverwrite);
      
      mean_projARFOR_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARFORD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARFORU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARFORMID_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARFORMIDD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARFORMIDU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARFORFOR_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARFORFORD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARFORFORU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARMID_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARMIDD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_projARMIDU_p[iter][iter2]->Write("", TObject::kOverwrite);
      
      mean_multAR_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multARD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multARU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multARCut_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multARCutD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multARCutU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multARCutEta_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multARCutEtaD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multARCutEtaU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multARCutPhi_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multARCutPhiD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multARCutPhiU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multAEta_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multAEtaD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multAEtaU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multAEtaCut_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multAEtaCutD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multAEtaCutU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multAPhi_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multAPhiD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multAPhiU_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multAPhiCut_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multAPhiCutD_p[iter][iter2]->Write("", TObject::kOverwrite);
      mean_multAPhiCutU_p[iter][iter2]->Write("", TObject::kOverwrite);
    }
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


int makeDiJetMeanHists(const std::string fList = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false, Int_t num = 0)
{
  TH1::SetDefaultSumw2();

  Bool_t montecarlo = isMonteCarlo(sType);

  std::cout << sType << std::endl;
  std::cout << montecarlo << std::endl;

  std::string buffer;
  std::vector<std::string> listOfFiles;
  int nLines = 0;
  ifstream inFile(fList.data());

  std::cout << fList << std::endl;
  std::cout << inFile.is_open() << std::endl;

  if(!inFile.is_open()){
    std::cout << "Error opening file. Exiting." <<std::endl;
    return 1;
  }
  else{
    while(true){
      inFile >> buffer;
      if(inFile.eof()) break;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }

  std::cout << "FileList Loaded" << std::endl;

  for(Int_t iter = 0; iter < (Int_t)(listOfFiles.size()); iter++){
    std::cout << listOfFiles[iter] << std::endl;
  }

  std::cout << "FileJob: " << listOfFiles[num] << std::endl;

  //  setFileTag(inName);
  inFile_p = new TFile(listOfFiles[num].data(), "READ");
  GetDiJetAnaSkim(inFile_p, sType);

  std::cout << "AnaSkim Loaded" << std::endl;

  if(montecarlo)
    trackTreeAna_p->AddFriend(genTreeAna_p);

  jetTreeAna_p->AddFriend(trackTreeAna_p);

  std::string outName = listOfFiles[num];
  const std::string cullString = "/";
  const std::string cutString[2] = {"AnaSkim", ".root"};
  const std::string repString[2] = {"MeanHists", ""};

  std::cout << "Cull string" << std::endl;

  while(true){
    std::size_t strIndex = outName.find(cullString);

    if(strIndex == std::string::npos) break;

    outName.replace(0, strIndex + 1, "");
  }

  std::cout << "Replace string" << std::endl;

  for(Int_t iter = 0; iter < 2; iter++){
    std::size_t strIndex = outName.find(cutString[iter]);
    if(!(strIndex == std::string::npos)){
      outName.replace(strIndex, cutString[iter].length(), repString[iter]);
    }
  }

  std::cout << outName << std::endl;

  const std::string Corr[2] = {"", "Corr"};
  const std::string Tight[2] = {"", "Tight"};

  for(Int_t corrIter = 1; corrIter < 2; corrIter++){
    for(Int_t tightIter = 0; tightIter < 1; tightIter++){
      for(Int_t algIter = 4; algIter < 8; algIter++){
      
	if(algIter == 4){	
	  makeImbAMeanHist(jetTreeAna_p, outName, "r", algIter, 0, 19, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);

	  if(isHI(sType)){
	    makeImbAMeanHist(jetTreeAna_p, outName, "r", algIter, 20, 59, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	    makeImbAMeanHist(jetTreeAna_p, outName, "r", algIter, 60, 99, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	    makeImbAMeanHist(jetTreeAna_p, outName, "r", algIter, 100, 199, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);	    
	  }
	}
	 
	if(montecarlo && (algIter == 4 || algIter == 7)){
	  makeImbAMeanHist(jetTreeAna_p, outName, "g", algIter, 0, 19, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	    
	  if(isHI(sType)){
	    makeImbAMeanHist(jetTreeAna_p, outName, "g", algIter, 20, 59, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	    makeImbAMeanHist(jetTreeAna_p, outName, "g", algIter, 60, 99, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	    makeImbAMeanHist(jetTreeAna_p, outName, "g", algIter, 100, 199, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	  }
	}	
	
	/*	
	makeImbACNCMeanHist(jetTreeAna_p, outName, algIter, 0, 59, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);
	if(isHI(sType))	makeImbACNCMeanHist(jetTreeAna_p, outName, algIter, 60, 199, Corr[corrIter], Tight[tightIter], sType, isHighPtTrk);      
	*/
	  
	if(algIter == 4){
	  makeImbARMeanHist(jetTreeAna_p, outName, "r", algIter, 0, 59, Corr[corrIter], sType, isHighPtTrk);
	  if(isHI(sType)) makeImbARMeanHist(jetTreeAna_p, outName, "r", algIter, 60, 199, Corr[corrIter], sType, isHighPtTrk);
	}
	
	if((algIter == 4 || algIter == 7) && montecarlo){
	  makeImbARMeanHist(jetTreeAna_p, outName, "g", algIter, 0, 59, Corr[corrIter], sType, isHighPtTrk);
	  if(isHI(sType)) makeImbARMeanHist(jetTreeAna_p, outName, "g", algIter, 60, 199, Corr[corrIter], sType, isHighPtTrk);
	}
      }
    }
  }
  
  return 1;
}


int main(int argc, char* argv[])
{
  if(argc != 5){
    std::cout << "Usage: makeDiJetMeanHists <inputFile> <sampleType> <isHITrk> <#>" << std::endl;
    return 1;
  }

  int rStatus = -1;

  rStatus = makeDiJetMeanHists(argv[1], sampleType(atoi(argv[2])), Bool_t(atoi(argv[3])), atoi(argv[4]));

  return rStatus;
}
