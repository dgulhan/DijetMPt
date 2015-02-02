//=============================================                      
// Author: Chris McGinn
//                     
// DiJet Histogram Maker, Missing Pt
//                                                                            
//=============================================     

#include <string>
#include <vector>

#include "TTree.h"
#include "TDatime.h"
#include "TFile.h"
#include "diJetFileTag.h"
#include "/net/hisrv0001/home/cfmcginn/DijetMPt/CMSSW_5_3_20/src/DijetInitialSkim/cfmVectFunc.h"
#include "/net/hisrv0001/home/cfmcginn/DijetMPt/CMSSW_5_3_20/src/DijetHists/getWeightedMean.h"
#include "cfmDijetMeanTree.h"

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


Bool_t isEventCut(Int_t setNum, sampleType sType, Bool_t isHighPtTrk = false)
{
  if(!eventSet_[setNum]) return true;
  if(AlgJtPt_[setNum][0] < leadJtPtCut[setNum] || AlgJtPt_[setNum][1] < subLeadJtPtCut[setNum]) return true;
  if(TMath::Abs(AlgJtEta_[setNum][0]) > 1.6 || TMath::Abs(AlgJtEta_[setNum][1]) > 1.6) return true;
  if(AlgJtDelPhi12_[setNum] < 5.0*TMath::Pi()/6.0) return true;

  //  if(AlgRefPt_[setNum][0] < AlgRefPt_[setNum][1]) return true;

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

void makeImbAMeanTree(TTree* anaTree_p, const std::string outName, sampleType sType = kHIDATA, Bool_t isHighPtTrk = false, Bool_t isCombined = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

  std::vector<float>* rProjAVal_p[2*nSumAlg][nPtBins][nAjCentBins][nAjBins];
  std::vector<float>* rProjAWeight_p[2*nSumAlg][nPtBins][nAjCentBins][nAjBins];

  std::vector<float>* rProjARVal_p[2*nSumAlg][nPtBins][nRCentBins][nRBins];
  std::vector<float>* rProjARWeight_p[2*nSumAlg][nPtBins][nRCentBins][nRBins];

  std::vector<float>* rProjARDVal_p[2*nSumAlg][nPtBins][nRCentBins][nRBins];
  std::vector<float>* rProjARDWeight_p[2*nSumAlg][nPtBins][nRCentBins][nRBins];

  std::vector<float>* rProjARUVal_p[2*nSumAlg][nPtBins][nRCentBins][nRBins];
  std::vector<float>* rProjARUWeight_p[2*nSumAlg][nPtBins][nRCentBins][nRBins];

  for(Int_t algIter = 0; algIter < 2*nSumAlg; algIter++){
    for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
      for(Int_t ajCentIter = 0; ajCentIter < nAjCentBins; ajCentIter++){
	for(Int_t ajIter = 0; ajIter < nAjBins; ajIter++){
	  rProjAVal_p[algIter][ptIter][ajCentIter][ajIter] = new std::vector<float>;
	  rProjAWeight_p[algIter][ptIter][ajCentIter][ajIter] = new std::vector<float>;
	}
      }

      for(Int_t rCentIter = 0; rCentIter < nRCentBins; rCentIter++){
	for(Int_t rIter = 0; rIter < nRBins; rIter++){
	  rProjARVal_p[algIter][ptIter][rCentIter][rIter] = new std::vector<float>;
	  rProjARWeight_p[algIter][ptIter][rCentIter][rIter] = new std::vector<float>;

	  rProjARDVal_p[algIter][ptIter][rCentIter][rIter] = new std::vector<float>;
	  rProjARDWeight_p[algIter][ptIter][rCentIter][rIter] = new std::vector<float>;

	  rProjARUVal_p[algIter][ptIter][rCentIter][rIter] = new std::vector<float>;
	  rProjARUWeight_p[algIter][ptIter][rCentIter][rIter] = new std::vector<float>;
	}
      }
    }
  }


  std::vector<float>* gProjAVal_p[nSumAlg][nPtBins][nAjCentBins][nAjBins];
  std::vector<float>* gProjAWeight_p[nSumAlg][nPtBins][nAjCentBins][nAjBins];

  std::vector<float>* gProjARVal_p[nSumAlg][nPtBins][nRCentBins][nRBins];
  std::vector<float>* gProjARWeight_p[nSumAlg][nPtBins][nRCentBins][nRBins];

  std::vector<float>* gProjARDVal_p[nSumAlg][nPtBins][nRCentBins][nRBins];
  std::vector<float>* gProjARDWeight_p[nSumAlg][nPtBins][nRCentBins][nRBins];

  std::vector<float>* gProjARUVal_p[nSumAlg][nPtBins][nRCentBins][nRBins];
  std::vector<float>* gProjARUWeight_p[nSumAlg][nPtBins][nRCentBins][nRBins];

  if(montecarlo){
    for(Int_t algIter = 0; algIter < nSumAlg; algIter++){
      for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
	for(Int_t ajCentIter = 0; ajCentIter < nAjCentBins; ajCentIter++){
	  for(Int_t ajIter = 0; ajIter < nAjBins; ajIter++){
	    gProjAVal_p[algIter][ptIter][ajCentIter][ajIter] = new std::vector<float>;
	    gProjAWeight_p[algIter][ptIter][ajCentIter][ajIter] = new std::vector<float>;
	  }
	}
	
	for(Int_t rCentIter = 0; rCentIter < nRCentBins; rCentIter++){
	  for(Int_t rIter = 0; rIter < nRBins; rIter++){
	    gProjARVal_p[algIter][ptIter][rCentIter][rIter] = new std::vector<float>;
	    gProjARWeight_p[algIter][ptIter][rCentIter][rIter] = new std::vector<float>;
	    
	    gProjARDVal_p[algIter][ptIter][rCentIter][rIter] = new std::vector<float>;
	    gProjARDWeight_p[algIter][ptIter][rCentIter][rIter] = new std::vector<float>;
	    
	    gProjARUVal_p[algIter][ptIter][rCentIter][rIter] = new std::vector<float>;
	    gProjARUWeight_p[algIter][ptIter][rCentIter][rIter] = new std::vector<float>;
	  }
	}
      }
    } 
  }

  for(Int_t jEntry = 0; jEntry < (Int_t)anaTree_p->GetEntries(); jEntry++){
    anaTree_p->GetEntry(jEntry);

    if(jEntry%1000 == 0) std::cout << jEntry << std::endl;

    Int_t ajCentPos = 0;
    Int_t rCentPos = 0;
    if(hi){
      for(Int_t ajCentIter = 0; ajCentIter < nAjCentBins; ajCentIter++){
	if(hiBin_ < ajCentBins[ajCentIter]){
	  ajCentPos = ajCentIter;
	  break;
	}
      }
      for(Int_t rCentIter = 0; rCentIter < nRCentBins; rCentIter++){
	if(hiBin_ < rCentBins[rCentIter]){
	  rCentPos = rCentIter;
	  break;
	}
      }
    }

    for(Int_t algIter = 0; algIter < nSumAlg; algIter++){
      if(isEventCut(algIter, sType, isHighPtTrk)) continue;
      
      Float_t weight = 1.0;
      if(montecarlo){
	weight = pthatWeight_;
	if(isHI(sType)) weight *= centWeight_[algIter];
      }

      Int_t ajPos = -1;
      for(Int_t ajIter = 0; ajIter < nAjBins; ajIter++){
	if(AlgJtAsymm12_[algIter] < ajBins[ajIter]){
	  ajPos = ajIter;
	  break;
	}
      }
      
      for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
	rProjAVal_p[algIter][ptIter][ajCentPos][ajPos]->push_back(rAlgImbProjA_[algIter][ptIter]);
	rProjAWeight_p[algIter][ptIter][ajCentPos][ajPos]->push_back(weight);

	rProjAVal_p[algIter+nSumAlg][ptIter][ajCentPos][ajPos]->push_back(rAlgImbProjA_[algIter+nSumAlg][ptIter]);
	rProjAWeight_p[algIter+nSumAlg][ptIter][ajCentPos][ajPos]->push_back(weight);

	if(montecarlo){
	  gProjAVal_p[algIter][ptIter][ajCentPos][ajPos]->push_back(gAlgImbProjA_[algIter][ptIter]);
	  gProjAWeight_p[algIter][ptIter][ajCentPos][ajPos]->push_back(weight);
	}
      }

      if(TMath::Abs(AlgJtEta_[algIter][0]) < 0.5 && TMath::Abs(AlgJtEta_[algIter][1]) < 0.5){
	for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
	  for(Int_t rIter = 0; rIter < nRBins; rIter++){
 	    if(isCombined) rProjARVal_p[algIter][ptIter][rCentPos][rIter]->push_back(rAlgImbProjAR_Comb_[algIter][ptIter][rIter]);
	    else rProjARVal_p[algIter][ptIter][rCentPos][rIter]->push_back(rAlgImbProjAR_[algIter][ptIter][rIter]);
	    rProjARWeight_p[algIter][ptIter][rCentPos][rIter]->push_back(weight);
	    
      if(isCombined) rProjARVal_p[algIter+nSumAlg][ptIter][rCentPos][rIter]->push_back(rAlgImbProjAR_Comb_[algIter+nSumAlg][ptIter][rIter]);
	    else rProjARVal_p[algIter+nSumAlg][ptIter][rCentPos][rIter]->push_back(rAlgImbProjAR_[algIter+nSumAlg][ptIter][rIter]);
	    rProjARWeight_p[algIter+nSumAlg][ptIter][rCentPos][rIter]->push_back(weight);

	    if(montecarlo){
 	      if(isCombined) gProjARVal_p[algIter][ptIter][rCentPos][rIter]->push_back(gAlgImbProjAR_Comb_[algIter][ptIter][rIter]);
        else gProjARVal_p[algIter][ptIter][rCentPos][rIter]->push_back(gAlgImbProjAR_[algIter][ptIter][rIter]);
	      gProjARWeight_p[algIter][ptIter][rCentPos][rIter]->push_back(weight);
	    }
	    
	    if(ajPos < 2){
	      if(isCombined) rProjARDVal_p[algIter][ptIter][rCentPos][rIter]->push_back(rAlgImbProjAR_Comb_[algIter][ptIter][rIter]);
	      else rProjARDVal_p[algIter][ptIter][rCentPos][rIter]->push_back(rAlgImbProjAR_[algIter][ptIter][rIter]);
	      rProjARDWeight_p[algIter][ptIter][rCentPos][rIter]->push_back(weight);
	      
	      if(isCombined) rProjARDVal_p[algIter+nSumAlg][ptIter][rCentPos][rIter]->push_back(rAlgImbProjAR_Comb_[algIter+nSumAlg][ptIter][rIter]);
 	      else rProjARDVal_p[algIter+nSumAlg][ptIter][rCentPos][rIter]->push_back(rAlgImbProjAR_[algIter+nSumAlg][ptIter][rIter]);
	      rProjARDWeight_p[algIter+nSumAlg][ptIter][rCentPos][rIter]->push_back(weight);

	      if(montecarlo){
		if(isCombined) gProjARDVal_p[algIter][ptIter][rCentPos][rIter]->push_back(gAlgImbProjAR_Comb_[algIter][ptIter][rIter]);
 		else gProjARDVal_p[algIter][ptIter][rCentPos][rIter]->push_back(gAlgImbProjAR_[algIter][ptIter][rIter]);
		gProjARDWeight_p[algIter][ptIter][rCentPos][rIter]->push_back(weight);
	      }
	    }
	    else{	      
        if(isCombined) rProjARUVal_p[algIter][ptIter][rCentPos][rIter]->push_back(rAlgImbProjAR_Comb_[algIter][ptIter][rIter]);
	      else rProjARUVal_p[algIter][ptIter][rCentPos][rIter]->push_back(rAlgImbProjAR_[algIter][ptIter][rIter]);
	      rProjARUWeight_p[algIter][ptIter][rCentPos][rIter]->push_back(weight);
	      	      
        if(isCombined) rProjARUVal_p[algIter+nSumAlg][ptIter][rCentPos][rIter]->push_back(rAlgImbProjAR_Comb_[algIter+nSumAlg][ptIter][rIter]);
	      else rProjARUVal_p[algIter+nSumAlg][ptIter][rCentPos][rIter]->push_back(rAlgImbProjAR_[algIter+nSumAlg][ptIter][rIter]);
	      rProjARUWeight_p[algIter+nSumAlg][ptIter][rCentPos][rIter]->push_back(weight);

	      if(montecarlo){
		     if(isCombined) gProjARUVal_p[algIter][ptIter][rCentPos][rIter]->push_back(gAlgImbProjAR_Comb_[algIter][ptIter][rIter]);
    		 else gProjARUVal_p[algIter][ptIter][rCentPos][rIter]->push_back(gAlgImbProjAR_[algIter][ptIter][rIter]);
		     gProjARUWeight_p[algIter][ptIter][rCentPos][rIter]->push_back(weight);
	      }
	    }
	  }
	}
      }


    }
  }

  std::cout << outName << std::endl;
  outFile_p = new TFile(Form("%s.root", outName.c_str()), "UPDATE");

  for(Int_t algIter = 0; algIter < 2*nSumAlg; algIter++){
    for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
      for(Int_t ajCentIter = 0; ajCentIter < nAjCentBins; ajCentIter++){
	for(Int_t ajIter = 0; ajIter < nAjBins; ajIter++){
	  getWeightedMean(rProjAVal_p[algIter][ptIter][ajCentIter][ajIter], rProjAWeight_p[algIter][ptIter][ajCentIter][ajIter], rProjA_MULT_[algIter][ptIter][ajCentIter][ajIter], rProjA_MEAN_[algIter][ptIter][ajCentIter][ajIter], rProjA_SIG_[algIter][ptIter][ajCentIter][ajIter], rProjA_WEIGHT_[algIter][ptIter][ajCentIter][ajIter], rProjA_WEIGHT2_[algIter][ptIter][ajCentIter][ajIter]);

	  if(montecarlo && algIter < nSumAlg) getWeightedMean(gProjAVal_p[algIter][ptIter][ajCentIter][ajIter], gProjAWeight_p[algIter][ptIter][ajCentIter][ajIter], gProjA_MULT_[algIter][ptIter][ajCentIter][ajIter], gProjA_MEAN_[algIter][ptIter][ajCentIter][ajIter], gProjA_SIG_[algIter][ptIter][ajCentIter][ajIter], gProjA_WEIGHT_[algIter][ptIter][ajCentIter][ajIter], gProjA_WEIGHT2_[algIter][ptIter][ajCentIter][ajIter]);
	}
      }

      for(Int_t rCentIter = 0; rCentIter < nRCentBins; rCentIter++){
	for(Int_t rIter = 0; rIter < nRBins; rIter++){
	  getWeightedMean(rProjARVal_p[algIter][ptIter][rCentIter][rIter], rProjARWeight_p[algIter][ptIter][rCentIter][rIter], rProjAR_MULT_[algIter][ptIter][rCentIter][rIter], rProjAR_MEAN_[algIter][ptIter][rCentIter][rIter], rProjAR_SIG_[algIter][ptIter][rCentIter][rIter], rProjAR_WEIGHT_[algIter][ptIter][rCentIter][rIter], rProjAR_WEIGHT2_[algIter][ptIter][rCentIter][rIter]);

	  getWeightedMean(rProjARDVal_p[algIter][ptIter][rCentIter][rIter], rProjARDWeight_p[algIter][ptIter][rCentIter][rIter], rProjARD_MULT_[algIter][ptIter][rCentIter][rIter], rProjARD_MEAN_[algIter][ptIter][rCentIter][rIter], rProjARD_SIG_[algIter][ptIter][rCentIter][rIter], rProjARD_WEIGHT_[algIter][ptIter][rCentIter][rIter], rProjARD_WEIGHT2_[algIter][ptIter][rCentIter][rIter]);

	  getWeightedMean(rProjARUVal_p[algIter][ptIter][rCentIter][rIter], rProjARUWeight_p[algIter][ptIter][rCentIter][rIter], rProjARU_MULT_[algIter][ptIter][rCentIter][rIter], rProjARU_MEAN_[algIter][ptIter][rCentIter][rIter], rProjARU_SIG_[algIter][ptIter][rCentIter][rIter], rProjARU_WEIGHT_[algIter][ptIter][rCentIter][rIter], rProjARU_WEIGHT2_[algIter][ptIter][rCentIter][rIter]);

	  if(montecarlo && algIter < nSumAlg){
	    getWeightedMean(gProjARVal_p[algIter][ptIter][rCentIter][rIter], gProjARWeight_p[algIter][ptIter][rCentIter][rIter], gProjAR_MULT_[algIter][ptIter][rCentIter][rIter], gProjAR_MEAN_[algIter][ptIter][rCentIter][rIter], gProjAR_SIG_[algIter][ptIter][rCentIter][rIter], gProjAR_WEIGHT_[algIter][ptIter][rCentIter][rIter], gProjAR_WEIGHT2_[algIter][ptIter][rCentIter][rIter]);

	    getWeightedMean(gProjARDVal_p[algIter][ptIter][rCentIter][rIter], gProjARDWeight_p[algIter][ptIter][rCentIter][rIter], gProjARD_MULT_[algIter][ptIter][rCentIter][rIter], gProjARD_MEAN_[algIter][ptIter][rCentIter][rIter], gProjARD_SIG_[algIter][ptIter][rCentIter][rIter], gProjARD_WEIGHT_[algIter][ptIter][rCentIter][rIter], gProjARD_WEIGHT2_[algIter][ptIter][rCentIter][rIter]);

	    getWeightedMean(gProjARUVal_p[algIter][ptIter][rCentIter][rIter], gProjARUWeight_p[algIter][ptIter][rCentIter][rIter], gProjARU_MULT_[algIter][ptIter][rCentIter][rIter], gProjARU_MEAN_[algIter][ptIter][rCentIter][rIter], gProjARU_SIG_[algIter][ptIter][rCentIter][rIter], gProjARU_WEIGHT_[algIter][ptIter][rCentIter][rIter], gProjARU_WEIGHT2_[algIter][ptIter][rCentIter][rIter]);
	  }

	}
      }
    }
  }

  InitDiJetMeanTree(sType);
  std::size_t pos = outName.find_last_of("_");
  fileName_ = outName.substr(pos-2);

  trackTreeMean_p->Fill();
  if(montecarlo) genTreeMean_p->Fill();

  trackTreeMean_p->Write();
  if(montecarlo) genTreeMean_p->Write();

  outFile_p->Close();
  delete outFile_p;
  outFile_p = 0;
  
  for(Int_t algIter = 0; algIter < 2*nSumAlg; algIter++){
    for(Int_t ptIter = 0; ptIter < nPtBins; ptIter++){
      for(Int_t ajCentIter = 0; ajCentIter < nAjCentBins; ajCentIter++){
	for(Int_t ajIter = 0; ajIter < nAjBins; ajIter++){

	  delete rProjAVal_p[algIter][ptIter][ajCentIter][ajIter];
	  delete rProjAWeight_p[algIter][ptIter][ajCentIter][ajIter];

	  if(montecarlo && algIter < nSumAlg){
	    delete gProjAVal_p[algIter][ptIter][ajCentIter][ajIter];
	    delete gProjAWeight_p[algIter][ptIter][ajCentIter][ajIter];
	  }

	}
      }

      for(Int_t rCentIter = 0; rCentIter < nRCentBins; rCentIter++){
	for(Int_t rIter = 0; rIter < nRBins; rIter++){

	  delete rProjARVal_p[algIter][ptIter][rCentIter][rIter];
	  delete rProjARWeight_p[algIter][ptIter][rCentIter][rIter];

	  delete rProjARDVal_p[algIter][ptIter][rCentIter][rIter];
	  delete rProjARDWeight_p[algIter][ptIter][rCentIter][rIter];

	  delete rProjARUVal_p[algIter][ptIter][rCentIter][rIter];
	  delete rProjARUWeight_p[algIter][ptIter][rCentIter][rIter];

	  if(montecarlo && algIter < nSumAlg){
	    delete gProjARVal_p[algIter][ptIter][rCentIter][rIter];
	    delete gProjARWeight_p[algIter][ptIter][rCentIter][rIter];

	    delete gProjARDVal_p[algIter][ptIter][rCentIter][rIter];
	    delete gProjARDWeight_p[algIter][ptIter][rCentIter][rIter];

	    delete gProjARUVal_p[algIter][ptIter][rCentIter][rIter];
	    delete gProjARUWeight_p[algIter][ptIter][rCentIter][rIter];
	  }

	}
      }
    }
  }

  return;
}

int makeDiJetMeanTree(const std::string fList = "", sampleType sType = kHIDATA, Bool_t isHighPtTrk = false, Int_t num = 0)
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
  const std::string repString[2] = {"MeanTree", ""};

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

  makeImbAMeanTree(jetTreeAna_p, outName, sType, isHighPtTrk);
  
  return 1;
}


int main(int argc, char* argv[])
{
  if(argc != 5){
    std::cout << "Usage: makeDiJetMeanTree <inputFile> <sampleType> <isHITrk> <#>" << std::endl;
    return 1;
  }

  int rStatus = -1;

  rStatus = makeDiJetMeanTree(argv[1], sampleType(atoi(argv[2])), Bool_t(atoi(argv[3])), atoi(argv[4]));

  return rStatus;
}
