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
#include "TH1F.h"

TFile* inFile_p = 0;
TFile* outFile_p = 0;

const Float_t leadJtCut = 120.;
const Float_t subLeadJtCut = 50.;
const Float_t trkMaxCut = 8.;

std::string getCentString(sampleType sType, Int_t centLow, Int_t centHi)
{
  if(isHI(sType)) return Form("%d%d", (Int_t)(centLow*.5), (Int_t)((centHi+1)*.5));
  else return "PP";
}


Bool_t isEventCut(Int_t setNum, sampleType sType, Bool_t isHighPtTrk = false)
{
  if(!eventSet_[setNum]) return true;

  if(AlgJtPt_[setNum][0] < leadJtCut || AlgJtPt_[setNum][1] < subLeadJtCut) return true;

  if(TMath::Abs(AlgJtEta_[setNum][0]) > 2.0 || TMath::Abs(AlgJtEta_[setNum][1]) > 2.0) return true;

  if(isHighPtTrk && (AlgJtTrkMax_[setNum][0] < trkMaxCut || AlgJtTrkMax_[setNum][1] < trkMaxCut)) return true;

  if(isMonteCarlo(sType) && pthat_ < 80) return true;

  return false;
}


void makeDijetHists_Jet(TTree* anaTree_p, const std::string outName, Int_t setNum, sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

  Int_t centMax = 5;
  if(!hi){
    centMax = 1;
    hiBin_ = 0;
  }

  const Int_t nAJBins = 20;
  const Int_t AJLow = 0.0;
  const Int_t AJHi = 1.0;

  const Int_t nDelPhiBins = 20;
  const Float_t DelPhiLow = 0.0001;
  const Float_t DelPhiHi = TMath::Pi();

  const Int_t nLeadPtBins = 20;
  const Int_t LeadPtLow = 119.5;
  const Int_t LeadPtHi = 699.5;

  const Int_t nSubleadPtBins = 20;
  const Int_t SubleadPtLow = 49.5;
  const Int_t SubleadPtHi = 699.5;

  const Int_t nThirdPtBins = 20;
  const Int_t ThirdPtLow = 29.5;
  const Int_t ThirdPtHi = 399.5;

  const Int_t nEtaBins = 20;
  const Float_t EtaLow = -2.0;
  const Float_t EtaHi = 2.0;

  const Int_t centLow[5] = {0, 20, 40, 60, 100};
  const Int_t centHi[5] = {19, 39, 59, 99, 199};

  std::string centString[5];
  TH1F* ajHist_p[5];
  TH1F* delPhi12Hist_p[5];
  TH1F* delPhi13Hist_p[5];
  TH1F* delPhi23Hist_p[5];
  TH1F* leadPtHist_p[5];
  TH1F* subleadPtHist_p[5];
  TH1F* thirdPtHist_p[5];
  TH1F* leadEtaHist_p[5];
  TH1F* subleadEtaHist_p[5];
  TH1F* thirdEtaHist_p[5];

  for(Int_t iter = 0; iter < centMax; iter++){
    centString[iter] = getCentString(sType, centLow[iter], centHi[iter]);

    ajHist_p[iter] = new TH1F(Form("%sAJ_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), Form("%sAJ_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), nAJBins, AJLow, AJHi);

    delPhi12Hist_p[iter] = new TH1F(Form("%sDelPhi12_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), Form("%sDelPhi12_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), nDelPhiBins, DelPhiLow, DelPhiHi);
    delPhi13Hist_p[iter] = new TH1F(Form("%sDelPhi13_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), Form("%sDelPhi13_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), nDelPhiBins, DelPhiLow, DelPhiHi);
    delPhi23Hist_p[iter] = new TH1F(Form("%sDelPhi23_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), Form("%sDelPhi23_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), nDelPhiBins, DelPhiLow, DelPhiHi);

    leadPtHist_p[iter] = new TH1F(Form("%sLeadPt_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), Form("%sLeadPt_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), nLeadPtBins, LeadPtLow, LeadPtHi);
    subleadPtHist_p[iter] = new TH1F(Form("%sSubleadPt_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), Form("%sSubleadPt_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), nSubleadPtBins, SubleadPtLow, SubleadPtHi);
    thirdPtHist_p[iter] = new TH1F(Form("%sThirdPt_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), Form("%sThirdPt_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), nThirdPtBins, ThirdPtLow, ThirdPtHi);

    leadEtaHist_p[iter] = new TH1F(Form("%sLeadEta_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), Form("%sLeadEta_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), nEtaBins, EtaLow, EtaHi);
    subleadEtaHist_p[iter] = new TH1F(Form("%sSubleadEta_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), Form("%sSubleadEta_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), nEtaBins, EtaLow, EtaHi);
    thirdEtaHist_p[iter] = new TH1F(Form("%sThirdEta_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), Form("%sThirdEta_%s_%s_h", algType[setNum].c_str(), centString[iter].c_str(), fileTag.c_str()), nEtaBins, EtaLow, EtaHi);
  }

  for(Int_t jEntry = 0; jEntry < (Int_t)anaTree_p->GetEntries(); jEntry++){
    anaTree_p->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(isEventCut(setNum, sType, isHighPtTrk)) continue;
   
    Float_t weight = 1.0;
    if(montecarlo){
      weight = pthatWeight_;
      if(hi) weight *= centWeight_[setNum];
    }

    for(Int_t centIter = 0; centIter < centMax; centIter++){
      if(hiBin_ >= centLow[centIter] && hiBin_ <= centHi[centIter]){
	leadPtHist_p[centIter]->Fill(AlgJtPt_[setNum][0], weight);
	subleadPtHist_p[centIter]->Fill(AlgJtPt_[setNum][1], weight);
	if(AlgJtPt_[setNum][2] > 30 && TMath::Abs(AlgJtEta_[setNum][2]) < 2.0) thirdPtHist_p[centIter]->Fill(AlgJtPt_[setNum][2], weight);


	leadEtaHist_p[centIter]->Fill(AlgJtEta_[setNum][0], weight);
	subleadEtaHist_p[centIter]->Fill(AlgJtEta_[setNum][1], weight);
	if(AlgJtPt_[setNum][2] > 30 && TMath::Abs(AlgJtEta_[setNum][2]) < 2.0) thirdEtaHist_p[centIter]->Fill(AlgJtEta_[setNum][2], weight);

	delPhi12Hist_p[centIter]->Fill(AlgJtDelPhi12_[setNum], weight);

	if(TMath::Abs(AlgJtEta_[setNum][0]) > 1.6 || TMath::Abs(AlgJtEta_[setNum][1]) > 1.6) continue;
	if(TMath::Abs(AlgJtDelPhi12_[setNum]) < 5.0*TMath::Pi()/6) continue;

	ajHist_p[centIter]->Fill(AlgJtAsymm_[setNum], weight);
	break;
      }
    }
  }

  outFile_p = new TFile(outName.c_str(), "UPDATE");
  std::cout << outName << std::endl;

  for(Int_t iter = 0; iter < centMax; iter++){
    leadPtHist_p[iter]->Scale(1./leadPtHist_p[iter]->Integral());
    leadPtHist_p[iter]->Write("", TObject::kOverwrite);
    subleadPtHist_p[iter]->Scale(1./subleadPtHist_p[iter]->Integral());
    subleadPtHist_p[iter]->Write("", TObject::kOverwrite);
    thirdPtHist_p[iter]->Scale(1./thirdPtHist_p[iter]->Integral());
    thirdPtHist_p[iter]->Write("", TObject::kOverwrite);

    leadEtaHist_p[iter]->Scale(1./leadEtaHist_p[iter]->Integral());
    leadEtaHist_p[iter]->Write("", TObject::kOverwrite);
    subleadEtaHist_p[iter]->Scale(1./subleadEtaHist_p[iter]->Integral());
    subleadEtaHist_p[iter]->Write("", TObject::kOverwrite);
    thirdEtaHist_p[iter]->Scale(1./thirdEtaHist_p[iter]->Integral());
    thirdEtaHist_p[iter]->Write("", TObject::kOverwrite);

    delPhi12Hist_p[iter]->Scale(1./delPhi12Hist_p[iter]->Integral());
    delPhi12Hist_p[iter]->Write("", TObject::kOverwrite);

    ajHist_p[iter]->Scale(1./ajHist_p[iter]->Integral());
    ajHist_p[iter]->Write("", TObject::kOverwrite);
  }

  outFile_p->Close();
  delete outFile_p;
  outFile_p = 0;

  for(Int_t iter = 0; iter < centMax; iter++){
    delete leadPtHist_p[iter];
    delete subleadPtHist_p[iter];
    delete thirdPtHist_p[iter];

    delete leadEtaHist_p[iter];
    delete subleadEtaHist_p[iter];
    delete thirdEtaHist_p[iter];

    delete delPhi12Hist_p[iter];
    delete delPhi13Hist_p[iter];
    delete delPhi23Hist_p[iter];
    delete ajHist_p[iter];
  }

  return;
}


void runMakeDiJetHists_Jet(const std::string inName, const std::string outName, sampleType sType = kHIDATA, Bool_t isHighPtTrk = false)
{
  TH1::SetDefaultSumw2();

  outFile_p = new TFile(outName.c_str(), "RECREATE");
  outFile_p->Close();
  delete outFile_p;

  setFileTag(inName);

  inFile_p = new TFile(inName.c_str(), "READ");
  GetDiJetAnaSkim(inFile_p, sType);

  std::cout << "AnaSkim Loaded" << std::endl;

  if(isMonteCarlo(sType))
    trackTreeAna_p->AddFriend(genTreeAna_p);

  jetTreeAna_p->AddFriend(trackTreeAna_p);

  for(Int_t setIter = 0; setIter < 7; setIter++){
    makeDijetHists_Jet(jetTreeAna_p, outName, setIter, sType, isHighPtTrk);
  }

  inFile_p->Close();
  delete inFile_p;

  return;
}
