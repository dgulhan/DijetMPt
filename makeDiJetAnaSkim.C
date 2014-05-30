//=============================================
// Author: Chris McGinn
// 
// DiJet Analysis Skim Class (MC)
//
//=============================================

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "commonUtility.h"
#include "cfmDiJetAnaSkim.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include "factorizedPtCorr.h"
#include "/net/hisrv0001/home/cfmcginn/emDiJet/CMSSW_5_3_12_patch3/tempHIFA/HiForestAnalysis/hiForest.h"
#include "TComplex.h"

const Float_t leadJtPtCut = 120.;
const Float_t subLeadJtPtCut = 50.;
const Float_t jtDelPhiCut = 0;
const Float_t jtEtaCut = 2.0; // Default Max at 2.4 to avoid transition junk, otherwise vary as needed

collisionType getCType(sampleType sType);

const char* algType[3] = {"PuCalo", "VsCalo", "T"};


Float_t getFlippedPhi(Float_t inPhi)
{
  Float_t outPhi;

  if(TMath::Abs(inPhi) > TMath::Pi()){
    std::cout << "getFlippedPhi: inPhi is outside accepted range, return -10" << std::endl;
    return -10;
  }
  else if(inPhi > 0)
    outPhi = inPhi - TMath::Pi();
  else
    outPhi = inPhi + TMath::Pi();

  return outPhi;
}


Bool_t sameSign(Float_t num1, Float_t num2){
  if((num1 > 0 && num2 > 0) || (num1 < 0 && num2 < 0)) return true;

  return false;
}


Float_t getAvePhi(Float_t inLeadPhi, Float_t inSubLeadPhi)
{
  Float_t flipPhi = getFlippedPhi(inSubLeadPhi);
  Float_t avePhi;

  if(sameSign(inLeadPhi, flipPhi) || (TMath::Abs(inLeadPhi) < TMath::Pi()/2 && TMath::Abs(flipPhi) < TMath::Pi()/2))
    avePhi = (flipPhi + inLeadPhi)/2;
  else if(TMath::Abs(inLeadPhi) > TMath::Pi()/2 && TMath::Abs(flipPhi) > TMath::Pi()/2){
    avePhi = (flipPhi + inLeadPhi)/2;
    if(avePhi > 0)
      avePhi = TMath::Pi() - avePhi;
    else
      avePhi = -TMath::Pi() - avePhi;
  }
  else{
    avePhi = 0.;
  }

  return avePhi;
}


void getJtVar(Int_t nJt, Float_t jtPt[], Float_t jtPhi[], Float_t jtEta[], Float_t refPt[], Float_t refPhi[], Float_t refEta[], Int_t algNum, Bool_t montecarlo = false)
{
  if(jtPt[0] < leadJtPtCut)
    return;

  eventSet_[algNum] = true;

  AlgLeadJtPt_[algNum] = jtPt[0];
  AlgLeadJtPhi_[algNum] = jtPhi[0];
  AlgLeadJtEta_[algNum] = jtEta[0];
  AlgSubLeadJtPt_[algNum] = jtPt[1];
  AlgSubLeadJtPhi_[algNum] = jtPhi[1];
  AlgSubLeadJtEta_[algNum] = jtEta[1];

  if(montecarlo && algNum != 2){
    AlgLeadRefPt_[algNum] = refPt[0];
    AlgLeadRefPhi_[algNum] = refPhi[0];
    AlgLeadRefEta_[algNum] = refEta[0];
    AlgSubLeadRefPt_[algNum] = refPt[1];
    AlgSubLeadRefPhi_[algNum] = refPhi[1];
    AlgSubLeadRefEta_[algNum] = refEta[1];
  }

  if(nJt > 2){
    AlgThirdJtPt_[algNum] = jtPt[2];
    AlgThirdJtPhi_[algNum] = jtPhi[2];
    AlgThirdJtEta_[algNum] = jtEta[2];

    if(montecarlo && algNum != 2){
      AlgThirdRefPt_[algNum] = refPt[2];
      AlgThirdRefPhi_[algNum] = refPhi[2];
      AlgThirdRefEta_[algNum] = refEta[2];
    }

    if(nJt > 3){
      AlgFourthJtPt_[algNum] = jtPt[3];
      AlgFourthJtPhi_[algNum] = jtPhi[3];
      AlgFourthJtEta_[algNum] = jtEta[3];

      if(montecarlo && algNum != 2){
        AlgFourthRefPt_[algNum] = refPt[3];
        AlgFourthRefPhi_[algNum] = refPhi[3];
        AlgFourthRefEta_[algNum] = refEta[3];
      }

    }
  }

  AlgJtAvePhi_[algNum] = getAvePhi(AlgLeadJtPhi_[algNum], AlgSubLeadJtPhi_[algNum]);
  AlgJtDelPhi_[algNum] = getAbsDphi(AlgLeadJtPhi_[algNum], AlgSubLeadJtPhi_[algNum]);
  AlgJtAsymm_[algNum] = (AlgLeadJtPt_[algNum] - AlgSubLeadJtPt_[algNum])/(AlgLeadJtPt_[algNum] + AlgSubLeadJtPt_[algNum]);

  if(montecarlo && algNum != 2){
    AlgRefAvePhi_[algNum] = getAvePhi(AlgLeadRefPhi_[algNum], AlgSubLeadRefPhi_[algNum]);
    AlgRefDelPhi_[algNum] = getAbsDphi(AlgLeadRefPhi_[algNum], AlgSubLeadRefPhi_[algNum]);
    AlgRefAsymm_[algNum] = (AlgLeadRefPt_[algNum] - AlgSubLeadRefPt_[algNum])/(AlgLeadRefPt_[algNum] + AlgSubLeadRefPt_[algNum]);
  }

  return;
}



void getPtProj(Float_t cutPt, Float_t inPt, Float_t phi, Float_t jtPhi, Float_t& ProjF, Float_t& Proj0_1, Float_t& Proj1_2, Float_t& Proj2_4, Float_t& Proj4_8, Float_t& Proj8_100)
{
  ProjF += -inPt*cos(getDPHI(phi, jtPhi));

  if(cutPt < 1)
    Proj0_1 += -inPt*cos(getDPHI(phi, jtPhi));
  else if(cutPt < 2)
    Proj1_2 += -inPt*cos(getDPHI(phi, jtPhi));
  else if(cutPt < 4)
    Proj2_4 += -inPt*cos(getDPHI(phi, jtPhi));
  else if(cutPt < 8)
    Proj4_8 += -inPt*cos(getDPHI(phi, jtPhi));
  else
    Proj8_100 += -inPt*cos(getDPHI(phi, jtPhi));

  return;
}


Float_t getTrkRMin(Float_t trkPhi, Float_t trkEta, Int_t nJt, Float_t jtPhi[], Float_t jtEta[])
{
  Float_t trkRMin = 199;

  if(nJt != 0){
    for(Int_t jtEntry = 0; jtEntry < nJt; jtEntry++){
      if(trkRMin > getDR(trkEta, trkPhi, jtEta[jtEntry], jtPhi[jtEntry]))
	trkRMin = getDR(trkEta, trkPhi, jtEta[jtEntry], jtPhi[jtEntry]);
    }
  }

  return trkRMin;
}


int makeDiJetAnaSkim(string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_CFMSKIM.root", Int_t num = 0)
{
  //Define MC or Data
  Bool_t montecarlo = false;
  if(sType == kPPMC || sType == kPAMC || sType == kHIMC)
    montecarlo = true;

  std::cout << sType << std::endl;
  std::cout << montecarlo << std::endl;

  string buffer;
  std::vector<string> listOfFiles;
  int nLines = 0;
  ifstream inFile(fList.data());

  std::cout << fList << std::endl;
  std::cout << inFile.is_open() << std::endl;

  if(!inFile.is_open()){
    std::cout << "Error opening file. Exiting." <<std::endl;
    return 1;
  }
  else{
    while(!inFile.eof()){
      inFile >> buffer;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }

  std::cout << "FileList Loaded" << std::endl;

  TFile* iniSkim_p = new TFile(listOfFiles[0].data(), "READ");

  GetDiJetIniSkim(iniSkim_p, montecarlo);

  std::cout << "IniSkim Loaded" << std::endl;

  //Setup correction tables

  InitCorrFiles();
  InitCorrHists();

  TFile *centHistFile_p = new TFile("centHist_eventSelect.root", "READ");
  TH1F *hist_DataOverMC_p[2];

  TFile *centHistFile_2pi3_p = new TFile("centHist_eventSelect_2pi3.root", "READ");
  TH1F *hist_DataOverMC_2pi3_p[2];

  TFile *centHistFile_120_5pi6_p = new TFile("centHist_eventSelect_120_5pi6.root", "READ");
  TH1F *hist_DataOverMC_120_5pi6_p[2];

  TFile *centHistFile_hatAll_0_p = new TFile("centHist_eventSelect_hatAll_0.root", "READ");
  TH1F *hist_DataOverMC_hatAll_0_p[2];

  TFile *centHistFile_hatAll_1_p = new TFile("centHist_eventSelect_hatAll_1.root", "READ");
  TH1F *hist_DataOverMC_hatAll_1_p[2];

  if(montecarlo){
    for(Int_t algIter = 0; algIter < 2; algIter++){
      hist_DataOverMC_p[algIter] = (TH1F*)centHistFile_p->Get(Form("hiBin_%s_DataOverMC_h", algType[algIter]));
      hist_DataOverMC_2pi3_p[algIter] = (TH1F*)centHistFile_2pi3_p->Get(Form("hiBin_%s_DataOverMC_h", algType[algIter]));
      hist_DataOverMC_120_5pi6_p[algIter] = (TH1F*)centHistFile_120_5pi6_p->Get(Form("hiBin_%s_DataOverMC_h", algType[algIter]));
      hist_DataOverMC_hatAll_0_p[algIter] = (TH1F*)centHistFile_hatAll_0_p->Get(Form("hiBin_%s_DataOverMC_h", algType[algIter]));
      hist_DataOverMC_hatAll_1_p[algIter] = (TH1F*)centHistFile_hatAll_1_p->Get(Form("hiBin_%s_DataOverMC_h", algType[algIter]));
    }
  }

  TFile *outFile = new TFile(Form("%s_%d.root", outName, num), "RECREATE");

  InitDiJetAnaSkim(montecarlo);

  Long64_t nentries = trackTreeIni_p->GetEntries();

  std::cout << nentries << std::endl;

  std::cout << "Cuts, Lead/Sublead Pt, delphi, eta: " << leadJtPtCut << ", " << subLeadJtPtCut << ", " << jtDelPhiCut << ", " << jtEtaCut << std::endl; 

  for(Long64_t jentry = 0; jentry < nentries; jentry++){
    trackTreeIni_p->GetEntry(jentry);
    jetTreeIni_p->GetEntry(jentry);

    if(montecarlo)
      genTreeIni_p->GetEntry(jentry);

    if(jentry%1000 == 0)
      std::cout << jentry << std::endl;

    InitJetVar(montecarlo);

    //Jet Edits here

    if(nPu3Calo_ > 1){
      getJtVar(nPu3Calo_, Pu3CaloPt_, Pu3CaloPhi_, Pu3CaloEta_, Pu3CaloRefPt_, Pu3CaloRefPhi_, Pu3CaloRefEta_, 0, montecarlo);
    }

    if(nVs3Calo_ > 1){
      getJtVar(nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_, Vs3CaloRefPt_, Vs3CaloRefPhi_, Vs3CaloRefEta_, 1, montecarlo);
    }

    Float_t dummyArray[nT3_];

    if(nT3_ > 1){
      getJtVar(nT3_, T3Pt_, T3Phi_, T3Eta_, dummyArray, dummyArray, dummyArray, 2, montecarlo);
    }
    
    if(eventSet_[PuCalo] == false && eventSet_[VsCalo] == false && eventSet_[T] == false){
      std::cout << "No event pass after IniSkim; Potential bug" << std::endl;
      continue;
    }

    if(montecarlo){
      Float_t pthatWeights[5] = {4.29284e-01, 2.99974e-02, 3.38946e-04, 1.06172e-04, 2.79631e-05};
      Float_t pthatCuts[6] = {30, 50, 80, 100, 120, 100000};

      for(Int_t hatIter = 0; hatIter < 5; hatIter++){
	if(pthat_ > pthatCuts[hatIter] && pthat_ < pthatCuts[hatIter + 1]){
	  pthatWeight_ = pthatWeights[hatIter];
	  break;
	}
      }
    }

    if(montecarlo){
      for(Int_t algIter = 0; algIter < 2; algIter++){
	centWeight_[algIter] = hist_DataOverMC_p[algIter]->GetBinContent(hist_DataOverMC_p[algIter]->FindBin(hiBin_));
	centWeight_2pi3_[algIter] = hist_DataOverMC_2pi3_p[algIter]->GetBinContent(hist_DataOverMC_2pi3_p[algIter]->FindBin(hiBin_));
	centWeight_120_5pi6_[algIter] = hist_DataOverMC_120_5pi6_p[algIter]->GetBinContent(hist_DataOverMC_120_5pi6_p[algIter]->FindBin(hiBin_));
	centWeight_hatAll_0_[algIter] = hist_DataOverMC_hatAll_0_p[algIter]->GetBinContent(hist_DataOverMC_hatAll_0_p[algIter]->FindBin(hiBin_));
	centWeight_hatAll_1_[algIter] = hist_DataOverMC_hatAll_1_p[algIter]->GetBinContent(hist_DataOverMC_hatAll_1_p[algIter]->FindBin(hiBin_));

	fullWeight_[algIter] = centWeight_hatAll_1_[algIter]*pthatWeight_;

	if(algIter == 3)
	  std::cout << "Hibin, pos, centWeight: " << hiBin_ << ", " << hist_DataOverMC_hatAll_1_p[algIter]->FindBin(hiBin_) << ", " << centWeight_hatAll_1_[algIter] << std::endl;  

      }
    }

    //Iterate over tracks

    InitProjPerp(montecarlo);

    for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
          
      //Grab proj. Pt Spectra For Tracks in each Event Subset
    
      for(Int_t jtIter = 0; jtIter < 3; jtIter++){
	if(eventSet_[jtIter]){
	  getPtProj(trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjAF_[jtIter], rAlgImbProjA0_1_[jtIter], rAlgImbProjA1_2_[jtIter], rAlgImbProjA2_4_[jtIter], rAlgImbProjA4_8_[jtIter], rAlgImbProjA8_100_[jtIter]);

	  Float_t tempLeadDelR = getDR(trkEta_[trkEntry], trkPhi_[trkEntry], AlgLeadJtEta_[jtIter], AlgLeadJtPhi_[jtIter]);
	  Float_t tempSubLeadDelR = getDR(trkEta_[trkEntry], trkPhi_[trkEntry], AlgSubLeadJtEta_[jtIter], AlgSubLeadJtPhi_[jtIter]);

	  if(tempLeadDelR > 0 && tempSubLeadDelR > 0){
	    if(tempLeadDelR < .8 || tempSubLeadDelR < .8)
	      getPtProj(trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjACF_[jtIter], rAlgImbProjAC0_1_[jtIter], rAlgImbProjAC1_2_[jtIter], rAlgImbProjAC2_4_[jtIter], rAlgImbProjAC4_8_[jtIter], rAlgImbProjAC8_100_[jtIter]); 
     	    else{
	      getPtProj(trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjANCF_[jtIter], rAlgImbProjANC0_1_[jtIter], rAlgImbProjANC1_2_[jtIter], rAlgImbProjANC2_4_[jtIter], rAlgImbProjANC4_8_[jtIter], rAlgImbProjANC8_100_[jtIter]);
	    }
	    
	    if(tempLeadDelR < .20 || tempSubLeadDelR < .20)
	      getPtProj(trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA1CF_[jtIter], rAlgImbProjA1C0_1_[jtIter], rAlgImbProjA1C1_2_[jtIter], rAlgImbProjA1C2_4_[jtIter], rAlgImbProjA1C4_8_[jtIter], rAlgImbProjA1C8_100_[jtIter]);
	    else if(tempLeadDelR < .40 || tempSubLeadDelR < .40)
	      getPtProj(trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA2CF_[jtIter], rAlgImbProjA2C0_1_[jtIter], rAlgImbProjA2C1_2_[jtIter], rAlgImbProjA2C2_4_[jtIter], rAlgImbProjA2C4_8_[jtIter], rAlgImbProjA2C8_100_[jtIter]);
	    else if(tempLeadDelR < .60 || tempSubLeadDelR < .60)
	      getPtProj(trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA3CF_[jtIter], rAlgImbProjA3C0_1_[jtIter], rAlgImbProjA3C1_2_[jtIter], rAlgImbProjA3C2_4_[jtIter], rAlgImbProjA3C4_8_[jtIter], rAlgImbProjA3C8_100_[jtIter]);
	    else if(tempLeadDelR < .80 || tempSubLeadDelR < .80)
	      getPtProj(trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA4CF_[jtIter], rAlgImbProjA4C0_1_[jtIter], rAlgImbProjA4C1_2_[jtIter], rAlgImbProjA4C2_4_[jtIter], rAlgImbProjA4C4_8_[jtIter], rAlgImbProjA4C8_100_[jtIter]);
	    else if(tempLeadDelR < 1.0 || tempSubLeadDelR < 1.0)
	      getPtProj(trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA5CF_[jtIter], rAlgImbProjA5C0_1_[jtIter], rAlgImbProjA5C1_2_[jtIter], rAlgImbProjA5C2_4_[jtIter], rAlgImbProjA5C4_8_[jtIter], rAlgImbProjA5C8_100_[jtIter]);
	    else if(tempLeadDelR < 1.2 || tempSubLeadDelR < 1.2)
	      getPtProj(trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA6CF_[jtIter], rAlgImbProjA6C0_1_[jtIter], rAlgImbProjA6C1_2_[jtIter], rAlgImbProjA6C2_4_[jtIter], rAlgImbProjA6C4_8_[jtIter], rAlgImbProjA6C8_100_[jtIter]);
	    else if(tempLeadDelR < 1.4 || tempSubLeadDelR < 1.4)
	      getPtProj(trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA7CF_[jtIter], rAlgImbProjA7C0_1_[jtIter], rAlgImbProjA7C1_2_[jtIter], rAlgImbProjA7C2_4_[jtIter], rAlgImbProjA7C4_8_[jtIter], rAlgImbProjA7C8_100_[jtIter]);
	    else if(tempLeadDelR < 1.6 || tempSubLeadDelR < 1.6)
	      getPtProj(trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA8CF_[jtIter], rAlgImbProjA8C0_1_[jtIter], rAlgImbProjA8C1_2_[jtIter], rAlgImbProjA8C2_4_[jtIter], rAlgImbProjA8C4_8_[jtIter], rAlgImbProjA8C8_100_[jtIter]);
	    else if(tempLeadDelR < 1.8 || tempSubLeadDelR < 1.8)
	      getPtProj(trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA9CF_[jtIter], rAlgImbProjA9C0_1_[jtIter], rAlgImbProjA9C1_2_[jtIter], rAlgImbProjA9C2_4_[jtIter], rAlgImbProjA9C4_8_[jtIter], rAlgImbProjA9C8_100_[jtIter]);
	    else 
	      getPtProj(trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], AlgJtAvePhi_[jtIter], rAlgImbProjA10CF_[jtIter], rAlgImbProjA10C0_1_[jtIter], rAlgImbProjA10C1_2_[jtIter], rAlgImbProjA10C2_4_[jtIter], rAlgImbProjA10C4_8_[jtIter], rAlgImbProjA10C8_100_[jtIter]); 
	    
	  }
	}
      }
    }
    
    Int_t hiBinDiv[5] = {20, 40, 60, 100, 200};
    Int_t hiSetEff[30] = {0, 5, 10, 15, 20, 25, 1, 6, 11, 16, 21, 26, 2, 7, 12, 17, 22, 27, 3, 8, 13, 18, 23, 27, 4, 9, 14, 19, 24, 27};
  
    for(Int_t hiBinIter = 0; hiBinIter < 5; hiBinIter++){
      if(hiBin_ < hiBinDiv[hiBinIter]){
	for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	  Int_t ptPos = getPtBin(trkPt_[trkEntry], hiSetEff[hiBinIter*6], hiSetEff[hiBinIter*6 + 1], hiSetEff[hiBinIter*6 + 2], hiSetEff[hiBinIter*6 + 3], hiSetEff[hiBinIter*6 + 4], hiSetEff[hiBinIter*6 + 5], 28);
	
	  //do temp rmin here

	  Float_t tempRMin[3];

	  tempRMin[PuCalo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nPu3Calo_, Pu3CaloPhi_, Pu3CaloEta_);
	  tempRMin[VsCalo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3Calo_, Vs3CaloPhi_, Vs3CaloEta_);
	  tempRMin[T] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nT3_, T3Phi_, T3Eta_);

	  Float_t tempFact[3];
	  Float_t tempCorr[3];


	  tempFact[PuCalo] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[PuCalo], VsCalocent_p[ptPos], VsCalophiEta_p[ptPos], VsCalopt_p[ptPos], VsCalodelR_p[ptPos]);
	  tempFact[PuCalo] = (1 - factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[PuCalo], FakeVsCalocent_p[ptPos], FakeVsCalophiEta_p[ptPos], FakeVsCalopt_p[ptPos], FakeVsCalodelR_p[ptPos], false))/tempFact[PuCalo];
	  tempCorr[PuCalo] = trkPt_[trkEntry]*tempFact[PuCalo];

	  tempFact[VsCalo] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[VsCalo], VsCalocent_p[ptPos], VsCalophiEta_p[ptPos], VsCalopt_p[ptPos], VsCalodelR_p[ptPos]);
	  tempFact[VsCalo] = (1 - factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[VsCalo], FakeVsCalocent_p[ptPos], FakeVsCalophiEta_p[ptPos], FakeVsCalopt_p[ptPos], FakeVsCalodelR_p[ptPos], false))/tempFact[VsCalo];
	  tempCorr[VsCalo] = trkPt_[trkEntry]*tempFact[VsCalo];


	  if(montecarlo){
	    tempFact[T] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[T], VsCalocent_p[ptPos], VsCalophiEta_p[ptPos], VsCalopt_p[ptPos], VsCalodelR_p[ptPos]);
	    tempCorr[T] = trkPt_[trkEntry]*(1 - factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[T], FakeVsCalocent_p[ptPos], FakeVsCalophiEta_p[ptPos], FakeVsCalopt_p[ptPos], FakeVsCalodelR_p[ptPos], false))/tempFact[T];
	  }

	  for(Int_t setIter = 0; setIter < 3; setIter++){
	    if(eventSet_[setIter]){

	      Float_t tempLeadR = getDR(trkEta_[trkEntry], trkPhi_[trkEntry], AlgLeadJtEta_[setIter], AlgLeadJtPhi_[setIter]);
	      Float_t tempSubLeadR = getDR(trkEta_[trkEntry], trkPhi_[trkEntry], AlgSubLeadJtEta_[setIter], AlgSubLeadJtPhi_[setIter]);


	      getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjAF_[setIter + 3], rAlgImbProjA0_1_[setIter + 3], rAlgImbProjA1_2_[setIter + 3], rAlgImbProjA2_4_[setIter + 3], rAlgImbProjA4_8_[setIter + 3], rAlgImbProjA8_100_[setIter + 3]);

	      if(tempLeadR > 0 && tempSubLeadR > 0){
		if(tempLeadR < .8 || tempSubLeadR < .8)
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjACF_[setIter + 3], rAlgImbProjAC0_1_[setIter + 3], rAlgImbProjAC1_2_[setIter + 3], rAlgImbProjAC2_4_[setIter + 3], rAlgImbProjAC4_8_[setIter + 3], rAlgImbProjAC8_100_[setIter + 3]);
		else
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjANCF_[setIter + 3], rAlgImbProjANC0_1_[setIter + 3], rAlgImbProjANC1_2_[setIter + 3], rAlgImbProjANC2_4_[setIter + 3], rAlgImbProjANC4_8_[setIter + 3], rAlgImbProjANC8_100_[setIter + 3]);		

		if(tempLeadR < .20 || tempSubLeadR < .20)
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA1CF_[setIter + 3], rAlgImbProjA1C0_1_[setIter + 3], rAlgImbProjA1C1_2_[setIter + 3], rAlgImbProjA1C2_4_[setIter + 3], rAlgImbProjA1C4_8_[setIter + 3], rAlgImbProjA1C8_100_[setIter + 3]);
		else if(tempLeadR < .40 || tempSubLeadR < .40)
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA2CF_[setIter + 3], rAlgImbProjA2C0_1_[setIter + 3], rAlgImbProjA2C1_2_[setIter + 3], rAlgImbProjA2C2_4_[setIter + 3], rAlgImbProjA2C4_8_[setIter + 3], rAlgImbProjA2C8_100_[setIter + 3]);
		else if(tempLeadR < .60 || tempSubLeadR < .60)
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA3CF_[setIter + 3], rAlgImbProjA3C0_1_[setIter + 3], rAlgImbProjA3C1_2_[setIter + 3], rAlgImbProjA3C2_4_[setIter + 3], rAlgImbProjA3C4_8_[setIter + 3], rAlgImbProjA3C8_100_[setIter + 3]);
		else if(tempLeadR < .80 || tempSubLeadR < .80)
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA4CF_[setIter + 3], rAlgImbProjA4C0_1_[setIter + 3], rAlgImbProjA4C1_2_[setIter + 3], rAlgImbProjA4C2_4_[setIter + 3], rAlgImbProjA4C4_8_[setIter + 3], rAlgImbProjA4C8_100_[setIter + 3]);
		else if(tempLeadR < 1.0 || tempSubLeadR < 1.0)
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA5CF_[setIter + 3], rAlgImbProjA5C0_1_[setIter + 3], rAlgImbProjA5C1_2_[setIter + 3], rAlgImbProjA5C2_4_[setIter + 3], rAlgImbProjA5C4_8_[setIter + 3], rAlgImbProjA5C8_100_[setIter + 3]);
		else if(tempLeadR < 1.2 || tempSubLeadR < 1.2)
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA6CF_[setIter + 3], rAlgImbProjA6C0_1_[setIter + 3], rAlgImbProjA6C1_2_[setIter + 3], rAlgImbProjA6C2_4_[setIter + 3], rAlgImbProjA6C4_8_[setIter + 3], rAlgImbProjA6C8_100_[setIter + 3]);
		else if(tempLeadR < 1.4 || tempSubLeadR < 1.4)
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA7CF_[setIter + 3], rAlgImbProjA7C0_1_[setIter + 3], rAlgImbProjA7C1_2_[setIter + 3], rAlgImbProjA7C2_4_[setIter + 3], rAlgImbProjA7C4_8_[setIter + 3], rAlgImbProjA7C8_100_[setIter + 3]);
		else if(tempLeadR < 1.6 || tempSubLeadR < 1.6)
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA8CF_[setIter + 3], rAlgImbProjA8C0_1_[setIter + 3], rAlgImbProjA8C1_2_[setIter + 3], rAlgImbProjA8C2_4_[setIter + 3], rAlgImbProjA8C4_8_[setIter + 3], rAlgImbProjA8C8_100_[setIter + 3]);
		else if(tempLeadR < 1.8 || tempSubLeadR < 1.8)
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA9CF_[setIter + 3], rAlgImbProjA9C0_1_[setIter + 3], rAlgImbProjA9C1_2_[setIter + 3], rAlgImbProjA9C2_4_[setIter + 3], rAlgImbProjA9C4_8_[setIter + 3], rAlgImbProjA9C8_100_[setIter + 3]);
		else
		  getPtProj(trkPt_[trkEntry], tempCorr[setIter], trkPhi_[trkEntry], AlgJtAvePhi_[setIter], rAlgImbProjA10CF_[setIter + 3], rAlgImbProjA10C0_1_[setIter + 3], rAlgImbProjA10C1_2_[setIter + 3], rAlgImbProjA10C2_4_[setIter + 3], rAlgImbProjA10C4_8_[setIter + 3], rAlgImbProjA10C8_100_[setIter + 3]);

		    
	      }
	    }
	  }  
	}
	break;
      }
    }
    
    if(montecarlo){
      //Iterate over truth

      for(Int_t genEntry = 0; genEntry < nGen_; genEntry++){            
      
	for(Int_t setIter = 0; setIter < 3; setIter++){
	  if(eventSet_[setIter]){
	    getPtProj(genPt_[genEntry], genPt_[genEntry], genPhi_[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjAF_[setIter], gAlgImbProjA0_1_[setIter], gAlgImbProjA1_2_[setIter], gAlgImbProjA2_4_[setIter], gAlgImbProjA4_8_[setIter], gAlgImbProjA8_100_[setIter]);

	  
	    Float_t tempLeadDelR = getDR(genEta_[genEntry], genPhi_[genEntry], AlgLeadJtEta_[setIter], AlgLeadJtPhi_[setIter]);
	    Float_t tempSubLeadDelR = getDR(genEta_[genEntry], genPhi_[genEntry], AlgSubLeadJtEta_[setIter], AlgSubLeadJtPhi_[setIter]);

	    if(tempLeadDelR > 0 && tempSubLeadDelR > 0){

	      if(tempLeadDelR < .8 || tempSubLeadDelR < .8)
		getPtProj(genPt_[genEntry], genPt_[genEntry], genPhi_[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjACF_[setIter], gAlgImbProjAC0_1_[setIter], gAlgImbProjAC1_2_[setIter], gAlgImbProjAC2_4_[setIter], gAlgImbProjAC4_8_[setIter], gAlgImbProjAC8_100_[setIter]);
	      else
		getPtProj(genPt_[genEntry], genPt_[genEntry], genPhi_[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjANCF_[setIter], gAlgImbProjANC0_1_[setIter], gAlgImbProjANC1_2_[setIter], gAlgImbProjANC2_4_[setIter], gAlgImbProjANC4_8_[setIter], gAlgImbProjANC8_100_[setIter]);		

	      
	      if(tempLeadDelR < .20 || tempSubLeadDelR < .20)
		getPtProj(genPt_[genEntry], genPt_[genEntry], genPhi_[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA1CF_[setIter], gAlgImbProjA1C0_1_[setIter], gAlgImbProjA1C1_2_[setIter], gAlgImbProjA1C2_4_[setIter], gAlgImbProjA1C4_8_[setIter], gAlgImbProjA1C8_100_[setIter]);
	      else if(tempLeadDelR < .40 || tempSubLeadDelR < .40)
		getPtProj(genPt_[genEntry], genPt_[genEntry], genPhi_[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA2CF_[setIter], gAlgImbProjA2C0_1_[setIter], gAlgImbProjA2C1_2_[setIter], gAlgImbProjA2C2_4_[setIter], gAlgImbProjA2C4_8_[setIter], gAlgImbProjA2C8_100_[setIter]);
	      else if(tempLeadDelR < .60 || tempSubLeadDelR < .60)
		getPtProj(genPt_[genEntry], genPt_[genEntry], genPhi_[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA3CF_[setIter], gAlgImbProjA3C0_1_[setIter], gAlgImbProjA3C1_2_[setIter], gAlgImbProjA3C2_4_[setIter], gAlgImbProjA3C4_8_[setIter], gAlgImbProjA3C8_100_[setIter]);
	      else if(tempLeadDelR < .80 || tempSubLeadDelR < .80)
		getPtProj(genPt_[genEntry], genPt_[genEntry], genPhi_[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA4CF_[setIter], gAlgImbProjA4C0_1_[setIter], gAlgImbProjA4C1_2_[setIter], gAlgImbProjA4C2_4_[setIter], gAlgImbProjA4C4_8_[setIter], gAlgImbProjA4C8_100_[setIter]);
	      else if(tempLeadDelR < 1.0 || tempSubLeadDelR < 1.0)
		getPtProj(genPt_[genEntry], genPt_[genEntry], genPhi_[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA5CF_[setIter], gAlgImbProjA5C0_1_[setIter], gAlgImbProjA5C1_2_[setIter], gAlgImbProjA5C2_4_[setIter], gAlgImbProjA5C4_8_[setIter], gAlgImbProjA5C8_100_[setIter]);
	      else if(tempLeadDelR < 1.2 || tempSubLeadDelR < 1.2)
		getPtProj(genPt_[genEntry], genPt_[genEntry], genPhi_[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA6CF_[setIter], gAlgImbProjA6C0_1_[setIter], gAlgImbProjA6C1_2_[setIter], gAlgImbProjA6C2_4_[setIter], gAlgImbProjA6C4_8_[setIter], gAlgImbProjA6C8_100_[setIter]);
	      else if(tempLeadDelR < 1.4 || tempSubLeadDelR < 1.4)
		getPtProj(genPt_[genEntry], genPt_[genEntry], genPhi_[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA7CF_[setIter], gAlgImbProjA7C0_1_[setIter], gAlgImbProjA7C1_2_[setIter], gAlgImbProjA7C2_4_[setIter], gAlgImbProjA7C4_8_[setIter], gAlgImbProjA7C8_100_[setIter]);
	      else if(tempLeadDelR < 1.6 || tempSubLeadDelR < 1.6)
		getPtProj(genPt_[genEntry], genPt_[genEntry], genPhi_[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA8CF_[setIter], gAlgImbProjA8C0_1_[setIter], gAlgImbProjA8C1_2_[setIter], gAlgImbProjA8C2_4_[setIter], gAlgImbProjA8C4_8_[setIter], gAlgImbProjA8C8_100_[setIter]);
	      else if(tempLeadDelR < 1.8 || tempSubLeadDelR < 1.8)
		getPtProj(genPt_[genEntry], genPt_[genEntry], genPhi_[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA9CF_[setIter], gAlgImbProjA9C0_1_[setIter], gAlgImbProjA9C1_2_[setIter], gAlgImbProjA9C2_4_[setIter], gAlgImbProjA9C4_8_[setIter], gAlgImbProjA9C8_100_[setIter]);
	      else
		getPtProj(genPt_[genEntry], genPt_[genEntry], genPhi_[genEntry], AlgJtAvePhi_[setIter], gAlgImbProjA10CF_[setIter], gAlgImbProjA10C0_1_[setIter], gAlgImbProjA10C1_2_[setIter], gAlgImbProjA10C2_4_[setIter], gAlgImbProjA10C4_8_[setIter], gAlgImbProjA10C8_100_[setIter]);


	    }
	  }
	}
      }
    }

    jetTreeAna_p->Fill();
    trackTreeAna_p->Fill();

    if(montecarlo)
      genTreeAna_p->Fill();
  }

  outFile->cd();
  jetTreeAna_p->Write();
  trackTreeAna_p->Write();

  if(montecarlo)
    genTreeAna_p->Write();

  outFile->Close();

  delete outFile;
  delete centHistFile_p;

  printf("Done.\n");
  return(0);
}


collisionType getCType(sampleType sType)
{
  switch (sType)
    {
    case kPPDATA:
    case kPPMC:
      return cPP;
    case kPADATA:
    case kPAMC:
      return cPPb;
    case kHIDATA:
    case kHIMC:
      return cPbPb;
    }
  return cPbPb; //probably a bad guess
}


int main(int argc, char *argv[])
{
  if(argc != 5)
    {
      std::cout << "Usage: sortForest <inputFile> <MCBool> <outputFile> <#>" << std::endl;
      return 1;
    }

  int rStatus = -1;

  rStatus = makeDiJetAnaSkim(argv[1], sampleType(atoi(argv[2])), argv[3], atoi(argv[4]));

  return rStatus;
}
