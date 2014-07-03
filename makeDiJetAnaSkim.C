//=============================================
// Author: Chris McGinn
// 
// DiJet Analysis Skim Class (MC)
//
//=============================================

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "cfmDiJetAnaSkim.h"
#include "stdlib.h"
#include <fstream>
#include "factorizedPtCorr.h"
#include "/net/hisrv0001/home/cfmcginn/emDiJet/CMSSW_5_3_12_patch3/tempHIFA/HiForestAnalysis/hiForest.h"
#include "TComplex.h"
#include "commonUtility.h"

const Float_t leadJtPtCut = 120.;
const Float_t subLeadJtPtCut = 50.;
const Float_t jtDelPhiCut = 0;
const Float_t jtEtaCut = 2.0; // Default Max at 2.4 to avoid transition junk, otherwise vary as needed

collisionType getCType(sampleType sType);

const char* algType[5] = {"PuCalo", "VsCalo", "T", "PuPF", "VsPF"};


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

void getJtVar(Int_t nJt, Float_t jtPt[], Float_t jtPhi[], Float_t jtEta[], Float_t trkMax[], Float_t rawPt[], Float_t refPt[], Float_t refPhi[], Float_t refEta[], Int_t algNum, Bool_t montecarlo = false)
{
  if(nJt == 0 || nJt == 1 || jtPt[0] < leadJtPtCut || jtPt[1] < subLeadJtPtCut)
    return;

  eventSet_[algNum] = true;

  Int_t iterMax = nJt;
  if(nJt > 4) iterMax = 4;

  for(Int_t leadIter = 0; leadIter < iterMax; leadIter++){
    AlgJtPt_[algNum][leadIter] = jtPt[leadIter];
    AlgJtPhi_[algNum][leadIter] = jtPhi[leadIter];
    AlgJtEta_[algNum][leadIter] = jtEta[leadIter];

    if(algNum!=2){
      AlgJtTrkMax_[algNum][leadIter] = trkMax[leadIter];
      AlgJtRawPt_[algNum][leadIter] = rawPt[leadIter];

      if(montecarlo){
	AlgRefPt_[algNum][leadIter] = refPt[leadIter];
	AlgRefPhi_[algNum][leadIter] = refPhi[leadIter];
	AlgRefEta_[algNum][leadIter] = refEta[leadIter];
      }
    }
  }

  AlgJtAvePhi_[algNum] = getAvePhi(AlgJtPhi_[algNum][0], AlgJtPhi_[algNum][1]);
  AlgJtDelPhi_[algNum] = getAbsDphi(AlgJtPhi_[algNum][0], AlgJtPhi_[algNum][1]);
  AlgJtAsymm_[algNum] = (AlgJtPt_[algNum][0] - AlgJtPt_[algNum][1])/(AlgJtPt_[algNum][0] + AlgJtPt_[algNum][1]);

  if(montecarlo && algNum != 2){
    if(AlgRefPhi_[algNum][0] > -10 && AlgRefPhi_[algNum][1] > -10){
      AlgRefAvePhi_[algNum] = getAvePhi(AlgRefPhi_[algNum][0], AlgRefPhi_[algNum][1]);
      AlgRefDelPhi_[algNum] = getAbsDphi(AlgRefPhi_[algNum][0], AlgRefPhi_[algNum][1]);
      AlgRefAsymm_[algNum] = (AlgRefPt_[algNum][0] - AlgRefPt_[algNum][1])/(AlgRefPt_[algNum][0] + AlgRefPt_[algNum][1]);
    }
    else{
      AlgRefAvePhi_[algNum] = -10;
      AlgRefDelPhi_[algNum] = -10;
      AlgRefAsymm_[algNum] = -10;
    }
  }

  return;
}


Int_t getPtRange(Float_t cutPt)
{
  if(cutPt < 1)
    return 0; 
  else if(cutPt < 2)
    return 1;
  else if(cutPt < 4)
    return 2;
  else if(cutPt < 8)
    return 3;
  else
    return 4;
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


int makeDiJetAnaSkim(string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_DIJETANASKIM.root", Int_t num = 0)
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

  GetDiJetIniSkim(iniSkim_p, montecarlo, sType);

  std::cout << "IniSkim Loaded" << std::endl;

  //Setup correction tables

  InitCorrFiles(sType);
  InitCorrHists(sType);

  TFile *centHistFile_80_p = new TFile("centHist_eventSet_80.root", "READ");
  TH1F *hist_DataOverMC_80_p[5];

  TFile *centHistFile_Merge_p = new TFile("centHist_eventSet_Merge.root", "READ");
  TH1F *hist_DataOverMC_Merge_p[5];

  if(sType == kHIMC){
    for(Int_t algIter = 0; algIter < 5; algIter++){
      if(algIter != 2){
	hist_DataOverMC_80_p[algIter] = (TH1F*)centHistFile_80_p->Get(Form("hiBin_%s_DataOverMC_h", algType[algIter]));
	if(algIter < 2)
	  hist_DataOverMC_Merge_p[algIter] = (TH1F*)centHistFile_Merge_p->Get(Form("hiBin_%s_DataOverMC_h", algType[algIter]));
      }
    }
  } 

  TFile *outFile = new TFile(Form("%s_%d.root", outName, num), "RECREATE");

  InitDiJetAnaSkim(sType);

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

    InitJetVar(sType);

    getJtVar(nPu3Calo_, Pu3CaloPt_, Pu3CaloPhi_, Pu3CaloEta_, Pu3CaloTrkMax_, Pu3CaloRawPt_, Pu3CaloRefPt_, Pu3CaloRefPhi_, Pu3CaloRefEta_, 0, montecarlo);
    getJtVar(nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_, Vs3CaloTrkMax_, Vs3CaloRawPt_, Vs3CaloRefPt_, Vs3CaloRefPhi_, Vs3CaloRefEta_, 1, montecarlo);
    Float_t dummyArray[nT3_];
    getJtVar(nT3_, T3Pt_, T3Phi_, T3Eta_, dummyArray, dummyArray, dummyArray, dummyArray, dummyArray, 2, montecarlo);
    getJtVar(nPu3PF_, Pu3PFPt_, Pu3PFPhi_, Pu3PFEta_, Pu3PFTrkMax_, Pu3PFRawPt_, Pu3PFRefPt_, Pu3PFRefPhi_, Pu3PFRefEta_, 3, montecarlo);
    getJtVar(nVs3PF_, Vs3PFPt_, Vs3PFPhi_, Vs3PFEta_, Vs3PFTrkMax_, Vs3PFRawPt_, Vs3PFRefPt_, Vs3PFRefPhi_, Vs3PFRefEta_, 4, montecarlo);
    
    if(eventSet_[PuCalo] == false && eventSet_[VsCalo] == false && eventSet_[T] == false && eventSet_[PuPF] == false && eventSet_[VsPF] == false){
      std::cout << "No event pass after IniSkim; Potential bug" << std::endl;

      std::cout << "PuPF Lead Jt: " << Pu3PFPt_[0] << ", " << Pu3PFPhi_[0] << ", " << Pu3PFEta_[0] << std::endl;
      std::cout << "PuPF subLead Jt: " << Pu3PFPt_[1] << ", " << Pu3PFPhi_[1] << ", " << Pu3PFEta_[1] << std::endl;

      std::cout << "VsPF Lead Jt: " << Vs3PFPt_[0] << ", " << Vs3PFPhi_[0] << ", " << Vs3PFEta_[0] << std::endl;
      std::cout << "VsPF subLead Jt: " << Vs3PFPt_[1] << ", " << Vs3PFPhi_[1] << ", " << Vs3PFEta_[1] << std::endl;

      continue;
    }

    run_ = runIni_;
    evt_ = evtIni_;
    lumi_ = lumiIni_;

    if(montecarlo)
      pthat_ = pthatIni_;
    
    if(sType == kHIDATA || sType == kHIMC){
      hiBin_ = hiBinIni_;
      hiEvtPlane_ = hiEvtPlaneIni_;
      psin_ = psinIni_;
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

  
    if(sType == kHIMC){
      for(Int_t algIter = 0; algIter < 5; algIter++){
	if(algIter != 2){
	  centWeight_80_[algIter] = hist_DataOverMC_80_p[algIter]->GetBinContent(hist_DataOverMC_80_p[algIter]->FindBin(hiBin_));
	  if(algIter < 2)
	    centWeight_Merge_[algIter] = hist_DataOverMC_Merge_p[algIter]->GetBinContent(hist_DataOverMC_Merge_p[algIter]->FindBin(hiBin_));
	}
      }
    }

    //Iterate over tracks

    InitProjPerp(montecarlo);

    //Switch below to iterated OR EDIT HERE

    if(eventSet_[PuCalo] || eventSet_[VsCalo] || eventSet_[T]){
      
      Float_t rBounds[10] = {.20, .40, .60, .80, 1.00, 1.20, 1.40, 1.60, 1.80, 100000};
      
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
        
	//Grab proj. Pt Spectra For Tracks in each Event Subset
	
	Int_t ptIter = getPtRange(trkPt_[trkEntry]);
	
	for(Int_t jtIter = 0; jtIter < 3; jtIter++){
	  if(eventSet_[jtIter]){
	    
	    rAlgImbProjA_[jtIter][5] += - trkPt_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[jtIter]));
	    rAlgImbProjA_[jtIter][ptIter] += - trkPt_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[jtIter]));
	    
	    Float_t tempLeadDelR = getDR(trkEta_[trkEntry], trkPhi_[trkEntry], AlgJtEta_[jtIter][0], AlgJtPhi_[jtIter][0]);
	    Float_t tempSubLeadDelR = getDR(trkEta_[trkEntry], trkPhi_[trkEntry], AlgJtEta_[jtIter][1], AlgJtPhi_[jtIter][1]);
	    
	    if(tempLeadDelR > 0 && tempSubLeadDelR > 0){
	      if(tempLeadDelR < .8 || tempSubLeadDelR < .8){
		rAlgImbProjAC_[jtIter][5] += - trkPt_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[jtIter]));
		rAlgImbProjAC_[jtIter][ptIter] += - trkPt_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[jtIter]));
	      }
	      else{
		rAlgImbProjANC_[jtIter][5] += - trkPt_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[jtIter]));
		rAlgImbProjANC_[jtIter][ptIter] += - trkPt_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[jtIter]));
	      }
	      
	      for(Int_t rIter = 0; rIter < 10; rIter++){
		if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
		  rAlgImbProjAR_[jtIter][5][rIter] += - trkPt_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[jtIter]));
		  rAlgImbProjAR_[jtIter][ptIter][rIter] += - trkPt_[trkEntry]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[jtIter]));
		  break;
		}
	      }
	      
	    }
	  }
	}
      }
      
      if(sType == kHIDATA || sType == kHIMC)
	InitPosArrPbPb(hiBin_);
      
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	Int_t ptPos = getPtBin(trkPt_[trkEntry], sType);
	
	Float_t tempRMin[3];
	
	tempRMin[PuCalo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nPu3Calo_, Pu3CaloPhi_, Pu3CaloEta_);
	tempRMin[VsCalo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3Calo_, Vs3CaloPhi_, Vs3CaloEta_);
	tempRMin[T] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nT3_, T3Phi_, T3Eta_);
	
	Float_t tempFact[3] = {0, 0, 0};
	Float_t tempCorr[3] = {0, 0, 0};
	
	tempFact[PuCalo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[PuCalo], sType);
	tempCorr[PuCalo] = trkPt_[trkEntry]*tempFact[PuCalo];
	
	tempFact[VsCalo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[VsCalo], sType);
	tempCorr[VsCalo] = trkPt_[trkEntry]*tempFact[VsCalo];
	
	if(montecarlo){
	  tempFact[T] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[T], sType);
	  tempCorr[T] = trkPt_[trkEntry]*tempFact[T];
	}

	Int_t ptIter = getPtRange(trkPt_[trkEntry]);
	
	for(Int_t setIter = 0; setIter < 3; setIter++){
	  if(eventSet_[setIter]){
	    
	    if(getAbsDphi(AlgJtAvePhi_[setIter], trkPhi_[trkEntry]) < TMath::Pi()/2)
	      AlgJtMult_[setIter][0] += tempFact[setIter];
	    else
	      AlgJtMult_[setIter][1] += tempFact[setIter];	  
	    
	    rAlgImbProjA_[setIter + 3][5] += - tempCorr[setIter]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[setIter]));
	    rAlgImbProjA_[setIter + 3][ptIter] += - tempCorr[setIter]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[setIter]));
	    
	    Float_t tempLeadR = getDR(trkEta_[trkEntry], trkPhi_[trkEntry], AlgJtEta_[setIter][0], AlgJtPhi_[setIter][0]);
	    Float_t tempSubLeadR = getDR(trkEta_[trkEntry], trkPhi_[trkEntry], AlgJtEta_[setIter][1], AlgJtPhi_[setIter][1]);
	    
	    if(tempLeadR > 0 && tempSubLeadR > 0){
	      if(tempLeadR < .8 || tempSubLeadR < .8){
		rAlgImbProjAC_[setIter + 3][5] += - tempCorr[setIter]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[setIter]));
		rAlgImbProjAC_[setIter + 3][ptIter] += - tempCorr[setIter]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[setIter]));
	      }
	      else{
		rAlgImbProjANC_[setIter + 3][5] += - tempCorr[setIter]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[setIter]));
		rAlgImbProjANC_[setIter + 3][ptIter] += - tempCorr[setIter]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[setIter]));
	      }
	      
	      for(Int_t rIter = 0; rIter < 10; rIter++){
		if(tempLeadR < rBounds[rIter] || tempSubLeadR < rBounds[rIter]){
		  rAlgImbProjAR_[setIter + 3][5][rIter] += - tempCorr[setIter]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[setIter]));
		  rAlgImbProjAR_[setIter + 3][ptIter][rIter] += - tempCorr[setIter]*cos(getDPHI(trkPhi_[trkEntry], AlgJtAvePhi_[setIter]));
		  break;
		}
	      }
	    }
	  }
	}  
      }
      
      
      if(montecarlo){
	//Iterate over truth
	
	for(Int_t genEntry = 0; genEntry < nGen_; genEntry++){            
	  
	  Int_t ptIter = getPtRange(genPt_[genEntry]);
	  
	  for(Int_t setIter = 0; setIter < 3; setIter++){
	    if(eventSet_[setIter]){
	      
	      gAlgImbProjA_[setIter][5] += -genPt_[genEntry]*cos(getDPHI(genPhi_[genEntry], AlgJtAvePhi_[setIter]));
	      gAlgImbProjA_[setIter][ptIter] += -genPt_[genEntry]*cos(getDPHI(genPhi_[genEntry], AlgJtAvePhi_[setIter]));
	      
	      Float_t tempLeadDelR = getDR(genEta_[genEntry], genPhi_[genEntry], AlgJtEta_[setIter][0], AlgJtPhi_[setIter][0]);
	      Float_t tempSubLeadDelR = getDR(genEta_[genEntry], genPhi_[genEntry], AlgJtEta_[setIter][1], AlgJtPhi_[setIter][1]);
	      
	      if(tempLeadDelR > 0 && tempSubLeadDelR > 0){
		
		if(tempLeadDelR < .8 || tempSubLeadDelR < .8){
		  gAlgImbProjAC_[setIter][5] += -genPt_[genEntry]*cos(getDPHI(genPhi_[genEntry], AlgJtAvePhi_[setIter]));
		  gAlgImbProjAC_[setIter][ptIter] += -genPt_[genEntry]*cos(getDPHI(genPhi_[genEntry], AlgJtAvePhi_[setIter]));
		}
		else{
		  gAlgImbProjANC_[setIter][5] += -genPt_[genEntry]*cos(getDPHI(genPhi_[genEntry], AlgJtAvePhi_[setIter]));
		  gAlgImbProjANC_[setIter][ptIter] += -genPt_[genEntry]*cos(getDPHI(genPhi_[genEntry], AlgJtAvePhi_[setIter]));
		}
		
		for(Int_t rIter = 0; rIter < 10; rIter++){
		  if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
		    gAlgImbProjAR_[setIter][5][rIter] += -genPt_[genEntry]*cos(getDPHI(genPhi_[genEntry], AlgJtAvePhi_[setIter]));
		    gAlgImbProjAR_[setIter][ptIter][rIter] += -genPt_[genEntry]*cos(getDPHI(genPhi_[genEntry], AlgJtAvePhi_[setIter]));
		    break;
		  }
		}	
	      }
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

  centHistFile_80_p->Close();
  delete centHistFile_80_p;

  centHistFile_Merge_p->Close();
  delete centHistFile_Merge_p;

  CleanupDiJetAnaSkim(montecarlo);

  iniSkim_p->Close();
  delete iniSkim_p;

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
