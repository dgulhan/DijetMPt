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
#include "TComplex.h"


int makeDiJetAnaSkim(std::string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_DIJETANASKIM.root", Int_t num = 0, Bool_t justJt = false)
{
  //Define MC or Data
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

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
    while(!inFile.eof()){
      inFile >> buffer;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }

  std::cout << "FileList Loaded" << std::endl;

  TFile* iniSkim_p = new TFile(listOfFiles[0].data(), "READ");

  GetDiJetIniSkim(iniSkim_p, sType, justJt);

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

  InitDiJetAnaSkim(sType, justJt);

  Long64_t nentries = jetTreeIni_p->GetEntries();

  std::cout << nentries << std::endl;

  std::cout << "Cuts, Lead/Sublead Pt, delphi, eta: " << leadJtPtCut << ", " << subLeadJtPtCut << ", " << jtDelPhiCut << ", " << jtEtaCut << std::endl; 

  for(Long64_t jentry = 0; jentry < nentries; jentry++){
    jetTreeIni_p->GetEntry(jentry);

    if(!justJt){
      trackTreeIni_p->GetEntry(jentry);
      
      if(montecarlo)
	genTreeIni_p->GetEntry(jentry);
    }      

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
    
    if(hi){
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

    InitProjPerp(sType);

    //Switch below to iterated OR EDIT HERE

    if((eventSet_[PuCalo] || eventSet_[VsCalo] || eventSet_[T]) && !justJt){
      
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
      
      if(hi)
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

    if(!justJt){
      trackTreeAna_p->Fill();
    
      if(montecarlo) genTreeAna_p->Fill();
    }
  }
  
  outFile->cd();
  jetTreeAna_p->Write("", TObject::kOverwrite);

  if(!justJt){
    trackTreeAna_p->Write("", TObject::kOverwrite);

    if(montecarlo) genTreeAna_p->Write("", TObject::kOverwrite);
  }

  centHistFile_80_p->Close();
  delete centHistFile_80_p;

  centHistFile_Merge_p->Close();
  delete centHistFile_Merge_p;

  CleanupDiJetAnaSkim();
  outFile->Close();
  delete outFile;

  iniSkim_p->Close();
  delete iniSkim_p;

  printf("Done.\n");
  return(0);
}


int main(int argc, char *argv[])
{
  if(argc != 6)
    {
      std::cout << "Usage: sortForest <inputFile> <sampleType> <outputFile> <#> <justJtBool>" << std::endl;
      return 1;
    }

  int rStatus = -1;

  rStatus = makeDiJetAnaSkim(argv[1], sampleType(atoi(argv[2])), argv[3], atoi(argv[4]), Bool_t(atoi(argv[5])));

  return rStatus;
}
