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
#include "TComplex.h"


Int_t pthatCuts_PYTH_HITrk[10] = {15, 30, 50, 80, 120, 170, 220, 280, 370, 10000000};
Float_t pthatWeights_PYTH_HITrk[9] = {.551019, .034814, .00242254, .000304825, .0000426788, .00000492814, .000000879673, .00000017353, .0000000292439};


int makeDiJetAnaSkim(std::string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_DIJETANASKIM.root", Int_t num = 0, Bool_t justJt = false, Bool_t isHITrk = false)
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

  TFile* iniSkim_p = new TFile(listOfFiles[num].data(), "READ");

  GetDiJetIniSkim(iniSkim_p, sType, justJt);

  std::cout << "IniSkim Loaded" << std::endl;

  //Setup correction tables

  InitCorrFiles(sType, isHITrk);
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

    if(montecarlo && isHITrk){
      for(Int_t hatIter = 0; hatIter < 9; hatIter++){
	if(pthat_ > pthatCuts_PYTH_HITrk[hatIter] && pthat_ < pthatCuts_PYTH_HITrk[hatIter + 1]){
	  pthatWeight_ = pthatWeights_PYTH_HITrk[hatIter];
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
      
      
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
        
	//Grab proj. Pt Spectra For Tracks in each Event Subset
	
	for(Int_t jtIter = 0; jtIter < 3; jtIter++){
	  if(eventSet_[jtIter])
	    GetTrkProjPerp(jtIter, jtIter, trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry]);
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

	Float_t tempFact[3] = {0., 0., 0.};
	Float_t tempCorr[3] = {0., 0., 0.};

	tempFact[PuCalo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[PuCalo], sType);
	tempCorr[PuCalo] = trkPt_[trkEntry]*tempFact[PuCalo];

	tempFact[VsCalo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[VsCalo], sType);
	tempCorr[VsCalo] = trkPt_[trkEntry]*tempFact[VsCalo];

	if(montecarlo){
	  tempFact[T] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[T], sType);
	  tempCorr[T] = trkPt_[trkEntry]*tempFact[T];
	}

	for(Int_t jtIter = 0; jtIter < 3; jtIter++){
	  if(eventSet_[jtIter]) 
            GetTrkProjPerp(jtIter, jtIter+3, trkPt_[trkEntry], tempCorr[jtIter], trkPhi_[trkEntry], trkEta_[trkEntry]);
	}

      }

      if(montecarlo){
	//Iterate over Truth

	for(Int_t genEntry = 0; genEntry < nGen_; genEntry++){
	  for(Int_t jtIter = 0; jtIter < 3; jtIter++){
	    GetGenProjPerp(jtIter, genPt_[genEntry], genPhi_[genEntry], genEta_[genEntry]);
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
  if(argc != 7)
    {
      std::cout << "Usage: sortForest <inputFile> <sampleType> <outputFile> <#> <justJtBool> <isHITrk>" << std::endl;
      return 1;
    }

  int rStatus = -1;

  rStatus = makeDiJetAnaSkim(argv[1], sampleType(atoi(argv[2])), argv[3], atoi(argv[4]), Bool_t(atoi(argv[5])), Bool_t(atoi(argv[6])));

  return rStatus;
}
