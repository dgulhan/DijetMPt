//=============================================
// Author: Chris McGinn
// 
// DiJet Analysis Skim Class (MC); Mixing mod
//
//=============================================

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "cfmDiJetAnaSkim.h"
#include "stdlib.h"
#include <fstream>
#include "TComplex.h"
#include "cfmTreeCentSort.h"

int makeDiJetAnaSkim_MIXING(std::string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_DIJETANASKIM.root", Int_t num = 0, Bool_t justJt = false)
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

  TFile* mixFile_p = new TFile("/mnt/hadoop/cms/store/user/dgulhan/HIMC/MB/Track8_Jet26_STARTHI53_LV1/merged2/HiForest_HYDJET_Track8_Jet26_STARTHI53_LV1_merged_forest_0.root", "READ");
  TTree* mixTree_p = (TTree*)mixFile_p->Get("hiEvtAnalyzer/HiTree");
  mixTree_p->AddFriend("skimanalysis/HltTree");
  RunFillCentEntryVect(mixTree_p);

  TTree* mixGenTree_p = (TTree*)mixFile_p->Get("HiGenParticleAna/hi");

  Int_t tempNGen_ = 0;
  Float_t tempGenPt_[maxEntrySim];
  Float_t tempGenPhi_[maxEntrySim];
  Float_t tempGenEta_[maxEntrySim];
  Int_t tempGenChg_[maxEntrySim];

  mixGenTree_p->SetBranchStatus("*", 0);
  mixGenTree_p->SetBranchStatus("mult", 1);
  mixGenTree_p->SetBranchStatus("pt", 1);
  mixGenTree_p->SetBranchStatus("phi", 1);
  mixGenTree_p->SetBranchStatus("eta", 1);
  mixGenTree_p->SetBranchStatus("chg", 1);

  mixGenTree_p->SetBranchAddress("mult", &tempNGen_);
  mixGenTree_p->SetBranchAddress("pt", tempGenPt_);
  mixGenTree_p->SetBranchAddress("phi", tempGenPhi_);
  mixGenTree_p->SetBranchAddress("eta", tempGenEta_);
  mixGenTree_p->SetBranchAddress("chg", tempGenChg_);


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

  
    //Iterate over tracks

    InitProjPerp(sType);

    //Switch below to iterated OR EDIT HERE

    if(eventSet_[VsCalo] && centEntryVect_p[hiBin_]->size() != 0){
      eventPt_[0] = 0;
      eventPt_[1] = 0;

      GetMixProjPerp(eventPt_, nLeadJtConst_, TrkLeadJtConstPt_, TrkLeadJtConstPhi_, TrkLeadJtConstEta_, TrkLeadJtConstCorr_);
      GetMixProjPerp(eventPt_, nSubLeadJtConst_, TrkSubLeadJtConstPt_, TrkSubLeadJtConstPhi_, TrkSubLeadJtConstEta_, TrkSubLeadJtConstCorr_);
      GetMixProjPerp(eventPt_, nThirdJtConst_, TrkThirdJtConstPt_, TrkThirdJtConstPhi_, TrkThirdJtConstEta_, TrkThirdJtConstCorr_);
      //      GetMixProjPerp(eventPt_, nFourthJtConst_, TrkFourthJtConstPt_, TrkFourthJtConstPhi_, TrkFourthJtConstEta_, TrkFourthJtConstCorr_);
      //      GetMixProjPerp(eventPt_, nFifthJtConst_, TrkFifthJtConstPt_, TrkFifthJtConstPhi_, TrkFifthJtConstEta_, TrkFifthJtConstCorr_);

      mixGenTree_p->GetEntry(centEntryVect_p[hiBin_]->at(0));
      mixTree_p->GetEntry(centEntryVect_p[hiBin_]->at(0));
      centEntryVect_p[hiBin_]->erase(centEntryVect_p[hiBin_]->begin());

      Int_t chgMult = 0;
      for(Int_t genIter = 0; genIter < tempNGen_; genIter++){
	if(tempGenChg_[genIter] == 0 || TMath::Abs(tempGenEta_[genIter]) > 2.4) continue;

        GetTrkProjPerp(1, 1, tempGenPt_[genIter], tempGenPt_[genIter], tempGenPhi_[genIter], tempGenEta_[genIter]);
	chgMult++;
      }

      eventPt_[0] = eventPt_[0]/chgMult;
      eventPt_[1] = eventPt_[1]/chgMult;

      for(Int_t genIter = 0; genIter < tempNGen_; genIter++){
	if(tempGenChg_[genIter] == 0 || TMath::Abs(tempGenEta_[genIter]) > 2.4) continue;

	Float_t tempGenPtCorr[2] = {tempGenPt_[genIter]*cos(tempGenPhi_[genIter]) + eventPt_[0], tempGenPt_[genIter]*sin(tempGenPhi_[genIter]) + eventPt_[1]};

	GetTrkProjPerp(1, 4, TMath::Sqrt(tempGenPtCorr[0]*tempGenPtCorr[0] + tempGenPtCorr[1]*tempGenPtCorr[1]), TMath::Sqrt(tempGenPtCorr[0]*tempGenPtCorr[0] + tempGenPtCorr[1]*tempGenPtCorr[1]), tempGenPhi_[genIter], tempGenEta_[genIter]);
      }

    }
    else if(centEntryVect_p[hiBin_]->size() == 0){
      std::cout << "Ran out of minBias events at centrality: " << hiBin_ << std::endl;
      continue;
    }

    jetTreeAna_p->Fill();
    trackTreeAna_p->Fill();
    if(montecarlo) genTreeAna_p->Fill();
  }
  
  outFile->cd();
  jetTreeAna_p->Write("", TObject::kOverwrite);
  trackTreeAna_p->Write("", TObject::kOverwrite);
  if(montecarlo) genTreeAna_p->Write("", TObject::kOverwrite);

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

  rStatus = makeDiJetAnaSkim_MIXING(argv[1], sampleType(atoi(argv[2])), argv[3], atoi(argv[4]), Bool_t(atoi(argv[5])));

  return rStatus;
}
