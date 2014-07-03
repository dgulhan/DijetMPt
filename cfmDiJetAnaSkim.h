//=============================================                                 
// Author: Chris McGinn                                                         
//                                                                              
// DiJet Analysis Skim Class (MC)                                                    
//                                                                              
//=============================================  
#ifndef cfmDiJetAnaSkim_h
#define cfmDiJetAnaSkim_h

#include "/net/hisrv0001/home/cfmcginn/DijetMPt/CMSSW_5_3_12_patch3/src/DijetInitialSkim/cfmDiJetIniSkim.h"

TTree* trackTreeAna_p = 0;
TTree* jetTreeAna_p = 0;
TTree* genTreeAna_p = 0;

//Track Tree Variables

//Tracks proj. onto Alg (enum ordered above, w/ corrected in back 5), All, Cone, and NotCone

//2nd Dim 0 == 0_1, 1 == 1_2 ... 5 == F


Float_t rAlgImbProjA_[6][6];
Float_t rAlgImbProjAC_[6][6];
Float_t rAlgImbProjANC_[6][6];

//DelR ProjA

Float_t rAlgImbProjAR_[6][6][10];

//Jet Tree Variables

Int_t run_;
Int_t evt_;
Int_t lumi_;
Int_t hiBin_;

Float_t pthat_;

Float_t hiEvtPlane_;
Float_t psin_;

//Set Bool

//Event Set Bool array, [0] == PuPF, [1] == PuCalo, .etc according to enum

Bool_t eventSet_[5];
Float_t centWeight_80_[5];
Float_t centWeight_Merge_[5];

Bool_t isQuarkJet_[5];
Bool_t isGluonJet_[5];

Float_t pthatWeight_;

//Jet Set, Array by algorithm, according to enum above, 2nd, 0 = lead, 1 =sublead etc.

Float_t AlgJtPt_[5][4];
Float_t AlgJtPhi_[5][4];
Float_t AlgJtEta_[5][4];
Float_t AlgJtTrkMax_[5][4];
Float_t AlgJtRawPt_[5][4];

Float_t AlgJtMult_[5][2];

Float_t AlgJtAvePhi_[5];
Float_t AlgJtDelPhi_[5];
Float_t AlgJtAsymm_[5];

Float_t AlgRefPt_[5][4];
Float_t AlgRefPhi_[5][4];
Float_t AlgRefEta_[5][4];

Float_t AlgRefAvePhi_[5];
Float_t AlgRefDelPhi_[5];
Float_t AlgRefAsymm_[5];

//Gen Tree Variables

//Gen. proj. onto Jets, ordered by algorithm according to enum, PuPF == [0], PuCalo == [1], etc.

Float_t gAlgImbProjA_[3][6];
Float_t gAlgImbProjAC_[3][6];
Float_t gAlgImbProjANC_[3][6];

//truth delRs
//ProjA delR

Float_t gAlgImbProjAR_[3][6][10];

void SetAnaBranches(sampleType sType = kHIDATA)
{
  Bool_t montecarlo = false;
  if(sType == kHIMC || sType == kPPMC || sType == kPAMC) montecarlo = true;

  //Track Tree Branches

  std::cout << "Branches Set" << std::endl;

  //Tracks proj. onto Alg, ordered according to enum above, All, Cone, and NotCone

  trackTreeAna_p->Branch("rAlgImbProjA", &rAlgImbProjA_, "rAlgImbProjA[6][6]/F");
  trackTreeAna_p->Branch("rAlgImbProjAC", &rAlgImbProjAC_, "rAlgImbProjAC[6][6]/F");
  trackTreeAna_p->Branch("rAlgImbProjANC", &rAlgImbProjANC_, "rAlgImbProjANC[6][6]/F");

  //ProjA DelRs

  trackTreeAna_p->Branch("rAlgImbProjAR", &rAlgImbProjAR_, "rAlgImbProjAR[6][6][10]/F");

  //Jet Tree Branches

  jetTreeAna_p->Branch("run", &run_, "run/I");
  jetTreeAna_p->Branch("evt", &evt_, "evt/I");
  jetTreeAna_p->Branch("lumi", &lumi_, "lumi/I");

  if(sType == kHIDATA || sType == kHIMC){
    jetTreeAna_p->Branch("hiBin", &hiBin_, "hiBin/I");
    jetTreeAna_p->Branch("hiEvtPlane", &hiEvtPlane_, "hiEvtPlane/F");
    jetTreeAna_p->Branch("psin", &psin_, "psin/F");
  }

  jetTreeAna_p->Branch("eventSet", &eventSet_, "eventSet[5]/O");

  if(sType == kHIMC){
    jetTreeAna_p->Branch("centWeight_80", &centWeight_80_, "centWeight_80[5]/F");
    jetTreeAna_p->Branch("centWeight_Merge", &centWeight_Merge_, "centWeight_Merge[5]/F");
  }

  jetTreeAna_p->Branch("AlgJtPt", &AlgJtPt_, "AlgJtPt[5][4]/F");
  jetTreeAna_p->Branch("AlgJtPhi", &AlgJtPhi_, "AlgJtPhi[5][4]/F");
  jetTreeAna_p->Branch("AlgJtEta", &AlgJtEta_, "AlgJtEta[5][4]/F");
  jetTreeAna_p->Branch("AlgJtTrkMax", &AlgJtTrkMax_, "AlgJtTrkMax[5][4]/F");
  jetTreeAna_p->Branch("AlgJtRawPt", &AlgJtRawPt_, "AlgJtRawPt[5][4]/F");

  jetTreeAna_p->Branch("AlgJtMult", &AlgJtMult_, "AlgJtMult[5][2]/F");

  jetTreeAna_p->Branch("AlgJtAvePhi", &AlgJtAvePhi_, "AlgJtAvePhi[5]/F");
  jetTreeAna_p->Branch("AlgJtDelPhi", &AlgJtDelPhi_, "AlgJtDelPhi[5]/F");
  jetTreeAna_p->Branch("AlgJtAsymm", &AlgJtAsymm_, "AlgJtAsymm[5]/F");

  if(montecarlo){
    jetTreeAna_p->Branch("isQuarkJet", &isQuarkJet_, "isQuarkJet[5]/O");
    jetTreeAna_p->Branch("isGluonJet", &isGluonJet_, "isGluonJet[5]/O");
  }

  if(sType == kHIMC)
    jetTreeAna_p->Branch("pthatWeight", &pthatWeight_, "pthatWeight/F");

  if(montecarlo){
    jetTreeAna_p->Branch("pthat", &pthat_, "pthat/F");

    //refpt for jets immediately above
    jetTreeAna_p->Branch("AlgRefPt", &AlgRefPt_, "AlgRefPt[5][4]/F");
    jetTreeAna_p->Branch("AlgRefPhi", &AlgRefPhi_, "AlgRefPhi[5][4]/F");
    jetTreeAna_p->Branch("AlgRefEta", &AlgRefEta_, "AlgRefEta[5][4]/F");

    jetTreeAna_p->Branch("AlgRefAvePhi", &AlgRefAvePhi_, "AlgRefAvePhi[5]/F");
    jetTreeAna_p->Branch("AlgRefDelPhi", &AlgRefDelPhi_, "AlgRefDelPhi[5]/F");
    jetTreeAna_p->Branch("AlgRefAsymm", &AlgRefAsymm_, "AlgRefAsymm[5]/F");


    //Gen. proj. onto jetAlg, array ordered according to enum

    genTreeAna_p->Branch("gAlgImbProjA", &gAlgImbProjA_, "gAlgImbProjA[3][6]/F");
    genTreeAna_p->Branch("gAlgImbProjAC", &gAlgImbProjAC_, "gAlgImbProjAC[3][6]/F");
    genTreeAna_p->Branch("gAlgImbProjANC", &gAlgImbProjANC_, "gAlgImbProjANC[3][6]/F");

    //Truth Del Rs, Proj
    //Proj A's

    genTreeAna_p->Branch("gAlgImbProjAR", &gAlgImbProjAR_, "gAlgImbProjAR[3][6][10]/F");
  }
}


void GetAnaBranches(sampleType sType = kHIDATA){
  Bool_t montecarlo = false;
  if(sType == kHIMC || sType == kPPMC || sType == kPAMC) montecarlo = true;

  std::cout << "Get Branches" << std::endl;

  //Track Tree Branches
   trackTreeAna_p->SetBranchAddress("rAlgImbProjA", rAlgImbProjA_);
   trackTreeAna_p->SetBranchAddress("rAlgImbProjAC", rAlgImbProjAC_);
   trackTreeAna_p->SetBranchAddress("rAlgImbProjANC", rAlgImbProjANC_);
   trackTreeAna_p->SetBranchAddress("rAlgImbProjAR", rAlgImbProjAR_);

   //Jet Tree Branches

   jetTreeAna_p->SetBranchAddress("run", &run_);
   jetTreeAna_p->SetBranchAddress("evt", &evt_);
   jetTreeAna_p->SetBranchAddress("lumi", &lumi_);

   if(sType == kHIDATA || sType == kHIMC)
     jetTreeAna_p->SetBranchAddress("hiBin", &hiBin_);

   jetTreeAna_p->SetBranchAddress("pthat", &pthat_);

   if(sType == kHIDATA || sType == kHIMC){
     jetTreeAna_p->SetBranchAddress("hiEvtPlane", &hiEvtPlane_);
     jetTreeAna_p->SetBranchAddress("psin", &psin_);
   }

   jetTreeAna_p->SetBranchAddress("eventSet", eventSet_);

   if(sType == kHIDATA){
     jetTreeAna_p->SetBranchAddress("centWeight_80", centWeight_80_);
     jetTreeAna_p->SetBranchAddress("centWeight_Merge", centWeight_Merge_);
   }

   if(montecarlo){
     jetTreeAna_p->SetBranchAddress("isQuarkJet", isQuarkJet_);
     jetTreeAna_p->SetBranchAddress("isGluonJet", isGluonJet_);
   }     

   if(sType == kHIDATA)
     jetTreeAna_p->SetBranchAddress("pthatWeight", &pthatWeight_);

   jetTreeAna_p->SetBranchAddress("AlgJtPt", AlgJtPt_);
   jetTreeAna_p->SetBranchAddress("AlgJtPhi", AlgJtPhi_);
   jetTreeAna_p->SetBranchAddress("AlgJtEta", AlgJtEta_);
   jetTreeAna_p->SetBranchAddress("AlgJtTrkMax", AlgJtTrkMax_);
   jetTreeAna_p->SetBranchAddress("AlgJtRawPt", AlgJtRawPt_);

   jetTreeAna_p->SetBranchAddress("AlgJtMult", AlgJtMult_);

   jetTreeAna_p->SetBranchAddress("AlgJtAvePhi", AlgJtAvePhi_);
   jetTreeAna_p->SetBranchAddress("AlgJtDelPhi", AlgJtDelPhi_);
   jetTreeAna_p->SetBranchAddress("AlgJtAsymm", AlgJtAsymm_);


   if(montecarlo){
     jetTreeAna_p->SetBranchAddress("AlgRefPt", AlgRefPt_);
     jetTreeAna_p->SetBranchAddress("AlgRefPhi", AlgRefPhi_);
     jetTreeAna_p->SetBranchAddress("AlgRefEta", AlgRefEta_);
     
     jetTreeAna_p->SetBranchAddress("AlgRefAvePhi", AlgRefAvePhi_);
     jetTreeAna_p->SetBranchAddress("AlgRefDelPhi", AlgRefDelPhi_);
     jetTreeAna_p->SetBranchAddress("AlgRefAsymm", AlgRefAsymm_);

     //Gen Tree Variables

     genTreeAna_p->SetBranchAddress("gAlgImbProjA", gAlgImbProjA_);     
     genTreeAna_p->SetBranchAddress("gAlgImbProjAC", gAlgImbProjAC_);     
     genTreeAna_p->SetBranchAddress("gAlgImbProjANC", gAlgImbProjANC_);     

     genTreeAna_p->SetBranchAddress("gAlgImbProjAR", gAlgImbProjAR_);     
   }
}


void InitDiJetAnaSkim(sampleType sType = kHIDATA)
{
  std::cout << "Init DiJet AnaSkim" << std::endl;

  trackTreeAna_p = new TTree("trackTreeAna", "trackTreeAna");
  jetTreeAna_p = new TTree("jetTreeAna", "jetTreeAna");

  if(sType == kHIMC || sType == kPPMC || sType == kPAMC)
    genTreeAna_p = new TTree("genTreeAna", "genTreeAna");

  SetAnaBranches(sType);
}


void CleanupDiJetAnaSkim()
{
  if(trackTreeAna_p != 0){
    delete trackTreeAna_p;
    trackTreeAna_p = 0;
  }
  if(jetTreeAna_p != 0){
    delete jetTreeAna_p;
    jetTreeAna_p = 0;
  }
  if(genTreeAna_p != 0){
    delete genTreeAna_p;
    genTreeAna_p = 0;
  }
}


void GetDiJetAnaSkim(TFile* anaFile_p, sampleType sType = kHIDATA){
  std::cout << "Get DiJet AnaSkim" << std::endl;

  trackTreeAna_p = (TTree*)anaFile_p->Get("trackTreeAna");
  jetTreeAna_p = (TTree*)anaFile_p->Get("jetTreeAna");

  if(sType == kHIMC || sType == kPPMC || sType == kPAMC)
    genTreeIni_p = (TTree*)anaFile_p->Get("genTreeAna");

  GetAnaBranches(sType);
}


void InitJetVar(sampleType sType = kHIDATA)
{
  Bool_t montecarlo = false;
  if(sType == kHIMC || sType == kPPMC || sType == kPAMC) montecarlo = true;

  for(Int_t initIter = 0; initIter < 5; initIter++){
    eventSet_[initIter] = false;

    if(sType == kHIMC){
      centWeight_80_[initIter] = 1;    
      centWeight_Merge_[initIter] = 1;    
    }

    for(Int_t initIter2 = 0; initIter2 < 4; initIter2++){
      AlgJtPt_[initIter][initIter2] = -10;
      AlgJtPhi_[initIter][initIter2] = -10;
      AlgJtEta_[initIter][initIter2] = -10;
      AlgJtTrkMax_[initIter][initIter2] = -10;
      AlgJtRawPt_[initIter][initIter2] = -10;

      if(initIter2 < 2)
	AlgJtMult_[initIter][initIter2] = 0;
    }

    AlgJtAvePhi_[initIter] = -10;
    AlgJtDelPhi_[initIter] = -10;
    AlgJtAsymm_[initIter] = -10;

    if(montecarlo){
      isQuarkJet_[initIter] = false;
      isGluonJet_[initIter] = false;

      pthatWeight_ = -10;

      for(Int_t initIter2 = 0; initIter2 < 4; initIter2++){
	AlgRefPt_[initIter][initIter2] = -10;
	AlgRefPhi_[initIter][initIter2] = -10;
	AlgRefEta_[initIter][initIter2] = -10;
      }

      AlgRefAvePhi_[initIter] = -10;
      AlgRefDelPhi_[initIter] = -10;
      AlgRefAsymm_[initIter] = -10;
    }
  }
}

void InitProjPerp(Bool_t montecarlo = false)
{
  //Tracks proj. onto Alg, ordered according to enum above, corr in the back 5, All, Cone, and NotCone

  for(Int_t initIter = 0; initIter < 6; initIter++){
    for(Int_t initIter2 = 0; initIter2 < 6; initIter2++){
      rAlgImbProjA_[initIter][initIter2] = 0;
      rAlgImbProjAC_[initIter][initIter2] = 0;
      rAlgImbProjANC_[initIter][initIter2] = 0;

      for(Int_t initIter3 = 0; initIter3 < 10; initIter3++){
	rAlgImbProjAR_[initIter][initIter2][initIter3] = 0;
      }
    }
  }

  if(montecarlo){
    //Gen. proj. onto Truth
    for(Int_t initIter = 0; initIter < 3; initIter++){
      for(Int_t initIter2 = 0; initIter2 < 6; initIter2++){
	gAlgImbProjA_[initIter][initIter2] = 0;
	gAlgImbProjAC_[initIter][initIter2] = 0;
	gAlgImbProjANC_[initIter][initIter2] = 0;

	for(Int_t initIter3 = 0; initIter3 < 10; initIter3++){
	  gAlgImbProjAR_[initIter][initIter2][initIter3] = 0;
	}
      }
    }    
  }
}


#endif
