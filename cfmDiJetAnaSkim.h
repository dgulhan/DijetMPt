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


//Cut var
const Float_t leadJtPtCut = 120.;
const Float_t subLeadJtPtCut = 50.;
const Float_t jtDelPhiCut = 0;
const Float_t jtEtaCut = 2.0; // Default Max at 2.4 to avoid transition junk, otherwise vary as needed                                                 
const char* algType[5] = {"PuCalo", "VsCalo", "T", "PuPF", "VsPF"};


void SetAnaBranches(sampleType sType = kHIDATA, Bool_t justJt = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

  //Track Tree Branches

  std::cout << "Branches Set" << std::endl;

  //Tracks proj. onto Alg, ordered according to enum above, All, Cone, and NotCone

  if(!justJt){
    trackTreeAna_p->Branch("rAlgImbProjA", rAlgImbProjA_, "rAlgImbProjA[6][6]/F");
    trackTreeAna_p->Branch("rAlgImbProjAC", rAlgImbProjAC_, "rAlgImbProjAC[6][6]/F");
    trackTreeAna_p->Branch("rAlgImbProjANC", rAlgImbProjANC_, "rAlgImbProjANC[6][6]/F");
    
    //ProjA DelRs
    
    trackTreeAna_p->Branch("rAlgImbProjAR", rAlgImbProjAR_, "rAlgImbProjAR[6][6][10]/F");
  }

  //Jet Tree Branches

  jetTreeAna_p->Branch("run", &run_, "run/I");
  jetTreeAna_p->Branch("evt", &evt_, "evt/I");
  jetTreeAna_p->Branch("lumi", &lumi_, "lumi/I");

  if(hi){
    jetTreeAna_p->Branch("hiBin", &hiBin_, "hiBin/I");
    jetTreeAna_p->Branch("hiEvtPlane", &hiEvtPlane_, "hiEvtPlane/F");
    jetTreeAna_p->Branch("psin", &psin_, "psin/F");
  }

  jetTreeAna_p->Branch("eventSet", eventSet_, "eventSet[5]/O");

  if(sType == kHIMC){
    jetTreeAna_p->Branch("centWeight_80", centWeight_80_, "centWeight_80[5]/F");
    jetTreeAna_p->Branch("centWeight_Merge", centWeight_Merge_, "centWeight_Merge[5]/F");
  }

  jetTreeAna_p->Branch("AlgJtPt", AlgJtPt_, "AlgJtPt[5][4]/F");
  jetTreeAna_p->Branch("AlgJtPhi", AlgJtPhi_, "AlgJtPhi[5][4]/F");
  jetTreeAna_p->Branch("AlgJtEta", AlgJtEta_, "AlgJtEta[5][4]/F");
  jetTreeAna_p->Branch("AlgJtTrkMax", AlgJtTrkMax_, "AlgJtTrkMax[5][4]/F");
  jetTreeAna_p->Branch("AlgJtRawPt", AlgJtRawPt_, "AlgJtRawPt[5][4]/F");

  jetTreeAna_p->Branch("AlgJtMult", AlgJtMult_, "AlgJtMult[5][2]/F");

  jetTreeAna_p->Branch("AlgJtAvePhi", AlgJtAvePhi_, "AlgJtAvePhi[5]/F");
  jetTreeAna_p->Branch("AlgJtDelPhi", AlgJtDelPhi_, "AlgJtDelPhi[5]/F");
  jetTreeAna_p->Branch("AlgJtAsymm", AlgJtAsymm_, "AlgJtAsymm[5]/F");

  if(montecarlo){
    jetTreeAna_p->Branch("isQuarkJet", isQuarkJet_, "isQuarkJet[5]/O");
    jetTreeAna_p->Branch("isGluonJet", isGluonJet_, "isGluonJet[5]/O");
  }

  if(sType == kHIMC)
    jetTreeAna_p->Branch("pthatWeight", &pthatWeight_, "pthatWeight/F");

  if(montecarlo){
    jetTreeAna_p->Branch("pthat", &pthat_, "pthat/F");

    //refpt for jets immediately above
    jetTreeAna_p->Branch("AlgRefPt", AlgRefPt_, "AlgRefPt[5][4]/F");
    jetTreeAna_p->Branch("AlgRefPhi", AlgRefPhi_, "AlgRefPhi[5][4]/F");
    jetTreeAna_p->Branch("AlgRefEta", AlgRefEta_, "AlgRefEta[5][4]/F");

    jetTreeAna_p->Branch("AlgRefAvePhi", AlgRefAvePhi_, "AlgRefAvePhi[5]/F");
    jetTreeAna_p->Branch("AlgRefDelPhi", AlgRefDelPhi_, "AlgRefDelPhi[5]/F");
    jetTreeAna_p->Branch("AlgRefAsymm", AlgRefAsymm_, "AlgRefAsymm[5]/F");


    //Gen. proj. onto jetAlg, array ordered according to enum

    if(!justJt){
      genTreeAna_p->Branch("gAlgImbProjA", gAlgImbProjA_, "gAlgImbProjA[3][6]/F");
      genTreeAna_p->Branch("gAlgImbProjAC", gAlgImbProjAC_, "gAlgImbProjAC[3][6]/F");
      genTreeAna_p->Branch("gAlgImbProjANC", gAlgImbProjANC_, "gAlgImbProjANC[3][6]/F");
      
      //Truth Del Rs, Proj
      //Proj A's
      
      genTreeAna_p->Branch("gAlgImbProjAR", gAlgImbProjAR_, "gAlgImbProjAR[3][6][10]/F");
    }
  }
}


void GetAnaBranches(sampleType sType = kHIDATA, Bool_t justJt = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

  std::cout << "Get Branches" << std::endl;

  //Track Tree Branches
  if(!justJt){
    trackTreeAna_p->SetBranchAddress("rAlgImbProjA", rAlgImbProjA_);
    trackTreeAna_p->SetBranchAddress("rAlgImbProjAC", rAlgImbProjAC_);
    trackTreeAna_p->SetBranchAddress("rAlgImbProjANC", rAlgImbProjANC_);
    trackTreeAna_p->SetBranchAddress("rAlgImbProjAR", rAlgImbProjAR_);
  }    
  
  //Jet Tree Branches

   jetTreeAna_p->SetBranchAddress("run", &run_);
   jetTreeAna_p->SetBranchAddress("evt", &evt_);
   jetTreeAna_p->SetBranchAddress("lumi", &lumi_);

   if(montecarlo)
     jetTreeAna_p->SetBranchAddress("pthat", &pthat_);

   if(hi){
     jetTreeAna_p->SetBranchAddress("hiBin", &hiBin_);
     jetTreeAna_p->SetBranchAddress("hiEvtPlane", &hiEvtPlane_);
     jetTreeAna_p->SetBranchAddress("psin", &psin_);
   }

   jetTreeAna_p->SetBranchAddress("eventSet", eventSet_);

   if(sType == kHIMC){
     jetTreeAna_p->SetBranchAddress("centWeight_80", centWeight_80_);
     jetTreeAna_p->SetBranchAddress("centWeight_Merge", centWeight_Merge_);
   }

   if(montecarlo){
     jetTreeAna_p->SetBranchAddress("isQuarkJet", isQuarkJet_);
     jetTreeAna_p->SetBranchAddress("isGluonJet", isGluonJet_);
   }     

   if(sType == kHIMC)
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

     if(!justJt){
       genTreeAna_p->SetBranchAddress("gAlgImbProjA", gAlgImbProjA_);     
       genTreeAna_p->SetBranchAddress("gAlgImbProjAC", gAlgImbProjAC_);     
       genTreeAna_p->SetBranchAddress("gAlgImbProjANC", gAlgImbProjANC_);     
       
       genTreeAna_p->SetBranchAddress("gAlgImbProjAR", gAlgImbProjAR_);     
     }
   }
}


void InitDiJetAnaSkim(sampleType sType = kHIDATA, Bool_t justJt = false)
{
  std::cout << "Init DiJet AnaSkim" << std::endl;

  if(!justJt) trackTreeAna_p = new TTree("trackTreeAna", "trackTreeAna");
  jetTreeAna_p = new TTree("jetTreeAna", "jetTreeAna");
  if(isMonteCarlo(sType) && !justJt) genTreeAna_p = new TTree("genTreeAna", "genTreeAna");

  SetAnaBranches(sType, justJt);
}


void CleanupDiJetAnaSkim()
{
  if(trackTreeAna_p != 0) delete trackTreeAna_p;
  if(jetTreeAna_p != 0) delete jetTreeAna_p;
  if(genTreeAna_p != 0) delete genTreeAna_p;
}


void GetDiJetAnaSkim(TFile* anaFile_p, sampleType sType = kHIDATA, Bool_t justJt = false){
  std::cout << "Get DiJet AnaSkim" << std::endl;

  if(!justJt) trackTreeAna_p = (TTree*)anaFile_p->Get("trackTreeAna");
  jetTreeAna_p = (TTree*)anaFile_p->Get("jetTreeAna");
  if(isMonteCarlo(sType) && !justJt) genTreeIni_p = (TTree*)anaFile_p->Get("genTreeAna");

  GetAnaBranches(sType, justJt);
}


Float_t getAbsDphi(Float_t phi1, Float_t phi2)
{
  return TMath::Abs(getDPHI(phi1, phi2));
}


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



void InitJetVar(sampleType sType = kHIDATA)
{
  Bool_t montecarlo = isMonteCarlo(sType);

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

void InitProjPerp(sampleType sType = kHIDATA)
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

  if(isMonteCarlo(sType)){
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



#endif
