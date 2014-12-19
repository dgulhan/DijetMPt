//=============================================                                 
// Author: Chris McGinn                                                         
//                                                                              
// DiJet Analysis Skim Class (MC)                                                    
//                                                                              
//=============================================  
#ifndef cfmDiJetAnaSkim_h
#define cfmDiJetAnaSkim_h

#include "/net/hisrv0001/home/cfmcginn/DijetMPt/CMSSW_5_3_20/src/DijetInitialSkim/cfmDiJetIniSkim.h"
#include <string>

TTree* trackTreeAna_p = 0;
TTree* jetTreeAna_p = 0;
TTree* genTreeAna_p = 0;

const Int_t nJtAlg = 18;
const Int_t nJtMax = 4;
const Int_t nSumAlg = 16;
const Int_t nPtBins = 6;
const Int_t nRBins = 10;

//Track Tree Variables

//Tracks proj. onto Alg (enum ordered above, w/ corrected in back 5), All, Cone, and NotCone

//2nd Dim 0 == 0_1, 1 == 1_2 ... 5 == F

Float_t rAlgJtMult_[2*nSumAlg][2];

Float_t rAlgImbProjA_[2*nSumAlg][nPtBins];
Float_t rAlgImbProjAC_[2*nSumAlg][nPtBins];
Float_t rAlgImbProjANC_[2*nSumAlg][nPtBins];

Float_t rAlgImbProjAC0_[2*nSumAlg][nPtBins];
Float_t rAlgImbProjAC1_[2*nSumAlg][nPtBins];
Float_t rAlgImbProjAC2_[2*nSumAlg][nPtBins];
Float_t rAlgImbProjAC3_[2*nSumAlg][nPtBins];

//DelR ProjA

Float_t rAlgImbProjAR_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgImbProjAEta_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgImbProjAPhi_[2*nSumAlg][nPtBins][nRBins];

Float_t rAlgImbProjAR_Cut_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgImbProjAR_CutEta_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgImbProjAR_CutPhi_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgImbProjAEta_Cut_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgImbProjAPhi_Cut_[2*nSumAlg][nPtBins][nRBins];

Float_t rAlgImbProjAR_FOR_[2*nSumAlg][nPtBins][nRBins];

//DelR Mult

Float_t rAlgMultAR_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgMultAEta_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgMultAPhi_[2*nSumAlg][nPtBins][nRBins];

Float_t rAlgMultAR_Cut_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgMultAR_CutEta_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgMultAR_CutPhi_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgMultAEta_Cut_[2*nSumAlg][nPtBins][nRBins];
Float_t rAlgMultAPhi_Cut_[2*nSumAlg][nPtBins][nRBins];

Float_t rAlgImbRawAR_[2*nSumAlg][nPtBins][nRBins];


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

Bool_t eventSet_[nJtAlg];
Float_t centWeight_[nJtAlg];
Float_t centWeight_Merge_[nJtAlg];

Bool_t isQuarkJet_[nJtAlg];
Bool_t isGluonJet_[nJtAlg];

Float_t pthatWeight_;
Float_t leadEtaWeight_;
Float_t subleadEtaWeight_;

//Jet Set, Array by algorithm, according to enum above, 2nd, 0 = lead, 1 =sublead etc.

Float_t AlgJtPt_[nJtAlg][nJtMax];
Float_t AlgJtPhi_[nJtAlg][nJtMax];
Float_t AlgJtEta_[nJtAlg][nJtMax];
Float_t AlgJtTrkMax_[nJtAlg][nJtMax];
Float_t AlgJtRawPt_[nJtAlg][nJtMax];

Float_t AlgJtAvePhi_[nJtAlg];
Float_t AlgJtDelPhi12_[nJtAlg];
Float_t AlgJtDelPhi13_[nJtAlg];
Float_t AlgJtDelPhi23_[nJtAlg];
Float_t AlgJtAsymm_[nJtAlg];

Float_t AlgRefPt_[nJtAlg][nJtMax];
Float_t AlgRefPhi_[nJtAlg][nJtMax];
Float_t AlgRefEta_[nJtAlg][nJtMax];

Float_t AlgRefAvePhi_[nJtAlg];
Float_t AlgRefDelPhi12_[nJtAlg];
Float_t AlgRefDelPhi13_[nJtAlg];
Float_t AlgRefDelPhi23_[nJtAlg];
Float_t AlgRefAsymm_[nJtAlg];

Float_t eventPt_[2];

//Gen Tree Variables

//Gen. proj. onto Jets, ordered by algorithm according to enum, PuPF == [0], PuCalo == [1], etc.

Float_t gAlgJtMult_[nSumAlg][2];

Float_t gAlgImbProjA_[nSumAlg][nPtBins];
Float_t gAlgImbProjAC_[nSumAlg][nPtBins];
Float_t gAlgImbProjANC_[nSumAlg][nPtBins];

//truth delRs
//ProjA delR

Float_t gAlgImbProjAR_[nSumAlg][nPtBins][nRBins];
Float_t gAlgImbProjAEta_[nSumAlg][nPtBins][nRBins];
Float_t gAlgImbProjAPhi_[nSumAlg][nPtBins][nRBins];

Float_t gAlgImbProjAR_Cut_[nSumAlg][nPtBins][nRBins];
Float_t gAlgImbProjAR_CutEta_[nSumAlg][nPtBins][nRBins];
Float_t gAlgImbProjAR_CutPhi_[nSumAlg][nPtBins][nRBins];
Float_t gAlgImbProjAEta_Cut_[nSumAlg][nPtBins][nRBins];
Float_t gAlgImbProjAPhi_Cut_[nSumAlg][nPtBins][nRBins];

Float_t gAlgImbProjAR_FOR_[nSumAlg][nPtBins][nRBins];

//delR Mult

Float_t gAlgMultAR_[nSumAlg][nPtBins][nRBins];
Float_t gAlgMultAEta_[nSumAlg][nPtBins][nRBins];
Float_t gAlgMultAPhi_[nSumAlg][nPtBins][nRBins];

Float_t gAlgMultAR_Cut_[nSumAlg][nPtBins][nRBins];
Float_t gAlgMultAR_CutEta_[nSumAlg][nPtBins][nRBins];
Float_t gAlgMultAR_CutPhi_[nSumAlg][nPtBins][nRBins];
Float_t gAlgMultAEta_Cut_[nSumAlg][nPtBins][nRBins];
Float_t gAlgMultAPhi_Cut_[nSumAlg][nPtBins][nRBins];

Float_t gAlgImbRawAR_[nSumAlg][nPtBins][nRBins];


//Cut var
const Float_t rBounds[nRBins] = {.20, .40, .60, .80, 1.00, 1.20, 1.40, 1.60, 1.80, 100000};

//const Float_t leadJtPtCut[nJtAlg] = {120., 120., 120., 117.448, 120., 123.66, 127.89, 120., 120., 120.};
//const Float_t subLeadJtPtCut[nJtAlg] = {50., 50., 50., 48.9367, 50., 51.5252, 53.2873, 50., 50., 50.};

const Float_t leadJtPtCut[nJtAlg] = {120., 120., 120., 120., 120., 120., 120., 120., 120., 120., 120., 120., 120., 120., 120., 120., 120., 120.};
const Float_t subLeadJtPtCut[nJtAlg] = {50., 50., 50., 50., 50., 50., 50., 50., 50., 50., 50., 50., 50., 50., 50., 50., 50., 50.};

const Float_t jtDelPhiCut = 0;
const Float_t jtEtaCut = 2.0; // Default Max at 2.4 to avoid transition junk, otherwise vary as needed                                                 
const std::string algType[nJtAlg] = {"Pu3Calo", "Pu4Calo", "Pu5Calo", "Vs2Calo", "Vs3Calo", "Vs4Calo", "Vs5Calo", "Vs2CaloFrag", "Vs3CaloFrag", "Vs4CaloFrag", "Vs5CaloFrag", "Vs3CaloRes", "T2", "T3", "T4", "T5", "PuPF", "VsPF"};


void SetAnaBranches(sampleType sType = kHIDATA, Bool_t justJt = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

  //Track Tree Branches

  std::cout << "Branches Set" << std::endl;

  //Tracks proj. onto Alg, ordered according to enum above, All, Cone, and NotCone

  trackTreeAna_p->Branch("rAlgJtMult", rAlgJtMult_, Form("rAlgJtMult[%d][2]/F", 2*nSumAlg));

  trackTreeAna_p->Branch("rAlgImbProjA", rAlgImbProjA_, Form("rAlgImbProjA[%d][%d]/F", 2*nSumAlg, nPtBins));
  trackTreeAna_p->Branch("rAlgImbProjAC", rAlgImbProjAC_, Form("rAlgImbProjAC[%d][%d]/F", 2*nSumAlg, nPtBins));
  trackTreeAna_p->Branch("rAlgImbProjANC", rAlgImbProjANC_, Form("rAlgImbProjANC[%d][%d]/F", 2*nSumAlg, nPtBins));
  
  trackTreeAna_p->Branch("rAlgImbProjAC0", rAlgImbProjAC0_, Form("rAlgImbProjAC0[%d][%d]/F", 2*nSumAlg, nPtBins));
  trackTreeAna_p->Branch("rAlgImbProjAC1", rAlgImbProjAC1_, Form("rAlgImbProjAC1[%d][%d]/F", 2*nSumAlg, nPtBins));
  trackTreeAna_p->Branch("rAlgImbProjAC2", rAlgImbProjAC2_, Form("rAlgImbProjAC2[%d][%d]/F", 2*nSumAlg, nPtBins));
  trackTreeAna_p->Branch("rAlgImbProjAC3", rAlgImbProjAC3_, Form("rAlgImbProjAC3[%d][%d]/F", 2*nSumAlg, nPtBins));

  //ProjA DelRs
  
  trackTreeAna_p->Branch("rAlgImbProjAR", rAlgImbProjAR_, Form("rAlgImbProjAR[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgImbProjAEta", rAlgImbProjAEta_, Form("rAlgImbProjAEta[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgImbProjAPhi", rAlgImbProjAPhi_, Form("rAlgImbProjAPhi[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));

  trackTreeAna_p->Branch("rAlgImbProjAR_Cut", rAlgImbProjAR_Cut_, Form("rAlgImbProjAR_Cut[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgImbProjAR_CutEta", rAlgImbProjAR_CutEta_, Form("rAlgImbProjAR_CutEta[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgImbProjAR_CutPhi", rAlgImbProjAR_CutPhi_, Form("rAlgImbProjAR_CutPhi[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgImbProjAEta_Cut", rAlgImbProjAEta_Cut_, Form("rAlgImbProjAEta_Cut[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgImbProjAPhi_Cut", rAlgImbProjAPhi_Cut_, Form("rAlgImbProjAPhi_Cut[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));

  trackTreeAna_p->Branch("rAlgImbProjAR_FOR", rAlgImbProjAR_FOR_, Form("rAlgImbProjAR_FOR[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));

  //Mult DelRs

  trackTreeAna_p->Branch("rAlgMultAR", rAlgMultAR_, Form("rAlgMultAR[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgMultAEta", rAlgMultAEta_, Form("rAlgMultAEta[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgMultAPhi", rAlgMultAPhi_, Form("rAlgMultAPhi[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));

  trackTreeAna_p->Branch("rAlgMultAR_Cut", rAlgMultAR_Cut_, Form("rAlgMultAR_Cut[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgMultAR_CutEta", rAlgMultAR_CutEta_, Form("rAlgMultAR_CutEta[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgMultAR_CutPhi", rAlgMultAR_CutPhi_, Form("rAlgMultAR_CutPhi[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgMultAEta_Cut", rAlgMultAEta_Cut_, Form("rAlgMultAEta_Cut[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));
  trackTreeAna_p->Branch("rAlgMultAPhi_Cut", rAlgMultAPhi_Cut_, Form("rAlgMultAPhi_Cut[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));

  trackTreeAna_p->Branch("rAlgImbRawAR", rAlgImbRawAR_, Form("rAlgImbRawAR[%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRBins));

  //Jet Tree Branches

  jetTreeAna_p->Branch("run", &run_, "run/I");
  jetTreeAna_p->Branch("evt", &evt_, "evt/I");
  jetTreeAna_p->Branch("lumi", &lumi_, "lumi/I");

  if(hi){
    jetTreeAna_p->Branch("hiBin", &hiBin_, "hiBin/I");
    jetTreeAna_p->Branch("hiEvtPlane", &hiEvtPlane_, "hiEvtPlane/F");
    jetTreeAna_p->Branch("psin", &psin_, "psin/F");
  }

  jetTreeAna_p->Branch("eventSet", eventSet_, Form("eventSet[%d]/O", nJtAlg));

  if(sType == kHIMC){
    jetTreeAna_p->Branch("centWeight", centWeight_, Form("centWeight[%d]/F", nJtAlg));
    jetTreeAna_p->Branch("centWeight_Merge", centWeight_Merge_, Form("centWeight_Merge[%d]/F", nJtAlg));
  }

  jetTreeAna_p->Branch("AlgJtPt", AlgJtPt_, Form("AlgJtPt[%d][%d]/F", nJtAlg, nJtMax));
  jetTreeAna_p->Branch("AlgJtPhi", AlgJtPhi_, Form("AlgJtPhi[%d][%d]/F", nJtAlg, nJtMax));
  jetTreeAna_p->Branch("AlgJtEta", AlgJtEta_, Form("AlgJtEta[%d][%d]/F", nJtAlg, nJtMax));
  jetTreeAna_p->Branch("AlgJtTrkMax", AlgJtTrkMax_, Form("AlgJtTrkMax[%d][%d]/F", nJtAlg, nJtMax));
  jetTreeAna_p->Branch("AlgJtRawPt", AlgJtRawPt_, Form("AlgJtRawPt[%d][%d]/F", nJtAlg, nJtMax));

  jetTreeAna_p->Branch("AlgJtAvePhi", AlgJtAvePhi_, Form("AlgJtAvePhi[%d]/F", nJtAlg));
  jetTreeAna_p->Branch("AlgJtDelPhi12", AlgJtDelPhi12_, Form("AlgJtDelPhi12[%d]/F", nJtAlg));
  jetTreeAna_p->Branch("AlgJtDelPhi13", AlgJtDelPhi13_, Form("AlgJtDelPhi13[%d]/F", nJtAlg));
  jetTreeAna_p->Branch("AlgJtDelPhi23", AlgJtDelPhi23_, Form("AlgJtDelPhi23[%d]/F", nJtAlg));
  jetTreeAna_p->Branch("AlgJtAsymm", AlgJtAsymm_, Form("AlgJtAsymm[%d]/F", nJtAlg));

  if(montecarlo){
    jetTreeAna_p->Branch("isQuarkJet", isQuarkJet_, Form("isQuarkJet[%d]/O", nJtAlg));
    jetTreeAna_p->Branch("isGluonJet", isGluonJet_, Form("isGluonJet[%d]/O", nJtAlg));
  }

  jetTreeAna_p->Branch("pthatWeight", &pthatWeight_, "pthatWeight/F");
  jetTreeAna_p->Branch("leadEtaWeight", &leadEtaWeight_, "leadEtaWeight/F");
  jetTreeAna_p->Branch("subleadEtaWeight", &subleadEtaWeight_, "subleadEtaWeight/F");

  if(justJt)
    jetTreeAna_p->Branch("eventPt", eventPt_, "eventPt[2]/F");

  if(montecarlo){
    jetTreeAna_p->Branch("pthat", &pthat_, "pthat/F");

    //refpt for jets immediately above
    jetTreeAna_p->Branch("AlgRefPt", AlgRefPt_, Form("AlgRefPt[%d][%d]/F", nJtAlg, nJtMax));
    jetTreeAna_p->Branch("AlgRefPhi", AlgRefPhi_, Form("AlgRefPhi[%d][%d]/F", nJtAlg, nJtMax));
    jetTreeAna_p->Branch("AlgRefEta", AlgRefEta_, Form("AlgRefEta[%d][%d]/F", nJtAlg, nJtMax));

    jetTreeAna_p->Branch("AlgRefAvePhi", AlgRefAvePhi_, Form("AlgRefAvePhi[%d]/F", nJtAlg));
    jetTreeAna_p->Branch("AlgRefDelPhi12", AlgRefDelPhi12_, Form("AlgRefDelPhi12[%d]/F", nJtAlg));
    jetTreeAna_p->Branch("AlgRefDelPhi13", AlgRefDelPhi13_, Form("AlgRefDelPhi13[%d]/F", nJtAlg));
    jetTreeAna_p->Branch("AlgRefDelPhi23", AlgRefDelPhi23_, Form("AlgRefDelPhi23[%d]/F", nJtAlg));
    jetTreeAna_p->Branch("AlgRefAsymm", AlgRefAsymm_, Form("AlgRefAsymm[%d]/F", nJtAlg));


    //Gen. proj. onto jetAlg, array ordered according to enum

    if(!justJt){
      genTreeAna_p->Branch("gAlgJtMult", gAlgJtMult_, Form("gAlgJtMult[%d][2]/F", nSumAlg));
      genTreeAna_p->Branch("gAlgImbProjA", gAlgImbProjA_, Form("gAlgImbProjA[%d][%d]/F", nSumAlg, nPtBins));
      genTreeAna_p->Branch("gAlgImbProjAC", gAlgImbProjAC_, Form("gAlgImbProjAC[%d][%d]/F", nSumAlg, nPtBins));
      genTreeAna_p->Branch("gAlgImbProjANC", gAlgImbProjANC_, Form("gAlgImbProjANC[%d][%d]/F", nSumAlg, nPtBins));
      
      //Truth Del Rs, Proj
      //Proj A's
      
      genTreeAna_p->Branch("gAlgImbProjAR", gAlgImbProjAR_, Form("gAlgImbProjAR[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgImbProjAEta", gAlgImbProjAEta_, Form("gAlgImbProjAEta[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgImbProjAPhi", gAlgImbProjAPhi_, Form("gAlgImbProjAPhi[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));

      genTreeAna_p->Branch("gAlgImbProjAR_Cut", gAlgImbProjAR_Cut_, Form("gAlgImbProjAR_Cut[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgImbProjAR_CutEta", gAlgImbProjAR_CutEta_, Form("gAlgImbProjAR_CutEta[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgImbProjAR_CutPhi", gAlgImbProjAR_CutPhi_, Form("gAlgImbProjAR_CutPhi[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgImbProjAEta_Cut", gAlgImbProjAEta_Cut_, Form("gAlgImbProjAEta_Cut[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgImbProjAPhi_Cut", gAlgImbProjAPhi_Cut_, Form("gAlgImbProjAPhi_Cut[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));

      genTreeAna_p->Branch("gAlgImbProjAR_FOR", gAlgImbProjAR_FOR_, Form("gAlgImbProjAR_FOR[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));

      //Mult DelRs

      genTreeAna_p->Branch("gAlgMultAR", gAlgMultAR_, Form("gAlgMultAR[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgMultAEta", gAlgMultAEta_, Form("gAlgMultAEta[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgMultAPhi", gAlgMultAPhi_, Form("gAlgMultAPhi[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));

      genTreeAna_p->Branch("gAlgMultAR_Cut", gAlgMultAR_Cut_, Form("gAlgMultAR_Cut[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgMultAR_CutEta", gAlgMultAR_CutEta_, Form("gAlgMultAR_CutEta[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgMultAR_CutPhi", gAlgMultAR_CutPhi_, Form("gAlgMultAR_CutPhi[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgMultAEta_Cut", gAlgMultAEta_Cut_, Form("gAlgMultAEta_Cut[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
      genTreeAna_p->Branch("gAlgMultAPhi_Cut", gAlgMultAPhi_Cut_, Form("gAlgMultAPhi_Cut[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));

      genTreeAna_p->Branch("gAlgImbRawAR", gAlgImbRawAR_, Form("gAlgImbRawAR[%d][%d][%d]/F", nSumAlg, nPtBins, nRBins));
    }
  }
}


void GetAnaBranches(sampleType sType = kHIDATA, Bool_t justJt = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);
  
  std::cout << "Get Branches" << std::endl;
  
  //Track Tree Branches
  
  trackTreeAna_p->SetBranchAddress("rAlgJtMult", rAlgJtMult_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjA", rAlgImbProjA_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAC", rAlgImbProjAC_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjANC", rAlgImbProjANC_);

  trackTreeAna_p->SetBranchAddress("rAlgImbProjAR", rAlgImbProjAR_);    
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAEta", rAlgImbProjAEta_);    
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAPhi", rAlgImbProjAPhi_);    

  trackTreeAna_p->SetBranchAddress("rAlgImbProjAR_Cut", rAlgImbProjAR_Cut_);    
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAR_CutEta", rAlgImbProjAR_CutEta_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAR_CutPhi", rAlgImbProjAR_CutPhi_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAEta_Cut", rAlgImbProjAEta_Cut_);    
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAPhi_Cut", rAlgImbProjAPhi_Cut_);    

  trackTreeAna_p->SetBranchAddress("rAlgImbProjAR_FOR", rAlgImbProjAR_FOR_);    

  //Mult DelRs

  trackTreeAna_p->SetBranchAddress("rAlgMultAR", rAlgMultAR_);    
  trackTreeAna_p->SetBranchAddress("rAlgMultAEta", rAlgMultAEta_);    
  trackTreeAna_p->SetBranchAddress("rAlgMultAPhi", rAlgMultAPhi_);    

  trackTreeAna_p->SetBranchAddress("rAlgMultAR_Cut", rAlgMultAR_Cut_);    
  trackTreeAna_p->SetBranchAddress("rAlgMultAR_CutEta", rAlgMultAR_CutEta_);    
  trackTreeAna_p->SetBranchAddress("rAlgMultAR_CutPhi", rAlgMultAR_CutPhi_);    
  trackTreeAna_p->SetBranchAddress("rAlgMultAEta_Cut", rAlgMultAEta_Cut_);    
  trackTreeAna_p->SetBranchAddress("rAlgMultAPhi_Cut", rAlgMultAPhi_Cut_);    

  trackTreeAna_p->SetBranchAddress("rAlgImbRawAR", rAlgImbRawAR_);    
  
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAC0", rAlgImbProjAC0_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAC1", rAlgImbProjAC1_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAC2", rAlgImbProjAC2_);
  trackTreeAna_p->SetBranchAddress("rAlgImbProjAC3", rAlgImbProjAC3_);
  
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
    jetTreeAna_p->SetBranchAddress("centWeight", centWeight_);
    jetTreeAna_p->SetBranchAddress("centWeight_Merge", centWeight_Merge_);
  }
  
  if(montecarlo){
    jetTreeAna_p->SetBranchAddress("isQuarkJet", isQuarkJet_);
    jetTreeAna_p->SetBranchAddress("isGluonJet", isGluonJet_);
  }     

  jetTreeAna_p->SetBranchAddress("pthatWeight", &pthatWeight_);
  jetTreeAna_p->SetBranchAddress("leadEtaWeight", &leadEtaWeight_);
  jetTreeAna_p->SetBranchAddress("subleadEtaWeight", &subleadEtaWeight_);
  
  jetTreeAna_p->SetBranchAddress("AlgJtPt", AlgJtPt_);
  jetTreeAna_p->SetBranchAddress("AlgJtPhi", AlgJtPhi_);
  jetTreeAna_p->SetBranchAddress("AlgJtEta", AlgJtEta_);
  jetTreeAna_p->SetBranchAddress("AlgJtTrkMax", AlgJtTrkMax_);
  jetTreeAna_p->SetBranchAddress("AlgJtRawPt", AlgJtRawPt_);
  
  jetTreeAna_p->SetBranchAddress("AlgJtAvePhi", AlgJtAvePhi_);
  jetTreeAna_p->SetBranchAddress("AlgJtDelPhi12", AlgJtDelPhi12_);
  jetTreeAna_p->SetBranchAddress("AlgJtDelPhi13", AlgJtDelPhi13_);
  jetTreeAna_p->SetBranchAddress("AlgJtDelPhi23", AlgJtDelPhi23_);
  jetTreeAna_p->SetBranchAddress("AlgJtAsymm", AlgJtAsymm_);
  
  if(justJt)
    jetTreeAna_p->SetBranchAddress("eventPt", eventPt_);
  
  if(montecarlo){
    jetTreeAna_p->SetBranchAddress("AlgRefPt", AlgRefPt_);
    jetTreeAna_p->SetBranchAddress("AlgRefPhi", AlgRefPhi_);
    jetTreeAna_p->SetBranchAddress("AlgRefEta", AlgRefEta_);
    
    jetTreeAna_p->SetBranchAddress("AlgRefAvePhi", AlgRefAvePhi_);
    jetTreeAna_p->SetBranchAddress("AlgRefDelPhi12", AlgRefDelPhi12_);
    jetTreeAna_p->SetBranchAddress("AlgRefDelPhi13", AlgRefDelPhi13_);
    jetTreeAna_p->SetBranchAddress("AlgRefDelPhi23", AlgRefDelPhi23_);
    jetTreeAna_p->SetBranchAddress("AlgRefAsymm", AlgRefAsymm_);
    
    //Gen Tree Variables
    
    if(!justJt){
      genTreeAna_p->SetBranchAddress("gAlgJtMult", gAlgJtMult_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjA", gAlgImbProjA_);     
      genTreeAna_p->SetBranchAddress("gAlgImbProjAC", gAlgImbProjAC_);     
      genTreeAna_p->SetBranchAddress("gAlgImbProjANC", gAlgImbProjANC_);     
      
      genTreeAna_p->SetBranchAddress("gAlgImbProjAR", gAlgImbProjAR_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjAEta", gAlgImbProjAEta_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjAPhi", gAlgImbProjAPhi_);

      genTreeAna_p->SetBranchAddress("gAlgImbProjAR_Cut", gAlgImbProjAR_Cut_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjAR_CutEta", gAlgImbProjAR_CutEta_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjAR_CutPhi", gAlgImbProjAR_CutPhi_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjAEta_Cut", gAlgImbProjAEta_Cut_);
      genTreeAna_p->SetBranchAddress("gAlgImbProjAPhi_Cut", gAlgImbProjAPhi_Cut_);

      genTreeAna_p->SetBranchAddress("gAlgImbProjAR_FOR", gAlgImbProjAR_FOR_);

      //Mult DelRs

      genTreeAna_p->SetBranchAddress("gAlgMultAR", gAlgMultAR_);
      genTreeAna_p->SetBranchAddress("gAlgMultAEta", gAlgMultAEta_);
      genTreeAna_p->SetBranchAddress("gAlgMultAPhi", gAlgMultAPhi_);

      genTreeAna_p->SetBranchAddress("gAlgMultAR_Cut", gAlgMultAR_Cut_);
      genTreeAna_p->SetBranchAddress("gAlgMultAR_CutEta", gAlgMultAR_CutEta_);
      genTreeAna_p->SetBranchAddress("gAlgMultAR_CutPhi", gAlgMultAR_CutPhi_);
      genTreeAna_p->SetBranchAddress("gAlgMultAEta_Cut", gAlgMultAEta_Cut_);
      genTreeAna_p->SetBranchAddress("gAlgMultAPhi_Cut", gAlgMultAPhi_Cut_);

      genTreeAna_p->SetBranchAddress("gAlgImbRawAR", gAlgImbRawAR_);
    }
  }
}


void InitDiJetAnaSkim(sampleType sType = kHIDATA, Bool_t justJt = false)
{
  std::cout << "Init DiJet AnaSkim" << std::endl;

  trackTreeAna_p = new TTree("trackTreeAna", "trackTreeAna");
  jetTreeAna_p = new TTree("jetTreeAna", "jetTreeAna");
  if(isMonteCarlo(sType)) genTreeAna_p = new TTree("genTreeAna", "genTreeAna");

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
  if(isMonteCarlo(sType) && !justJt) genTreeAna_p = (TTree*)anaFile_p->Get("genTreeAna");

  GetAnaBranches(sType, justJt);
}


Float_t getAbsDphi(Float_t phi1, Float_t phi2, Int_t tag, Int_t tag2 = -1, Int_t tag3 = -1)
{
  return TMath::Abs(getDPHI(phi1, phi2, tag, tag2, tag3));
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



void InitJetVar(sampleType sType = kHIDATA)
{
  Bool_t montecarlo = isMonteCarlo(sType);

  for(Int_t initIter = 0; initIter < nJtAlg; initIter++){
    eventSet_[initIter] = false;

    if(sType == kHIMC){
      centWeight_[initIter] = 1;    
      centWeight_Merge_[initIter] = 1;    
    }

    for(Int_t initIter2 = 0; initIter2 < nJtMax; initIter2++){
      AlgJtPt_[initIter][initIter2] = -999;
      AlgJtPhi_[initIter][initIter2] = -999;
      AlgJtEta_[initIter][initIter2] = -999;
      AlgJtTrkMax_[initIter][initIter2] = -999;
      AlgJtRawPt_[initIter][initIter2] = -999;
    }

    AlgJtAvePhi_[initIter] = -999;
    AlgJtDelPhi12_[initIter] = -999;
    AlgJtDelPhi13_[initIter] = -999;
    AlgJtDelPhi23_[initIter] = -999;
    AlgJtAsymm_[initIter] = -999;

    if(montecarlo){
      isQuarkJet_[initIter] = false;
      isGluonJet_[initIter] = false;

      pthatWeight_ = -999;
      leadEtaWeight_ = -999;
      subleadEtaWeight_ = -999;

      for(Int_t initIter2 = 0; initIter2 < nJtMax; initIter2++){
	AlgRefPt_[initIter][initIter2] = -999;
	AlgRefPhi_[initIter][initIter2] = -999;
	AlgRefEta_[initIter][initIter2] = -999;
      }

      AlgRefAvePhi_[initIter] = -999;
      AlgRefDelPhi12_[initIter] = -999;
      AlgRefDelPhi13_[initIter] = -999;
      AlgRefDelPhi23_[initIter] = -999;
      AlgRefAsymm_[initIter] = -999;
    }
  }
}

void InitProjPerp(sampleType sType = kHIDATA)
{
  //Tracks proj. onto Alg, ordered according to enum above, corr in the back 5, All, Cone, and NotCone

  for(Int_t initIter = 0; initIter < 2*nSumAlg; initIter++){
    for(Int_t initIter2 = 0; initIter2 < nPtBins; initIter2++){
      if(initIter2 < 2) rAlgJtMult_[initIter][initIter2] = 0;

      rAlgImbProjA_[initIter][initIter2] = 0;
      rAlgImbProjAC_[initIter][initIter2] = 0;
      rAlgImbProjANC_[initIter][initIter2] = 0;

      rAlgImbProjAC0_[initIter][initIter2] = 0;
      rAlgImbProjAC1_[initIter][initIter2] = 0;
      rAlgImbProjAC2_[initIter][initIter2] = 0;
      rAlgImbProjAC3_[initIter][initIter2] = 0;

      for(Int_t initIter3 = 0; initIter3 < nRBins; initIter3++){
	rAlgImbProjAR_[initIter][initIter2][initIter3] = 0;
	rAlgImbProjAEta_[initIter][initIter2][initIter3] = 0;
	rAlgImbProjAPhi_[initIter][initIter2][initIter3] = 0;

	rAlgImbProjAR_Cut_[initIter][initIter2][initIter3] = 0;
	rAlgImbProjAR_CutEta_[initIter][initIter2][initIter3] = 0;
	rAlgImbProjAR_CutPhi_[initIter][initIter2][initIter3] = 0;
	rAlgImbProjAEta_Cut_[initIter][initIter2][initIter3] = 0;
	rAlgImbProjAPhi_Cut_[initIter][initIter2][initIter3] = 0;

	rAlgImbProjAR_FOR_[initIter][initIter2][initIter3] = 0;

	//Init Mults

	rAlgMultAR_[initIter][initIter2][initIter3] = 0;
	rAlgMultAEta_[initIter][initIter2][initIter3] = 0;
	rAlgMultAPhi_[initIter][initIter2][initIter3] = 0;

	rAlgMultAR_Cut_[initIter][initIter2][initIter3] = 0;
	rAlgMultAR_CutEta_[initIter][initIter2][initIter3] = 0;
	rAlgMultAR_CutPhi_[initIter][initIter2][initIter3] = 0;
	rAlgMultAEta_Cut_[initIter][initIter2][initIter3] = 0;
	rAlgMultAPhi_Cut_[initIter][initIter2][initIter3] = 0;

	rAlgImbRawAR_[initIter][initIter2][initIter3] = 0;
      }
    }
  }

  if(isMonteCarlo(sType)){
    //Gen. proj. onto Truth
    for(Int_t initIter = 0; initIter < nSumAlg; initIter++){
      for(Int_t initIter2 = 0; initIter2 < nPtBins; initIter2++){
	if(initIter2 < 2) gAlgJtMult_[initIter][initIter2] = 0;

	gAlgImbProjA_[initIter][initIter2] = 0;
	gAlgImbProjAC_[initIter][initIter2] = 0;
	gAlgImbProjANC_[initIter][initIter2] = 0;

	for(Int_t initIter3 = 0; initIter3 < nRBins; initIter3++){
	  gAlgImbProjAR_[initIter][initIter2][initIter3] = 0;
	  gAlgImbProjAEta_[initIter][initIter2][initIter3] = 0;
	  gAlgImbProjAPhi_[initIter][initIter2][initIter3] = 0;

	  gAlgImbProjAR_Cut_[initIter][initIter2][initIter3] = 0;
	  gAlgImbProjAR_CutEta_[initIter][initIter2][initIter3] = 0;
	  gAlgImbProjAR_CutPhi_[initIter][initIter2][initIter3] = 0;
	  gAlgImbProjAEta_Cut_[initIter][initIter2][initIter3] = 0;
	  gAlgImbProjAPhi_Cut_[initIter][initIter2][initIter3] = 0;

	  gAlgImbProjAR_FOR_[initIter][initIter2][initIter3] = 0;

	  //Init Mults

	  gAlgMultAR_[initIter][initIter2][initIter3] = 0;
	  gAlgMultAEta_[initIter][initIter2][initIter3] = 0;
	  gAlgMultAPhi_[initIter][initIter2][initIter3] = 0;

	  gAlgMultAR_Cut_[initIter][initIter2][initIter3] = 0;
	  gAlgMultAR_CutEta_[initIter][initIter2][initIter3] = 0;
	  gAlgMultAR_CutPhi_[initIter][initIter2][initIter3] = 0;
	  gAlgMultAEta_Cut_[initIter][initIter2][initIter3] = 0;
	  gAlgMultAPhi_Cut_[initIter][initIter2][initIter3] = 0;

	  gAlgImbRawAR_[initIter][initIter2][initIter3] = 0;
	}
      }
    }    
  }
}


void getJtVar(Int_t nJt, Float_t jtPt[], Float_t jtPhi[], Float_t jtEta[], Float_t trkMax[], Float_t rawPt[], Float_t refPt[], Float_t refPhi[], Float_t refEta[], Int_t algNum, Bool_t montecarlo = false, Bool_t truth = false)
{
  if(nJt == 0 || nJt == 1 || jtPt[0] < leadJtPtCut[algNum] || jtPt[1] < subLeadJtPtCut[algNum]) return;

  eventSet_[algNum] = true;

  Int_t iterMax = nJt;
  if(nJt > nJtMax) iterMax = nJtMax;

  for(Int_t leadIter = 0; leadIter < iterMax; leadIter++){
    AlgJtPt_[algNum][leadIter] = jtPt[leadIter];
    AlgJtPhi_[algNum][leadIter] = jtPhi[leadIter];
    AlgJtEta_[algNum][leadIter] = jtEta[leadIter];

    if(!truth){
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
  AlgJtDelPhi12_[algNum] = getAbsDphi(AlgJtPhi_[algNum][0], AlgJtPhi_[algNum][1], 5);
  
  if(AlgJtPt_[algNum][2] > 30 && TMath::Abs(AlgJtEta_[algNum][2]) < 2.0){
    AlgJtDelPhi13_[algNum] = getAbsDphi(AlgJtPhi_[algNum][0], AlgJtPhi_[algNum][2], 6);
    AlgJtDelPhi23_[algNum] = getAbsDphi(AlgJtPhi_[algNum][1], AlgJtPhi_[algNum][2], 7);
  }
  
  AlgJtAsymm_[algNum] = (AlgJtPt_[algNum][0] - AlgJtPt_[algNum][1])/(AlgJtPt_[algNum][0] + AlgJtPt_[algNum][1]);
    
  if(montecarlo && !truth){
    if(AlgRefPhi_[algNum][0] > -9 && AlgRefPhi_[algNum][1] > -9){
      AlgRefAvePhi_[algNum] = getAvePhi(AlgRefPhi_[algNum][0], AlgRefPhi_[algNum][1]);

      AlgRefDelPhi12_[algNum] = getAbsDphi(AlgRefPhi_[algNum][0], AlgRefPhi_[algNum][1], 8);

      if(AlgRefPhi_[algNum][2] > -9 && nJt > 2){
	AlgRefDelPhi13_[algNum] = getAbsDphi(AlgRefPhi_[algNum][0], AlgRefPhi_[algNum][2], 9);
	AlgRefDelPhi23_[algNum] = getAbsDphi(AlgRefPhi_[algNum][1], AlgRefPhi_[algNum][2], 10);
      }

      AlgRefAsymm_[algNum] = (AlgRefPt_[algNum][0] - AlgRefPt_[algNum][1])/(AlgRefPt_[algNum][0] + AlgRefPt_[algNum][1]);
    }
    else{
      AlgRefAvePhi_[algNum] = -999;
      AlgRefDelPhi12_[algNum] = -999;
      AlgRefDelPhi13_[algNum] = -999;
      AlgRefDelPhi23_[algNum] = -999;
      AlgRefAsymm_[algNum] = -999;
    }
  }

  return;
}

void GetTrkProjPerp(Int_t jtAlg, Int_t jtAlgCorr, Float_t rawPt, Float_t corrPt, Float_t phi, Float_t eta)
{
  if(rawPt < 0.5) return;

  Int_t multVal;
  if(jtAlg == jtAlgCorr) multVal = 1;
  else multVal = corrPt/rawPt;

  Int_t signVal = 1;
  if(cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2)) < 0) signVal = -1;

  if(getAbsDphi(AlgJtAvePhi_[jtAlg], phi, 2) < TMath::Pi()/2) rAlgJtMult_[jtAlgCorr][0] += multVal;
  else rAlgJtMult_[jtAlgCorr][1] += multVal;

  Int_t ptIter = getPtRange(rawPt);
  Float_t tempLeadDelR = getDR(eta, phi, AlgJtEta_[jtAlg][0], AlgJtPhi_[jtAlg][0], 2);
  Float_t tempSubLeadDelR = getDR(eta, phi, AlgJtEta_[jtAlg][1], AlgJtPhi_[jtAlg][1], 2);
  Float_t tempLeadDelEta = TMath::Abs(eta - AlgJtEta_[jtAlg][0]);
  Float_t tempSubLeadDelEta = TMath::Abs(eta - AlgJtEta_[jtAlg][1]);
  Float_t tempLeadDelPhi = TMath::Abs(getDPHI(phi, AlgJtPhi_[jtAlg][0], 2));
  Float_t tempSubLeadDelPhi = TMath::Abs(getDPHI(phi, AlgJtPhi_[jtAlg][1], 2));

  rAlgImbProjA_[jtAlgCorr][5] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
  rAlgImbProjA_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));

  if(tempLeadDelR > 0 && tempSubLeadDelR > 0){
    if(tempLeadDelR < 0.8 || tempSubLeadDelR < 0.8){
      rAlgImbProjAC_[jtAlgCorr][5] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
      rAlgImbProjAC_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
    }
    else{
      rAlgImbProjANC_[jtAlgCorr][5] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
      rAlgImbProjANC_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
    }

    if(tempLeadDelR < 0.5 || tempSubLeadDelR < 0.5){
      rAlgImbProjAC0_[jtAlgCorr][5] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
      rAlgImbProjAC0_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
    }
    else if(tempLeadDelR < 1.0 || tempSubLeadDelR < 1.0){
      rAlgImbProjAC1_[jtAlgCorr][5] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
      rAlgImbProjAC1_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
    }
    else if(tempLeadDelR < 1.5 || tempSubLeadDelR < 1.5){
      rAlgImbProjAC2_[jtAlgCorr][5] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
      rAlgImbProjAC2_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
    }
    else{
      rAlgImbProjAC3_[jtAlgCorr][5] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
      rAlgImbProjAC3_[jtAlgCorr][ptIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	rAlgImbProjAR_[jtAlgCorr][5][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));	
	rAlgImbProjAR_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));

	rAlgMultAR_[jtAlgCorr][5][rIter] -= signVal*multVal;	
	rAlgMultAR_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	rAlgImbRawAR_[jtAlgCorr][5][rIter] += TMath::Abs(corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2)));	
	rAlgImbRawAR_[jtAlgCorr][ptIter][rIter] += TMath::Abs(corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2)));
	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta > 1.0 && tempSubLeadDelEta > 1.0) break;
      if(tempLeadDelPhi > 1.0 && tempSubLeadDelPhi > 1.0) break;

      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	rAlgImbProjAR_Cut_[jtAlgCorr][5][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));	
	rAlgImbProjAR_Cut_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));

	rAlgMultAR_Cut_[jtAlgCorr][5][rIter] -= signVal*multVal;	
	rAlgMultAR_Cut_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta > 1.0 && tempSubLeadDelEta > 1.0) break;

      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	rAlgImbProjAR_CutEta_[jtAlgCorr][5][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));	
	rAlgImbProjAR_CutEta_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));

	rAlgMultAR_CutEta_[jtAlgCorr][5][rIter] -= signVal*multVal;	
	rAlgMultAR_CutEta_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelPhi > 1.0 && tempSubLeadDelPhi > 1.0) break;

      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	rAlgImbProjAR_CutPhi_[jtAlgCorr][5][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));	
	rAlgImbProjAR_CutPhi_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));

	rAlgMultAR_CutPhi_[jtAlgCorr][5][rIter] -= signVal*multVal;	
	rAlgMultAR_CutPhi_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelPhi < rBounds[rIter] || tempSubLeadDelPhi < rBounds[rIter]){
	rAlgImbProjAPhi_[jtAlgCorr][5][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));	
	rAlgImbProjAPhi_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));

	rAlgMultAPhi_[jtAlgCorr][5][rIter] -= signVal*multVal;
	rAlgMultAPhi_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta > 1.0 && tempSubLeadDelEta > 1.0) break;
      if(tempLeadDelPhi < rBounds[rIter] || tempSubLeadDelPhi < rBounds[rIter]){
	rAlgImbProjAPhi_Cut_[jtAlgCorr][5][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));	
	rAlgImbProjAPhi_Cut_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));

	rAlgMultAPhi_Cut_[jtAlgCorr][5][rIter] -= signVal*multVal;	
	rAlgMultAPhi_Cut_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta < rBounds[rIter] || tempSubLeadDelEta < rBounds[rIter]){
	rAlgImbProjAEta_[jtAlgCorr][5][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));	
	rAlgImbProjAEta_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));

	rAlgMultAEta_[jtAlgCorr][5][rIter] -= signVal*multVal;	
	rAlgMultAEta_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelPhi > 1.0 && tempSubLeadDelPhi > 1.0) break;
      if(tempLeadDelEta < rBounds[rIter] || tempSubLeadDelEta < rBounds[rIter]){
	rAlgImbProjAEta_Cut_[jtAlgCorr][5][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));	
	rAlgImbProjAEta_Cut_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));

	rAlgMultAEta_Cut_[jtAlgCorr][5][rIter] -= signVal*multVal;	
	rAlgMultAEta_Cut_[jtAlgCorr][ptIter][rIter] -= signVal*multVal;

	break;
      }
    }

    if((AlgJtEta_[jtAlg][0] < 0.6 && AlgJtEta_[jtAlg][1] < 0.6)){
      if(eta > TMath::Min(AlgJtEta_[jtAlg][0], AlgJtEta_[jtAlg][1]) - 0.3){
	for(Int_t rIter = 0; rIter < 10; rIter++){
	  if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	    rAlgImbProjAR_FOR_[jtAlgCorr][5][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
	    rAlgImbProjAR_FOR_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
	    break;
	  }
	}
      }
    }
    else if((AlgJtEta_[jtAlg][0] > -0.6 && AlgJtEta_[jtAlg][1] > -0.6)){
      if(eta < TMath::Max(AlgJtEta_[jtAlg][0], AlgJtEta_[jtAlg][1]) + 0.3){
	for(Int_t rIter = 0; rIter < 10; rIter++){
	  if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	    rAlgImbProjAR_FOR_[jtAlgCorr][5][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
	    rAlgImbProjAR_FOR_[jtAlgCorr][ptIter][rIter] -= corrPt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 2));
	    break;
	  }
	}
      }
    }
  
  }


  return;
}


void GetMixProjPerp(Float_t evPt[2], Int_t nConst, Float_t constPt[], Float_t constPhi[], Float_t constEta[], Float_t constCorr[])
{
  for(Int_t constIter = 0; constIter < nConst; constIter++){
    evPt[0] -= constPt[constIter]*constCorr[constIter]*cos(constPhi[constIter]);
    evPt[1] -= constPt[constIter]*constCorr[constIter]*sin(constPhi[constIter]);
    GetTrkProjPerp(1, 1, constPt[constIter], constPt[constIter]*constCorr[constIter], constPhi[constIter], constEta[constIter]);
    GetTrkProjPerp(1, 4, constPt[constIter], constPt[constIter]*constCorr[constIter], constPhi[constIter], constEta[constIter]);
  }

  return;
}


void GetGenProjPerp(Int_t jtAlg, Float_t pt, Float_t phi, Float_t eta)
{
  if(pt < 0.5) return;
  if(TMath::IsNaN(pt)) return;

  Int_t signVal = 1;
  if(cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3)) < 0) signVal = -1;

  if(getAbsDphi(AlgJtAvePhi_[jtAlg], phi, 3) < TMath::Pi()/2) gAlgJtMult_[jtAlg][0] += 1;
  else gAlgJtMult_[jtAlg][1] += 1;

  Int_t ptIter = getPtRange(pt);
  Float_t tempLeadDelR = getDR(eta, phi, AlgJtEta_[jtAlg][0], AlgJtPhi_[jtAlg][0], 3);
  Float_t tempSubLeadDelR = getDR(eta, phi, AlgJtEta_[jtAlg][1], AlgJtPhi_[jtAlg][1], 3);
  Float_t tempLeadDelEta = TMath::Abs(eta - AlgJtEta_[jtAlg][0]);
  Float_t tempSubLeadDelEta = TMath::Abs(eta - AlgJtEta_[jtAlg][1]);
  Float_t tempLeadDelPhi = TMath::Abs(getDPHI(phi, AlgJtPhi_[jtAlg][0], 3));
  Float_t tempSubLeadDelPhi = TMath::Abs(getDPHI(phi, AlgJtPhi_[jtAlg][1], 3));

  gAlgImbProjA_[jtAlg][5] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
  gAlgImbProjA_[jtAlg][ptIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));

  if(tempLeadDelR > 0 && tempSubLeadDelR > 0){
    if(tempLeadDelR < 0.8 || tempSubLeadDelR < 0.8){
      gAlgImbProjAC_[jtAlg][5] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
      gAlgImbProjAC_[jtAlg][ptIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
    }
    else{
      gAlgImbProjANC_[jtAlg][5] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
      gAlgImbProjANC_[jtAlg][ptIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	gAlgImbProjAR_[jtAlg][5][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
	gAlgImbProjAR_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));

	gAlgMultAR_[jtAlg][5][rIter] -= signVal;
	gAlgMultAR_[jtAlg][ptIter][rIter] -= signVal;

	gAlgImbRawAR_[jtAlg][5][rIter] += TMath::Abs(pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3)));
	gAlgImbRawAR_[jtAlg][ptIter][rIter] += TMath::Abs(pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3)));
	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta > 1.0 && tempSubLeadDelEta > 1.0) break;
      if(tempLeadDelPhi > 1.0 && tempSubLeadDelPhi > 1.0) break;

      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	gAlgImbProjAR_Cut_[jtAlg][5][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
	gAlgImbProjAR_Cut_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));

	gAlgMultAR_Cut_[jtAlg][5][rIter] -= signVal;
	gAlgMultAR_Cut_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta > 1.0 && tempSubLeadDelEta > 1.0) break;

      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	gAlgImbProjAR_CutEta_[jtAlg][5][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
	gAlgImbProjAR_CutEta_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));

	gAlgMultAR_CutEta_[jtAlg][5][rIter] -= signVal;
	gAlgMultAR_CutEta_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelPhi > 1.0 && tempSubLeadDelPhi > 1.0) break;

      if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	gAlgImbProjAR_CutPhi_[jtAlg][5][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
	gAlgImbProjAR_CutPhi_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));

	gAlgMultAR_CutPhi_[jtAlg][5][rIter] -= signVal;
	gAlgMultAR_CutPhi_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }


    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta < rBounds[rIter] || tempSubLeadDelEta < rBounds[rIter]){
	gAlgImbProjAEta_[jtAlg][5][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
	gAlgImbProjAEta_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));

	gAlgMultAEta_[jtAlg][5][rIter] -= signVal;
	gAlgMultAEta_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }


    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelPhi > 1.0 && tempSubLeadDelPhi > 1.0) break;

      if(tempLeadDelEta < rBounds[rIter] || tempSubLeadDelEta < rBounds[rIter]){
	gAlgImbProjAEta_Cut_[jtAlg][5][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
	gAlgImbProjAEta_Cut_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));

	gAlgMultAEta_Cut_[jtAlg][5][rIter] -= signVal;
	gAlgMultAEta_Cut_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelPhi < rBounds[rIter] || tempSubLeadDelPhi < rBounds[rIter]){
	gAlgImbProjAPhi_[jtAlg][5][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
	gAlgImbProjAPhi_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));

	gAlgMultAPhi_[jtAlg][5][rIter] -= signVal;
	gAlgMultAPhi_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }


    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(tempLeadDelEta > 1.0 && tempSubLeadDelEta > 1.0) break;

      if(tempLeadDelPhi < rBounds[rIter] || tempSubLeadDelPhi < rBounds[rIter]){
	gAlgImbProjAPhi_Cut_[jtAlg][5][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
	gAlgImbProjAPhi_Cut_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));

	gAlgMultAPhi_Cut_[jtAlg][5][rIter] -= signVal;
	gAlgMultAPhi_Cut_[jtAlg][ptIter][rIter] -= signVal;

	break;
      }
    }


    if((AlgJtEta_[jtAlg][0] < 0.6 && AlgJtEta_[jtAlg][1] < 0.6)){
      if(eta > TMath::Min(AlgJtEta_[jtAlg][0], AlgJtEta_[jtAlg][1]) - 0.3){
	for(Int_t rIter = 0; rIter < 10; rIter++){
	  if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	    gAlgImbProjAR_FOR_[jtAlg][5][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
	    gAlgImbProjAR_FOR_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
	    break;
	  }
	}
      }
    }
    else if((AlgJtEta_[jtAlg][0] > -0.6 && AlgJtEta_[jtAlg][1] > -0.6)){
      if(eta < TMath::Max(AlgJtEta_[jtAlg][0], AlgJtEta_[jtAlg][1]) + 0.3){
	for(Int_t rIter = 0; rIter < 10; rIter++){
	  if(tempLeadDelR < rBounds[rIter] || tempSubLeadDelR < rBounds[rIter]){
	    gAlgImbProjAR_FOR_[jtAlg][5][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
	    gAlgImbProjAR_FOR_[jtAlg][ptIter][rIter] -= pt*cos(getDPHI(phi, AlgJtAvePhi_[jtAlg], 3));
	    break;
	  }
	}
      }
    }


  }

  return;
}


#endif
