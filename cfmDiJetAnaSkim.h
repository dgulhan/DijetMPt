//=============================================                                 
// Author: Chris McGinn                                                         
//                                                                              
// DiJet Analysis Skim Class (MC)                                                    
//                                                                              
//=============================================  
#ifndef cfmDiJetAnaSkim_h
#define cfmDiJetAnaSkim_h

#include "cfmDiJetIniSkim.h"

TTree* trackTreeAna_p;
TTree* jetTreeAna_p;
TTree* genTreeAna_p;

//Track Tree Variables

//Tracks proj. onto Alg (enum ordered above, w/ corrected in back 5), All, Cone, and NotCone

Float_t rAlgImbProjAF_[6];
Float_t rAlgImbProjA0_1_[6];
Float_t rAlgImbProjA1_2_[6];
Float_t rAlgImbProjA2_4_[6];
Float_t rAlgImbProjA4_8_[6];
Float_t rAlgImbProjA8_100_[6];
Float_t rAlgImbProjACF_[6];
Float_t rAlgImbProjAC0_1_[6];
Float_t rAlgImbProjAC1_2_[6];
Float_t rAlgImbProjAC2_4_[6];
Float_t rAlgImbProjAC4_8_[6];
Float_t rAlgImbProjAC8_100_[6];
Float_t rAlgImbProjANCF_[6];
Float_t rAlgImbProjANC0_1_[6];
Float_t rAlgImbProjANC1_2_[6];
Float_t rAlgImbProjANC2_4_[6];
Float_t rAlgImbProjANC4_8_[6];
Float_t rAlgImbProjANC8_100_[6];

//DelR ProjA

Float_t rAlgImbProjA1CF_[6];
Float_t rAlgImbProjA1C0_1_[6];
Float_t rAlgImbProjA1C1_2_[6];
Float_t rAlgImbProjA1C2_4_[6];
Float_t rAlgImbProjA1C4_8_[6];
Float_t rAlgImbProjA1C8_100_[6];

Float_t rAlgImbProjA2CF_[6];
Float_t rAlgImbProjA2C0_1_[6];
Float_t rAlgImbProjA2C1_2_[6];
Float_t rAlgImbProjA2C2_4_[6];
Float_t rAlgImbProjA2C4_8_[6];
Float_t rAlgImbProjA2C8_100_[6];

Float_t rAlgImbProjA3CF_[6];
Float_t rAlgImbProjA3C0_1_[6];
Float_t rAlgImbProjA3C1_2_[6];
Float_t rAlgImbProjA3C2_4_[6];
Float_t rAlgImbProjA3C4_8_[6];
Float_t rAlgImbProjA3C8_100_[6];

Float_t rAlgImbProjA4CF_[6];
Float_t rAlgImbProjA4C0_1_[6];
Float_t rAlgImbProjA4C1_2_[6];
Float_t rAlgImbProjA4C2_4_[6];
Float_t rAlgImbProjA4C4_8_[6];
Float_t rAlgImbProjA4C8_100_[6];

Float_t rAlgImbProjA5CF_[6];
Float_t rAlgImbProjA5C0_1_[6];
Float_t rAlgImbProjA5C1_2_[6];
Float_t rAlgImbProjA5C2_4_[6];
Float_t rAlgImbProjA5C4_8_[6];
Float_t rAlgImbProjA5C8_100_[6];

Float_t rAlgImbProjA6CF_[6];
Float_t rAlgImbProjA6C0_1_[6];
Float_t rAlgImbProjA6C1_2_[6];
Float_t rAlgImbProjA6C2_4_[6];
Float_t rAlgImbProjA6C4_8_[6];
Float_t rAlgImbProjA6C8_100_[6];

Float_t rAlgImbProjA7CF_[6];
Float_t rAlgImbProjA7C0_1_[6];
Float_t rAlgImbProjA7C1_2_[6];
Float_t rAlgImbProjA7C2_4_[6];
Float_t rAlgImbProjA7C4_8_[6];
Float_t rAlgImbProjA7C8_100_[6];

Float_t rAlgImbProjA8CF_[6];
Float_t rAlgImbProjA8C0_1_[6];
Float_t rAlgImbProjA8C1_2_[6];
Float_t rAlgImbProjA8C2_4_[6];
Float_t rAlgImbProjA8C4_8_[6];
Float_t rAlgImbProjA8C8_100_[6];

Float_t rAlgImbProjA9CF_[6];
Float_t rAlgImbProjA9C0_1_[6];
Float_t rAlgImbProjA9C1_2_[6];
Float_t rAlgImbProjA9C2_4_[6];
Float_t rAlgImbProjA9C4_8_[6];
Float_t rAlgImbProjA9C8_100_[6];

Float_t rAlgImbProjA10CF_[6];
Float_t rAlgImbProjA10C0_1_[6];
Float_t rAlgImbProjA10C1_2_[6];
Float_t rAlgImbProjA10C2_4_[6];
Float_t rAlgImbProjA10C4_8_[6];
Float_t rAlgImbProjA10C8_100_[6];

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

Bool_t eventSet_[3];
Float_t centWeight_[3];
Float_t centWeight_2pi3_[3];
Float_t centWeight_120_5pi6_[3];
Float_t centWeight_hatAll_0_[3];
Float_t centWeight_hatAll_1_[3];

Float_t fullWeight_[3];

Bool_t isQuarkJet_[3];
Bool_t isGluonJet_[3];

Float_t pthatWeight_;

//Jet Set, Array by algorithm, according to enum above

Float_t AlgLeadJtPt_[3];
Float_t AlgLeadJtPhi_[3];
Float_t AlgLeadJtEta_[3];
Float_t AlgSubLeadJtPt_[3];
Float_t AlgSubLeadJtPhi_[3];
Float_t AlgSubLeadJtEta_[3];
Float_t AlgThirdJtPt_[3];
Float_t AlgThirdJtPhi_[3];
Float_t AlgThirdJtEta_[3];
Float_t AlgFourthJtPt_[3];
Float_t AlgFourthJtPhi_[3];
Float_t AlgFourthJtEta_[3];

Float_t AlgJtAvePhi_[3];
Float_t AlgJtDelPhi_[3];
Float_t AlgJtAsymm_[3];

Float_t AlgLeadRefPt_[3];
Float_t AlgLeadRefPhi_[3];
Float_t AlgLeadRefEta_[3];
Float_t AlgSubLeadRefPt_[3];
Float_t AlgSubLeadRefPhi_[3];
Float_t AlgSubLeadRefEta_[3];
Float_t AlgThirdRefPt_[3];
Float_t AlgThirdRefPhi_[3];
Float_t AlgThirdRefEta_[3];
Float_t AlgFourthRefPt_[3];
Float_t AlgFourthRefPhi_[3];
Float_t AlgFourthRefEta_[3];

Float_t AlgRefAvePhi_[3];
Float_t AlgRefDelPhi_[3];
Float_t AlgRefAsymm_[3];

//Gen Tree Variables

//Gen. proj. onto Jets, ordered by algorithm according to enum, PuPF == [0], PuCalo == [1], etc.

Float_t gAlgImbProjAF_[3];
Float_t gAlgImbProjA0_1_[3];
Float_t gAlgImbProjA1_2_[3];
Float_t gAlgImbProjA2_4_[3];
Float_t gAlgImbProjA4_8_[3];
Float_t gAlgImbProjA8_100_[3];
Float_t gAlgImbProjACF_[3];
Float_t gAlgImbProjAC0_1_[3];
Float_t gAlgImbProjAC1_2_[3];
Float_t gAlgImbProjAC2_4_[3];
Float_t gAlgImbProjAC4_8_[3];
Float_t gAlgImbProjAC8_100_[3];
Float_t gAlgImbProjANCF_[3];
Float_t gAlgImbProjANC0_1_[3];
Float_t gAlgImbProjANC1_2_[3];
Float_t gAlgImbProjANC2_4_[3];
Float_t gAlgImbProjANC4_8_[3];
Float_t gAlgImbProjANC8_100_[3];

//truth delRs
//ProjA delR

Float_t gAlgImbProjA1CF_[3];
Float_t gAlgImbProjA1C0_1_[3];
Float_t gAlgImbProjA1C1_2_[3];
Float_t gAlgImbProjA1C2_4_[3];
Float_t gAlgImbProjA1C4_8_[3];
Float_t gAlgImbProjA1C8_100_[3];

Float_t gAlgImbProjA2CF_[3];
Float_t gAlgImbProjA2C0_1_[3];
Float_t gAlgImbProjA2C1_2_[3];
Float_t gAlgImbProjA2C2_4_[3];
Float_t gAlgImbProjA2C4_8_[3];
Float_t gAlgImbProjA2C8_100_[3];

Float_t gAlgImbProjA3CF_[3];
Float_t gAlgImbProjA3C0_1_[3];
Float_t gAlgImbProjA3C1_2_[3];
Float_t gAlgImbProjA3C2_4_[3];
Float_t gAlgImbProjA3C4_8_[3];
Float_t gAlgImbProjA3C8_100_[3];

Float_t gAlgImbProjA4CF_[3];
Float_t gAlgImbProjA4C0_1_[3];
Float_t gAlgImbProjA4C1_2_[3];
Float_t gAlgImbProjA4C2_4_[3];
Float_t gAlgImbProjA4C4_8_[3];
Float_t gAlgImbProjA4C8_100_[3];

Float_t gAlgImbProjA5CF_[3];
Float_t gAlgImbProjA5C0_1_[3];
Float_t gAlgImbProjA5C1_2_[3];
Float_t gAlgImbProjA5C2_4_[3];
Float_t gAlgImbProjA5C4_8_[3];
Float_t gAlgImbProjA5C8_100_[3];

Float_t gAlgImbProjA6CF_[3];
Float_t gAlgImbProjA6C0_1_[3];
Float_t gAlgImbProjA6C1_2_[3];
Float_t gAlgImbProjA6C2_4_[3];
Float_t gAlgImbProjA6C4_8_[3];
Float_t gAlgImbProjA6C8_100_[3];

Float_t gAlgImbProjA7CF_[3];
Float_t gAlgImbProjA7C0_1_[3];
Float_t gAlgImbProjA7C1_2_[3];
Float_t gAlgImbProjA7C2_4_[3];
Float_t gAlgImbProjA7C4_8_[3];
Float_t gAlgImbProjA7C8_100_[3];

Float_t gAlgImbProjA8CF_[3];
Float_t gAlgImbProjA8C0_1_[3];
Float_t gAlgImbProjA8C1_2_[3];
Float_t gAlgImbProjA8C2_4_[3];
Float_t gAlgImbProjA8C4_8_[3];
Float_t gAlgImbProjA8C8_100_[3];

Float_t gAlgImbProjA9CF_[3];
Float_t gAlgImbProjA9C0_1_[3];
Float_t gAlgImbProjA9C1_2_[3];
Float_t gAlgImbProjA9C2_4_[3];
Float_t gAlgImbProjA9C4_8_[3];
Float_t gAlgImbProjA9C8_100_[3];

Float_t gAlgImbProjA10CF_[3];
Float_t gAlgImbProjA10C0_1_[3];
Float_t gAlgImbProjA10C1_2_[3];
Float_t gAlgImbProjA10C2_4_[3];
Float_t gAlgImbProjA10C4_8_[3];
Float_t gAlgImbProjA10C8_100_[3];


void SetAnaBranches(Bool_t montecarlo = false)
{
  //Track Tree Branches

  std::cout << "Branches Set" << std::endl;

  //Tracks proj. onto Alg, ordered according to enum above, All, Cone, and NotCone

  trackTreeAna_p->Branch("rAlgImbProjAF", &rAlgImbProjAF_, "rAlgImbProjAF[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA0_1", &rAlgImbProjA0_1_, "rAlgImbProjA0_1[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA1_2", &rAlgImbProjA1_2_, "rAlgImbProjA1_2[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA2_4", &rAlgImbProjA2_4_, "rAlgImbProjA2_4[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA4_8", &rAlgImbProjA4_8_, "rAlgImbProjA4_8[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA8_100", &rAlgImbProjA8_100_, "rAlgImbProjA8_100[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjACF", &rAlgImbProjACF_, "rAlgImbProjACF[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjAC0_1", &rAlgImbProjAC0_1_, "rAlgImbProjAC0_1[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjAC1_2", &rAlgImbProjAC1_2_, "rAlgImbProjAC1_2[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjAC2_4", &rAlgImbProjAC2_4_, "rAlgImbProjAC2_4[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjAC4_8", &rAlgImbProjAC4_8_, "rAlgImbProjAC4_8[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjAC8_100", &rAlgImbProjAC8_100_, "rAlgImbProjAC8_100[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjANCF", &rAlgImbProjANCF_, "rAlgImbProjANCF[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjANC0_1", &rAlgImbProjANC0_1_, "rAlgImbProjANC0_1[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjANC1_2", &rAlgImbProjANC1_2_, "rAlgImbProjANC1_2[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjANC2_4", &rAlgImbProjANC2_4_, "rAlgImbProjANC2_4[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjANC4_8", &rAlgImbProjANC4_8_, "rAlgImbProjANC4_8[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjANC8_100", &rAlgImbProjANC8_100_, "rAlgImbProjANC8_100[6]/F");

  //ProjA DelRs

  trackTreeAna_p->Branch("rAlgImbProjA1CF", &rAlgImbProjA1CF_, "rAlgImbProjA1CF[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA1C0_1", &rAlgImbProjA1C0_1_, "rAlgImbProjA1C0_1[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA1C1_2", &rAlgImbProjA1C1_2_, "rAlgImbProjA1C1_2[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA1C2_4", &rAlgImbProjA1C2_4_, "rAlgImbProjA1C2_4[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA1C4_8", &rAlgImbProjA1C4_8_, "rAlgImbProjA1C4_8[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA1C8_100", &rAlgImbProjA1C8_100_, "rAlgImbProjA1C8_100[6]/F");

  trackTreeAna_p->Branch("rAlgImbProjA2CF", &rAlgImbProjA2CF_, "rAlgImbProjA2CF[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA2C0_1", &rAlgImbProjA2C0_1_, "rAlgImbProjA2C0_1[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA2C1_2", &rAlgImbProjA2C1_2_, "rAlgImbProjA2C1_2[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA2C2_4", &rAlgImbProjA2C2_4_, "rAlgImbProjA2C2_4[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA2C4_8", &rAlgImbProjA2C4_8_, "rAlgImbProjA2C4_8[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA2C8_100", &rAlgImbProjA2C8_100_, "rAlgImbProjA2C8_100[6]/F");

  trackTreeAna_p->Branch("rAlgImbProjA3CF", &rAlgImbProjA3CF_, "rAlgImbProjA3CF[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA3C0_1", &rAlgImbProjA3C0_1_, "rAlgImbProjA3C0_1[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA3C1_2", &rAlgImbProjA3C1_2_, "rAlgImbProjA3C1_2[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA3C2_4", &rAlgImbProjA3C2_4_, "rAlgImbProjA3C2_4[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA3C4_8", &rAlgImbProjA3C4_8_, "rAlgImbProjA3C4_8[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA3C8_100", &rAlgImbProjA3C8_100_, "rAlgImbProjA3C8_100[6]/F");  

  trackTreeAna_p->Branch("rAlgImbProjA4CF", &rAlgImbProjA4CF_, "rAlgImbProjA4CF[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA4C0_1", &rAlgImbProjA4C0_1_, "rAlgImbProjA4C0_1[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA4C1_2", &rAlgImbProjA4C1_2_, "rAlgImbProjA4C1_2[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA4C2_4", &rAlgImbProjA4C2_4_, "rAlgImbProjA4C2_4[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA4C4_8", &rAlgImbProjA4C4_8_, "rAlgImbProjA4C4_8[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA4C8_100", &rAlgImbProjA4C8_100_, "rAlgImbProjA4C8_100[6]/F");  

  trackTreeAna_p->Branch("rAlgImbProjA5CF", &rAlgImbProjA5CF_, "rAlgImbProjA5CF[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA5C0_1", &rAlgImbProjA5C0_1_, "rAlgImbProjA5C0_1[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA5C1_2", &rAlgImbProjA5C1_2_, "rAlgImbProjA5C1_2[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA5C2_4", &rAlgImbProjA5C2_4_, "rAlgImbProjA5C2_4[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA5C4_8", &rAlgImbProjA5C4_8_, "rAlgImbProjA5C4_8[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA5C8_100", &rAlgImbProjA5C8_100_, "rAlgImbProjA5C8_100[6]/F");  

  trackTreeAna_p->Branch("rAlgImbProjA6CF", &rAlgImbProjA6CF_, "rAlgImbProjA6CF[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA6C0_1", &rAlgImbProjA6C0_1_, "rAlgImbProjA6C0_1[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA6C1_2", &rAlgImbProjA6C1_2_, "rAlgImbProjA6C1_2[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA6C2_4", &rAlgImbProjA6C2_4_, "rAlgImbProjA6C2_4[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA6C4_8", &rAlgImbProjA6C4_8_, "rAlgImbProjA6C4_8[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA6C8_100", &rAlgImbProjA6C8_100_, "rAlgImbProjA6C8_100[6]/F");  

  trackTreeAna_p->Branch("rAlgImbProjA7CF", &rAlgImbProjA7CF_, "rAlgImbProjA7CF[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA7C0_1", &rAlgImbProjA7C0_1_, "rAlgImbProjA7C0_1[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA7C1_2", &rAlgImbProjA7C1_2_, "rAlgImbProjA7C1_2[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA7C2_4", &rAlgImbProjA7C2_4_, "rAlgImbProjA7C2_4[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA7C4_8", &rAlgImbProjA7C4_8_, "rAlgImbProjA7C4_8[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA7C8_100", &rAlgImbProjA7C8_100_, "rAlgImbProjA7C8_100[6]/F");  

  trackTreeAna_p->Branch("rAlgImbProjA8CF", &rAlgImbProjA8CF_, "rAlgImbProjA8CF[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA8C0_1", &rAlgImbProjA8C0_1_, "rAlgImbProjA8C0_1[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA8C1_2", &rAlgImbProjA8C1_2_, "rAlgImbProjA8C1_2[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA8C2_4", &rAlgImbProjA8C2_4_, "rAlgImbProjA8C2_4[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA8C4_8", &rAlgImbProjA8C4_8_, "rAlgImbProjA8C4_8[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA8C8_100", &rAlgImbProjA8C8_100_, "rAlgImbProjA8C8_100[6]/F");  

  trackTreeAna_p->Branch("rAlgImbProjA9CF", &rAlgImbProjA9CF_, "rAlgImbProjA9CF[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA9C0_1", &rAlgImbProjA9C0_1_, "rAlgImbProjA9C0_1[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA9C1_2", &rAlgImbProjA9C1_2_, "rAlgImbProjA9C1_2[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA9C2_4", &rAlgImbProjA9C2_4_, "rAlgImbProjA9C2_4[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA9C4_8", &rAlgImbProjA9C4_8_, "rAlgImbProjA9C4_8[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA9C8_100", &rAlgImbProjA9C8_100_, "rAlgImbProjA9C8_100[6]/F");  

  trackTreeAna_p->Branch("rAlgImbProjA10CF", &rAlgImbProjA10CF_, "rAlgImbProjA10CF[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA10C0_1", &rAlgImbProjA10C0_1_, "rAlgImbProjA10C0_1[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA10C1_2", &rAlgImbProjA10C1_2_, "rAlgImbProjA10C1_2[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA10C2_4", &rAlgImbProjA10C2_4_, "rAlgImbProjA10C2_4[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA10C4_8", &rAlgImbProjA10C4_8_, "rAlgImbProjA10C4_8[6]/F");
  trackTreeAna_p->Branch("rAlgImbProjA10C8_100", &rAlgImbProjA10C8_100_, "rAlgImbProjA10C8_100[6]/F");  

  //Jet Tree Branches

  jetTreeAna_p->Branch("run", &run_, "run/I");
  jetTreeAna_p->Branch("evt", &evt_, "evt/I");
  jetTreeAna_p->Branch("lumi", &lumi_, "lumi/I");
  jetTreeAna_p->Branch("hiBin", &hiBin_, "hiBin/I");

  jetTreeAna_p->Branch("hiEvtPlane", &hiEvtPlane_, "hiEvtPlane/F");
  jetTreeAna_p->Branch("psin", &psin_, "psin/F");

  jetTreeAna_p->Branch("eventSet", &eventSet_, "eventSet[3]/O");
  jetTreeAna_p->Branch("centWeight", &centWeight_, "centWeight[3]/F");
  jetTreeAna_p->Branch("centWeight_2pi3", &centWeight_2pi3_, "centWeight_2pi3[3]/F");
  jetTreeAna_p->Branch("centWeight_120_5pi6", &centWeight_120_5pi6_, "centWeight_120_5pi6[3]/F");
  jetTreeAna_p->Branch("centWeight_hatAll_0", &centWeight_hatAll_0_, "centWeight_hatAll_0[3]/F");
  jetTreeAna_p->Branch("centWeight_hatAll_1", &centWeight_hatAll_1_, "centWeight_hatAll_1[3]/F");

  jetTreeAna_p->Branch("fullWeight", &fullWeight_, "fullWeight[3]/F");

  jetTreeAna_p->Branch("AlgLeadJtPt", &AlgLeadJtPt_, "AlgLeadJtPt[3]/F");
  jetTreeAna_p->Branch("AlgLeadJtPhi", &AlgLeadJtPhi_, "AlgLeadJtPhi[3]/F");
  jetTreeAna_p->Branch("AlgLeadJtEta", &AlgLeadJtEta_, "AlgLeadJtEta[3]/F");
  jetTreeAna_p->Branch("AlgSubLeadJtPt", &AlgSubLeadJtPt_, "AlgSubLeadJtPt[3]/F");
  jetTreeAna_p->Branch("AlgSubLeadJtPhi", &AlgSubLeadJtPhi_, "AlgSubLeadJtPhi[3]/F");
  jetTreeAna_p->Branch("AlgSubLeadJtEta", &AlgSubLeadJtEta_, "AlgSubLeadJtEta[3]/F");
  jetTreeAna_p->Branch("AlgThirdJtPt", &AlgThirdJtPt_, "AlgThirdJtPt[3]/F");
  jetTreeAna_p->Branch("AlgThirdJtPhi", &AlgThirdJtPhi_, "AlgThirdJtPhi[3]/F");
  jetTreeAna_p->Branch("AlgThirdJtEta", &AlgThirdJtEta_, "AlgThirdJtEta[3]/F");
  jetTreeAna_p->Branch("AlgFourthJtPt", &AlgFourthJtPt_, "AlgFourthJtPt[3]/F");
  jetTreeAna_p->Branch("AlgFourthJtPhi", &AlgFourthJtPhi_, "AlgFourthJtPhi[3]/F");
  jetTreeAna_p->Branch("AlgFourthJtEta", &AlgFourthJtEta_, "AlgFourthJtEta[3]/F");
  jetTreeAna_p->Branch("AlgJtAvePhi", &AlgJtAvePhi_, "AlgJtAvePhi[3]/F");
  jetTreeAna_p->Branch("AlgJtDelPhi", &AlgJtDelPhi_, "AlgJtDelPhi[3]/F");
  jetTreeAna_p->Branch("AlgJtAsymm", &AlgJtAsymm_, "AlgJtAsymm[3]/F");

  if(montecarlo){
    jetTreeAna_p->Branch("isQuarkJet", &isQuarkJet_, "isQuarkJet[3]/O");
    jetTreeAna_p->Branch("isGluonJet", &isGluonJet_, "isGluonJet[3]/O");

    jetTreeAna_p->Branch("pthatWeight", &pthatWeight_, "pthatWeight/F");
    jetTreeAna_p->Branch("pthat", &pthat_, "pthat/F");

    //refpt for jets immediately above
    jetTreeAna_p->Branch("AlgLeadRefPt", &AlgLeadRefPt_, "AlgLeadRefPt[3]/F");
    jetTreeAna_p->Branch("AlgLeadRefPhi", &AlgLeadRefPhi_, "AlgLeadRefPhi[3]/F");
    jetTreeAna_p->Branch("AlgLeadRefEta", &AlgLeadRefEta_, "AlgLeadRefEta[3]/F");
    jetTreeAna_p->Branch("AlgSubLeadRefPt", &AlgSubLeadRefPt_, "AlgSubLeadRefPt[3]/F");
    jetTreeAna_p->Branch("AlgSubLeadRefPhi", &AlgSubLeadRefPhi_, "AlgSubLeadRefPhi[3]/F");
    jetTreeAna_p->Branch("AlgSubLeadRefEta", &AlgSubLeadRefEta_, "AlgSubLeadRefEta[3]/F");
    jetTreeAna_p->Branch("AlgThirdRefPt", &AlgThirdRefPt_, "AlgThirdRefPt[3]/F");
    jetTreeAna_p->Branch("AlgThirdRefPhi", &AlgThirdRefPhi_, "AlgThirdRefPhi[3]/F");
    jetTreeAna_p->Branch("AlgThirdRefEta", &AlgThirdRefEta_, "AlgThirdRefEta[3]/F");
    jetTreeAna_p->Branch("AlgFourthRefPt", &AlgFourthRefPt_, "AlgFourthRefPt[3]/F");
    jetTreeAna_p->Branch("AlgFourthRefPhi", &AlgFourthRefPhi_, "AlgFourthRefPhi[3]/F");
    jetTreeAna_p->Branch("AlgFourthRefEta", &AlgFourthRefEta_, "AlgFourthRefEta[3]/F");

    jetTreeAna_p->Branch("AlgRefAvePhi", &AlgRefAvePhi_, "AlgRefAvePhi[3]/F");
    jetTreeAna_p->Branch("AlgRefDelPhi", &AlgRefDelPhi_, "AlgRefDelPhi[3]/F");
    jetTreeAna_p->Branch("AlgRefAsymm", &AlgRefAsymm_, "AlgRefAsymm[3]/F");


    //Gen. proj. onto jetAlg, array ordered according to enum

    genTreeAna_p->Branch("gAlgImbProjAF", &gAlgImbProjAF_, "gAlgImbProjAF[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA0_1", &gAlgImbProjA0_1_, "gAlgImbProjA0_1[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA1_2", &gAlgImbProjA1_2_, "gAlgImbProjA1_2[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA2_4", &gAlgImbProjA2_4_, "gAlgImbProjA2_4[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA4_8", &gAlgImbProjA4_8_, "gAlgImbProjA4_8[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA8_100", &gAlgImbProjA8_100_, "gAlgImbProjA8_100[3]/F");
    genTreeAna_p->Branch("gAlgImbProjACF", &gAlgImbProjACF_, "gAlgImbProjACF[3]/F");
    genTreeAna_p->Branch("gAlgImbProjAC0_1", &gAlgImbProjAC0_1_, "gAlgImbProjAC0_1[3]/F");
    genTreeAna_p->Branch("gAlgImbProjAC1_2", &gAlgImbProjAC1_2_, "gAlgImbProjAC1_2[3]/F");
    genTreeAna_p->Branch("gAlgImbProjAC2_4", &gAlgImbProjAC2_4_, "gAlgImbProjAC2_4[3]/F");
    genTreeAna_p->Branch("gAlgImbProjAC4_8", &gAlgImbProjAC4_8_, "gAlgImbProjAC4_8[3]/F");
    genTreeAna_p->Branch("gAlgImbProjAC8_100", &gAlgImbProjAC8_100_, "gAlgImbProjAC8_100[3]/F");
    genTreeAna_p->Branch("gAlgImbProjANCF", &gAlgImbProjANCF_, "gAlgImbProjANCF[3]/F");
    genTreeAna_p->Branch("gAlgImbProjANC0_1", &gAlgImbProjANC0_1_, "gAlgImbProjANC0_1[3]/F");
    genTreeAna_p->Branch("gAlgImbProjANC1_2", &gAlgImbProjANC1_2_, "gAlgImbProjANC1_2[3]/F");
    genTreeAna_p->Branch("gAlgImbProjANC2_4", &gAlgImbProjANC2_4_, "gAlgImbProjANC2_4[3]/F");
    genTreeAna_p->Branch("gAlgImbProjANC4_8", &gAlgImbProjANC4_8_, "gAlgImbProjANC4_8[3]/F");
    genTreeAna_p->Branch("gAlgImbProjANC8_100", &gAlgImbProjANC8_100_, "gAlgImbProjANC8_100[3]/F");


    //Truth Del Rs, Proj
    //Proj A's

    genTreeAna_p->Branch("gAlgImbProjA1CF", &gAlgImbProjA1CF_, "gAlgImbProjA1CF[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA1C0_1", &gAlgImbProjA1C0_1_, "gAlgImbProjA1C0_1[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA1C1_2", &gAlgImbProjA1C1_2_, "gAlgImbProjA1C1_2[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA1C2_4", &gAlgImbProjA1C2_4_, "gAlgImbProjA1C2_4[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA1C4_8", &gAlgImbProjA1C4_8_, "gAlgImbProjA1C4_8[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA1C8_100", &gAlgImbProjA1C8_100_, "gAlgImbProjA1C8_100[3]/F");

    genTreeAna_p->Branch("gAlgImbProjA2CF", &gAlgImbProjA2CF_, "gAlgImbProjA2CF[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA2C0_1", &gAlgImbProjA2C0_1_, "gAlgImbProjA2C0_1[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA2C1_2", &gAlgImbProjA2C1_2_, "gAlgImbProjA2C1_2[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA2C2_4", &gAlgImbProjA2C2_4_, "gAlgImbProjA2C2_4[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA2C4_8", &gAlgImbProjA2C4_8_, "gAlgImbProjA2C4_8[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA2C8_100", &gAlgImbProjA2C8_100_, "gAlgImbProjA2C8_100[3]/F");

    genTreeAna_p->Branch("gAlgImbProjA3CF", &gAlgImbProjA3CF_, "gAlgImbProjA3CF[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA3C0_1", &gAlgImbProjA3C0_1_, "gAlgImbProjA3C0_1[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA3C1_2", &gAlgImbProjA3C1_2_, "gAlgImbProjA3C1_2[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA3C2_4", &gAlgImbProjA3C2_4_, "gAlgImbProjA3C2_4[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA3C4_8", &gAlgImbProjA3C4_8_, "gAlgImbProjA3C4_8[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA3C8_100", &gAlgImbProjA3C8_100_, "gAlgImbProjA3C8_100[3]/F");

    genTreeAna_p->Branch("gAlgImbProjA4CF", &gAlgImbProjA4CF_, "gAlgImbProjA4CF[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA4C0_1", &gAlgImbProjA4C0_1_, "gAlgImbProjA4C0_1[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA4C1_2", &gAlgImbProjA4C1_2_, "gAlgImbProjA4C1_2[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA4C2_4", &gAlgImbProjA4C2_4_, "gAlgImbProjA4C2_4[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA4C4_8", &gAlgImbProjA4C4_8_, "gAlgImbProjA4C4_8[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA4C8_100", &gAlgImbProjA4C8_100_, "gAlgImbProjA4C8_100[3]/F");

    genTreeAna_p->Branch("gAlgImbProjA5CF", &gAlgImbProjA5CF_, "gAlgImbProjA5CF[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA5C0_1", &gAlgImbProjA5C0_1_, "gAlgImbProjA5C0_1[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA5C1_2", &gAlgImbProjA5C1_2_, "gAlgImbProjA5C1_2[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA5C2_4", &gAlgImbProjA5C2_4_, "gAlgImbProjA5C2_4[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA5C4_8", &gAlgImbProjA5C4_8_, "gAlgImbProjA5C4_8[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA5C8_100", &gAlgImbProjA5C8_100_, "gAlgImbProjA5C8_100[3]/F");

    genTreeAna_p->Branch("gAlgImbProjA6CF", &gAlgImbProjA6CF_, "gAlgImbProjA6CF[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA6C0_1", &gAlgImbProjA6C0_1_, "gAlgImbProjA6C0_1[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA6C1_2", &gAlgImbProjA6C1_2_, "gAlgImbProjA6C1_2[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA6C2_4", &gAlgImbProjA6C2_4_, "gAlgImbProjA6C2_4[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA6C4_8", &gAlgImbProjA6C4_8_, "gAlgImbProjA6C4_8[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA6C8_100", &gAlgImbProjA6C8_100_, "gAlgImbProjA6C8_100[3]/F");

    genTreeAna_p->Branch("gAlgImbProjA7CF", &gAlgImbProjA7CF_, "gAlgImbProjA7CF[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA7C0_1", &gAlgImbProjA7C0_1_, "gAlgImbProjA7C0_1[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA7C1_2", &gAlgImbProjA7C1_2_, "gAlgImbProjA7C1_2[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA7C2_4", &gAlgImbProjA7C2_4_, "gAlgImbProjA7C2_4[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA7C4_8", &gAlgImbProjA7C4_8_, "gAlgImbProjA7C4_8[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA7C8_100", &gAlgImbProjA7C8_100_, "gAlgImbProjA7C8_100[3]/F");

    genTreeAna_p->Branch("gAlgImbProjA8CF", &gAlgImbProjA8CF_, "gAlgImbProjA8CF[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA8C0_1", &gAlgImbProjA8C0_1_, "gAlgImbProjA8C0_1[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA8C1_2", &gAlgImbProjA8C1_2_, "gAlgImbProjA8C1_2[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA8C2_4", &gAlgImbProjA8C2_4_, "gAlgImbProjA8C2_4[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA8C4_8", &gAlgImbProjA8C4_8_, "gAlgImbProjA8C4_8[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA8C8_100", &gAlgImbProjA8C8_100_, "gAlgImbProjA8C8_100[3]/F");

    genTreeAna_p->Branch("gAlgImbProjA9CF", &gAlgImbProjA9CF_, "gAlgImbProjA9CF[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA9C0_1", &gAlgImbProjA9C0_1_, "gAlgImbProjA9C0_1[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA9C1_2", &gAlgImbProjA9C1_2_, "gAlgImbProjA9C1_2[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA9C2_4", &gAlgImbProjA9C2_4_, "gAlgImbProjA9C2_4[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA9C4_8", &gAlgImbProjA9C4_8_, "gAlgImbProjA9C4_8[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA9C8_100", &gAlgImbProjA9C8_100_, "gAlgImbProjA9C8_100[3]/F");

    genTreeAna_p->Branch("gAlgImbProjA10CF", &gAlgImbProjA10CF_, "gAlgImbProjA10CF[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA10C0_1", &gAlgImbProjA10C0_1_, "gAlgImbProjA10C0_1[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA10C1_2", &gAlgImbProjA10C1_2_, "gAlgImbProjA10C1_2[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA10C2_4", &gAlgImbProjA10C2_4_, "gAlgImbProjA10C2_4[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA10C4_8", &gAlgImbProjA10C4_8_, "gAlgImbProjA10C4_8[3]/F");
    genTreeAna_p->Branch("gAlgImbProjA10C8_100", &gAlgImbProjA10C8_100_, "gAlgImbProjA10C8_100[3]/F");
  }
}


void InitDiJetAnaSkim(Bool_t montecarlo = false)
{
  std::cout << "Init DiJet AnaSkim" << std::endl;

  trackTreeAna_p = new TTree("trackTreeAna", "trackTreeAna");
  jetTreeAna_p = new TTree("jetTreeAna", "jetTreeAna");

  if(montecarlo)
    genTreeAna_p = new TTree("genTreeAna", "genTreeAna");

  SetAnaBranches(montecarlo);
}


void InitJetVar(Bool_t montecarlo = false)
{
  for(Int_t initIter = 0; initIter < 3; initIter++){
    eventSet_[initIter] = false;
    centWeight_[initIter] = -1;
    centWeight_2pi3_[initIter] = -1;
    centWeight_120_5pi6_[initIter] = -1;
    centWeight_hatAll_0_[initIter] = -1;
    centWeight_hatAll_1_[initIter] = -1;

    fullWeight_[initIter] = -1;

    AlgLeadJtPt_[initIter] = -10;
    AlgSubLeadJtPt_[initIter] = -10;
    AlgThirdJtPt_[initIter] = -10;
    AlgFourthJtPt_[initIter] = -10;
    AlgLeadJtPhi_[initIter] = -10;
    AlgSubLeadJtPhi_[initIter] = -10;
    AlgThirdJtPhi_[initIter] = -10;
    AlgFourthJtPhi_[initIter] = -10;
    AlgLeadJtEta_[initIter] = -10;
    AlgSubLeadJtEta_[initIter] = -10;
    AlgThirdJtEta_[initIter] = -10;
    AlgFourthJtEta_[initIter] = -10;

    AlgJtAvePhi_[initIter] = -10;
    AlgJtDelPhi_[initIter] = -10;
    AlgJtAsymm_[initIter] = -10;

    if(montecarlo){
      isQuarkJet_[initIter] = false;
      isGluonJet_[initIter] = false;

      pthatWeight_ = -10;

      AlgLeadRefPt_[initIter] = -10;
      AlgSubLeadRefPt_[initIter] = -10;
      AlgThirdRefPt_[initIter] = -10;
      AlgFourthRefPt_[initIter] = -10;

      AlgLeadRefPhi_[initIter] = -10;
      AlgSubLeadRefPhi_[initIter] = -10;
      AlgThirdRefPhi_[initIter] = -10;
      AlgFourthRefPhi_[initIter] = -10;

      AlgLeadRefEta_[initIter] = -10;
      AlgSubLeadRefEta_[initIter] = -10;
      AlgThirdRefEta_[initIter] = -10;
      AlgFourthRefEta_[initIter] = -10;

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
    rAlgImbProjAF_[initIter] = 0;
    rAlgImbProjA0_1_[initIter] = 0;
    rAlgImbProjA1_2_[initIter] = 0;
    rAlgImbProjA2_4_[initIter] = 0;
    rAlgImbProjA4_8_[initIter] = 0;
    rAlgImbProjA8_100_[initIter] = 0;
    rAlgImbProjACF_[initIter] = 0;
    rAlgImbProjAC0_1_[initIter] = 0;
    rAlgImbProjAC1_2_[initIter] = 0;
    rAlgImbProjAC2_4_[initIter] = 0;
    rAlgImbProjAC4_8_[initIter] = 0;
    rAlgImbProjAC8_100_[initIter] = 0;
    rAlgImbProjANCF_[initIter] = 0;
    rAlgImbProjANC0_1_[initIter] = 0;
    rAlgImbProjANC1_2_[initIter] = 0;
    rAlgImbProjANC2_4_[initIter] = 0;
    rAlgImbProjANC4_8_[initIter] = 0;
    rAlgImbProjANC8_100_[initIter] = 0;

    //DelRs Proj

    //ProjA

    rAlgImbProjA1CF_[initIter] = 0;
    rAlgImbProjA1C0_1_[initIter] = 0;
    rAlgImbProjA1C1_2_[initIter] = 0;
    rAlgImbProjA1C2_4_[initIter] = 0;
    rAlgImbProjA1C4_8_[initIter] = 0;
    rAlgImbProjA1C8_100_[initIter] = 0;

    rAlgImbProjA2CF_[initIter] = 0;
    rAlgImbProjA2C0_1_[initIter] = 0;
    rAlgImbProjA2C1_2_[initIter] = 0;
    rAlgImbProjA2C2_4_[initIter] = 0;
    rAlgImbProjA2C4_8_[initIter] = 0;
    rAlgImbProjA2C8_100_[initIter] = 0;

    rAlgImbProjA3CF_[initIter] = 0;
    rAlgImbProjA3C0_1_[initIter] = 0;
    rAlgImbProjA3C1_2_[initIter] = 0;
    rAlgImbProjA3C2_4_[initIter] = 0;
    rAlgImbProjA3C4_8_[initIter] = 0;
    rAlgImbProjA3C8_100_[initIter] = 0;

    rAlgImbProjA4CF_[initIter] = 0;
    rAlgImbProjA4C0_1_[initIter] = 0;
    rAlgImbProjA4C1_2_[initIter] = 0;
    rAlgImbProjA4C2_4_[initIter] = 0;
    rAlgImbProjA4C4_8_[initIter] = 0;
    rAlgImbProjA4C8_100_[initIter] = 0;

    rAlgImbProjA5CF_[initIter] = 0;
    rAlgImbProjA5C0_1_[initIter] = 0;
    rAlgImbProjA5C1_2_[initIter] = 0;
    rAlgImbProjA5C2_4_[initIter] = 0;
    rAlgImbProjA5C4_8_[initIter] = 0;
    rAlgImbProjA5C8_100_[initIter] = 0;

    rAlgImbProjA6CF_[initIter] = 0;
    rAlgImbProjA6C0_1_[initIter] = 0;
    rAlgImbProjA6C1_2_[initIter] = 0;
    rAlgImbProjA6C2_4_[initIter] = 0;
    rAlgImbProjA6C4_8_[initIter] = 0;
    rAlgImbProjA6C8_100_[initIter] = 0;

    rAlgImbProjA7CF_[initIter] = 0;
    rAlgImbProjA7C0_1_[initIter] = 0;
    rAlgImbProjA7C1_2_[initIter] = 0;
    rAlgImbProjA7C2_4_[initIter] = 0;
    rAlgImbProjA7C4_8_[initIter] = 0;
    rAlgImbProjA7C8_100_[initIter] = 0;

    rAlgImbProjA8CF_[initIter] = 0;
    rAlgImbProjA8C0_1_[initIter] = 0;
    rAlgImbProjA8C1_2_[initIter] = 0;
    rAlgImbProjA8C2_4_[initIter] = 0;
    rAlgImbProjA8C4_8_[initIter] = 0;
    rAlgImbProjA8C8_100_[initIter] = 0;

    rAlgImbProjA9CF_[initIter] = 0;
    rAlgImbProjA9C0_1_[initIter] = 0;
    rAlgImbProjA9C1_2_[initIter] = 0;
    rAlgImbProjA9C2_4_[initIter] = 0;
    rAlgImbProjA9C4_8_[initIter] = 0;
    rAlgImbProjA9C8_100_[initIter] = 0;

    rAlgImbProjA10CF_[initIter] = 0;
    rAlgImbProjA10C0_1_[initIter] = 0;
    rAlgImbProjA10C1_2_[initIter] = 0;
    rAlgImbProjA10C2_4_[initIter] = 0;
    rAlgImbProjA10C4_8_[initIter] = 0;
    rAlgImbProjA10C8_100_[initIter] = 0;
  }

  if(montecarlo){
    //Gen. proj. onto Truth
    for(Int_t initIter = 0; initIter < 3; initIter++){

      gAlgImbProjAF_[initIter] = 0;
      gAlgImbProjA0_1_[initIter] = 0;
      gAlgImbProjA1_2_[initIter] = 0;
      gAlgImbProjA2_4_[initIter] = 0;
      gAlgImbProjA4_8_[initIter] = 0;
      gAlgImbProjA8_100_[initIter] = 0;
      gAlgImbProjACF_[initIter] = 0;
      gAlgImbProjAC0_1_[initIter] = 0;
      gAlgImbProjAC1_2_[initIter] = 0;
      gAlgImbProjAC2_4_[initIter] = 0;
      gAlgImbProjAC4_8_[initIter] = 0;
      gAlgImbProjAC8_100_[initIter] = 0;
      gAlgImbProjANCF_[initIter] = 0;
      gAlgImbProjANC0_1_[initIter] = 0;
      gAlgImbProjANC1_2_[initIter] = 0;
      gAlgImbProjANC2_4_[initIter] = 0;
      gAlgImbProjANC4_8_[initIter] = 0;
      gAlgImbProjANC8_100_[initIter] = 0;

      //DelRs
      //proj As

      gAlgImbProjA1CF_[initIter] = 0;
      gAlgImbProjA1C0_1_[initIter] = 0;
      gAlgImbProjA1C1_2_[initIter] = 0;
      gAlgImbProjA1C2_4_[initIter] = 0;
      gAlgImbProjA1C4_8_[initIter] = 0;
      gAlgImbProjA1C8_100_[initIter] = 0;

      gAlgImbProjA2CF_[initIter] = 0;
      gAlgImbProjA2C0_1_[initIter] = 0;
      gAlgImbProjA2C1_2_[initIter] = 0;
      gAlgImbProjA2C2_4_[initIter] = 0;
      gAlgImbProjA2C4_8_[initIter] = 0;
      gAlgImbProjA2C8_100_[initIter] = 0;

      gAlgImbProjA3CF_[initIter] = 0;
      gAlgImbProjA3C0_1_[initIter] = 0;
      gAlgImbProjA3C1_2_[initIter] = 0;
      gAlgImbProjA3C2_4_[initIter] = 0;
      gAlgImbProjA3C4_8_[initIter] = 0;
      gAlgImbProjA3C8_100_[initIter] = 0;

      gAlgImbProjA4CF_[initIter] = 0;
      gAlgImbProjA4C0_1_[initIter] = 0;
      gAlgImbProjA4C1_2_[initIter] = 0;
      gAlgImbProjA4C2_4_[initIter] = 0;
      gAlgImbProjA4C4_8_[initIter] = 0;
      gAlgImbProjA4C8_100_[initIter] = 0;

      gAlgImbProjA5CF_[initIter] = 0;
      gAlgImbProjA5C0_1_[initIter] = 0;
      gAlgImbProjA5C1_2_[initIter] = 0;
      gAlgImbProjA5C2_4_[initIter] = 0;
      gAlgImbProjA5C4_8_[initIter] = 0;
      gAlgImbProjA5C8_100_[initIter] = 0;

      gAlgImbProjA6CF_[initIter] = 0;
      gAlgImbProjA6C0_1_[initIter] = 0;
      gAlgImbProjA6C1_2_[initIter] = 0;
      gAlgImbProjA6C2_4_[initIter] = 0;
      gAlgImbProjA6C4_8_[initIter] = 0;
      gAlgImbProjA6C8_100_[initIter] = 0;

      gAlgImbProjA7CF_[initIter] = 0;
      gAlgImbProjA7C0_1_[initIter] = 0;
      gAlgImbProjA7C1_2_[initIter] = 0;
      gAlgImbProjA7C2_4_[initIter] = 0;
      gAlgImbProjA7C4_8_[initIter] = 0;
      gAlgImbProjA7C8_100_[initIter] = 0;

      gAlgImbProjA8CF_[initIter] = 0;
      gAlgImbProjA8C0_1_[initIter] = 0;
      gAlgImbProjA8C1_2_[initIter] = 0;
      gAlgImbProjA8C2_4_[initIter] = 0;
      gAlgImbProjA8C4_8_[initIter] = 0;
      gAlgImbProjA8C8_100_[initIter] = 0;

      gAlgImbProjA9CF_[initIter] = 0;
      gAlgImbProjA9C0_1_[initIter] = 0;
      gAlgImbProjA9C1_2_[initIter] = 0;
      gAlgImbProjA9C2_4_[initIter] = 0;
      gAlgImbProjA9C4_8_[initIter] = 0;
      gAlgImbProjA9C8_100_[initIter] = 0;

      gAlgImbProjA10CF_[initIter] = 0;
      gAlgImbProjA10C0_1_[initIter] = 0;
      gAlgImbProjA10C1_2_[initIter] = 0;
      gAlgImbProjA10C2_4_[initIter] = 0;
      gAlgImbProjA10C4_8_[initIter] = 0;
      gAlgImbProjA10C8_100_[initIter] = 0;
    }    
  }

}


#endif
