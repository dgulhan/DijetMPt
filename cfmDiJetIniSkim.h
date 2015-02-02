//=============================================                                 
// Author: Chris McGinn                                                         
//                                                                              
// DiJet Initial Skim Class (MC)                                                    
//                                                                              
//=============================================  
#ifndef cfmDiJetIniSkim_h
#define cfmDiJetIniSkim_h

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include "/net/hisrv0001/home/cfmcginn/emDiJet/CMSSW_5_3_12_patch3/tempHIFA/HiForestAnalysis/commonSetup.h"
#include "TMath.h"
#include "factorizedPtCorr.h"
#include "effCorrNPF.h"
#include "jecPtCorr.h"
#include "resPtCorr.h"

enum AlgoType_PbPb{
  Pu3Calo,     //0
  Pu4Calo,     //1
  Pu5Calo,     //2
  Vs2Calo,     //3
  Vs3Calo,     //4
  Vs4Calo,     //5
  Vs5Calo,     //6
  Vs2CaloFrag, //7
  Vs3CaloFrag, //8
  Vs4CaloFrag, //9
  Vs5CaloFrag, //10
  Vs2CaloRes,  //11
  Vs3CaloRes,  //12
  Vs4CaloRes,  //13
  Vs5CaloRes,  //14
  T2,          //15
  T3,          //16
  T4,          //17
  T5,          //18
  PuPF,        //19
  VsPF         //20
};


Bool_t isMonteCarlo(sampleType sType = kHIDATA){
  if(sType == kHIMC || sType == kPPMC || sType == kPAMC) return true;
  else return false;
}


Bool_t isHI(sampleType sType = kHIDATA){
  if(sType == kHIDATA || sType == kHIMC) return true;
  else return false;
}


TString getSampleName ( sampleType colli) {
  if (colli == kHIDATA) return "pbpbDATA";
  if (colli == kHIMC) return "pbpbMC";
  if (colli == kPPDATA) return "ppDATA";
  if (colli == kPPMC) return "ppMC";
  if (colli == kPADATA) return "ppbDATA";
  if (colli == kPAMC) return "ppbMC";
  return "NULL";
}
TString getSampleName ( int colli) {
  if (colli == kHIDATA) return "pbpbDATA";
  if (colli == kHIMC) return "pbpbMC";
  if (colli == kPPDATA) return "ppDATA";
  if (colli == kPPMC) return "ppMC";
  if (colli == kPADATA) return "ppbDATA";
  if (colli == kPAMC) return "ppbMC";
  return "NULL";
}

TTree* trackTreeIni_p;
TTree* pfCandTreeIni_p;
TTree* jetTreeIni_p;
TTree* genTreeIni_p;

//Track Tree Variables

Int_t nTrk_;
Float_t trkPt_[maxTracks];
Float_t trkPhi_[maxTracks];
Float_t trkEta_[maxTracks];

//PFCand Tree Variables

Int_t nPF_;
Float_t pfPt_[maxPF];
Float_t pfPhi_[maxPF];
Float_t pfEta_[maxPF];

//Jet Tree Variables

Int_t runIni_;
Int_t evtIni_;
Int_t lumiIni_;
Int_t hiBinIni_;
Float_t pthatIni_;

Float_t hiEvtPlaneIni_;
Float_t psinIni_;

Int_t nPu3Calo_;
Float_t Pu3CaloPt_[maxJets];
Float_t Pu3CaloPhi_[maxJets];
Float_t Pu3CaloEta_[maxJets];
Float_t Pu3CaloTrkMax_[maxJets];
Float_t Pu3CaloRawPt_[maxJets];
Float_t Pu3CaloRefPt_[maxJets];
Float_t Pu3CaloRefPhi_[maxJets];
Float_t Pu3CaloRefEta_[maxJets];
Int_t Pu3CaloRefPart_[maxJets];

Int_t nPu4Calo_;
Float_t Pu4CaloPt_[maxJets];
Float_t Pu4CaloPhi_[maxJets];
Float_t Pu4CaloEta_[maxJets];
Float_t Pu4CaloTrkMax_[maxJets];
Float_t Pu4CaloRawPt_[maxJets];
Float_t Pu4CaloRefPt_[maxJets];
Float_t Pu4CaloRefPhi_[maxJets];
Float_t Pu4CaloRefEta_[maxJets];
Int_t Pu4CaloRefPart_[maxJets];

Int_t nPu5Calo_;
Float_t Pu5CaloPt_[maxJets];
Float_t Pu5CaloPhi_[maxJets];
Float_t Pu5CaloEta_[maxJets];
Float_t Pu5CaloTrkMax_[maxJets];
Float_t Pu5CaloRawPt_[maxJets];
Float_t Pu5CaloRefPt_[maxJets];
Float_t Pu5CaloRefPhi_[maxJets];
Float_t Pu5CaloRefEta_[maxJets];
Int_t Pu5CaloRefPart_[maxJets];

Int_t nVs2Calo_;
Float_t Vs2CaloPt_[maxJets];
Float_t Vs2CaloPhi_[maxJets];
Float_t Vs2CaloEta_[maxJets];
Float_t Vs2CaloTrkMax_[maxJets];
Float_t Vs2CaloRawPt_[maxJets];
Float_t Vs2CaloRefPt_[maxJets];
Float_t Vs2CaloRefPhi_[maxJets];
Float_t Vs2CaloRefEta_[maxJets];
Int_t Vs2CaloRefPart_[maxJets];

Int_t nVs3Calo_;
Float_t Vs3CaloPt_[maxJets];
Float_t Vs3CaloPhi_[maxJets];
Float_t Vs3CaloEta_[maxJets];
Float_t Vs3CaloTrkMax_[maxJets];
Float_t Vs3CaloRawPt_[maxJets];
Float_t Vs3CaloRefPt_[maxJets];
Float_t Vs3CaloRefPhi_[maxJets];
Float_t Vs3CaloRefEta_[maxJets];
Int_t Vs3CaloRefPart_[maxJets];

Int_t nVs4Calo_;
Float_t Vs4CaloPt_[maxJets];
Float_t Vs4CaloPhi_[maxJets];
Float_t Vs4CaloEta_[maxJets];
Float_t Vs4CaloTrkMax_[maxJets];
Float_t Vs4CaloRawPt_[maxJets];
Float_t Vs4CaloRefPt_[maxJets];
Float_t Vs4CaloRefPhi_[maxJets];
Float_t Vs4CaloRefEta_[maxJets];
Int_t Vs4CaloRefPart_[maxJets];

Int_t nVs5Calo_;
Float_t Vs5CaloPt_[maxJets];
Float_t Vs5CaloPhi_[maxJets];
Float_t Vs5CaloEta_[maxJets];
Float_t Vs5CaloTrkMax_[maxJets];
Float_t Vs5CaloRawPt_[maxJets];
Float_t Vs5CaloRefPt_[maxJets];
Float_t Vs5CaloRefPhi_[maxJets];
Float_t Vs5CaloRefEta_[maxJets];
Int_t Vs5CaloRefPart_[maxJets];


Int_t nVs2CaloFrag_;
Float_t Vs2CaloFragPt_[maxJets];
Float_t Vs2CaloFragPhi_[maxJets];
Float_t Vs2CaloFragEta_[maxJets];
Float_t Vs2CaloFragTrkMax_[maxJets];
Float_t Vs2CaloFragRawPt_[maxJets];
Float_t Vs2CaloFragRefPt_[maxJets];
Float_t Vs2CaloFragRefPhi_[maxJets];
Float_t Vs2CaloFragRefEta_[maxJets];
Int_t Vs2CaloFragRefPart_[maxJets];

Int_t nVs3CaloFrag_;
Float_t Vs3CaloFragPt_[maxJets];
Float_t Vs3CaloFragPhi_[maxJets];
Float_t Vs3CaloFragEta_[maxJets];
Float_t Vs3CaloFragTrkMax_[maxJets];
Float_t Vs3CaloFragRawPt_[maxJets];
Float_t Vs3CaloFragRefPt_[maxJets];
Float_t Vs3CaloFragRefPhi_[maxJets];
Float_t Vs3CaloFragRefEta_[maxJets];
Int_t Vs3CaloFragRefPart_[maxJets];

Int_t nVs4CaloFrag_;
Float_t Vs4CaloFragPt_[maxJets];
Float_t Vs4CaloFragPhi_[maxJets];
Float_t Vs4CaloFragEta_[maxJets];
Float_t Vs4CaloFragTrkMax_[maxJets];
Float_t Vs4CaloFragRawPt_[maxJets];
Float_t Vs4CaloFragRefPt_[maxJets];
Float_t Vs4CaloFragRefPhi_[maxJets];
Float_t Vs4CaloFragRefEta_[maxJets];
Int_t Vs4CaloFragRefPart_[maxJets];

Int_t nVs5CaloFrag_;
Float_t Vs5CaloFragPt_[maxJets];
Float_t Vs5CaloFragPhi_[maxJets];
Float_t Vs5CaloFragEta_[maxJets];
Float_t Vs5CaloFragTrkMax_[maxJets];
Float_t Vs5CaloFragRawPt_[maxJets];
Float_t Vs5CaloFragRefPt_[maxJets];
Float_t Vs5CaloFragRefPhi_[maxJets];
Float_t Vs5CaloFragRefEta_[maxJets];
Int_t Vs5CaloFragRefPart_[maxJets];


Int_t nVs2CaloRes_;
Float_t Vs2CaloResPt_[maxJets];
Float_t Vs2CaloResPhi_[maxJets];
Float_t Vs2CaloResEta_[maxJets];
Float_t Vs2CaloResTrkMax_[maxJets];
Float_t Vs2CaloResRawPt_[maxJets];
Float_t Vs2CaloResRefPt_[maxJets];
Float_t Vs2CaloResRefPhi_[maxJets];
Float_t Vs2CaloResRefEta_[maxJets];
Int_t Vs2CaloResRefPart_[maxJets];


Int_t nVs3CaloRes_;
Float_t Vs3CaloResPt_[maxJets];
Float_t Vs3CaloResPhi_[maxJets];
Float_t Vs3CaloResEta_[maxJets];
Float_t Vs3CaloResTrkMax_[maxJets];
Float_t Vs3CaloResRawPt_[maxJets];
Float_t Vs3CaloResRefPt_[maxJets];
Float_t Vs3CaloResRefPhi_[maxJets];
Float_t Vs3CaloResRefEta_[maxJets];
Int_t Vs3CaloResRefPart_[maxJets];


Int_t nVs4CaloRes_;
Float_t Vs4CaloResPt_[maxJets];
Float_t Vs4CaloResPhi_[maxJets];
Float_t Vs4CaloResEta_[maxJets];
Float_t Vs4CaloResTrkMax_[maxJets];
Float_t Vs4CaloResRawPt_[maxJets];
Float_t Vs4CaloResRefPt_[maxJets];
Float_t Vs4CaloResRefPhi_[maxJets];
Float_t Vs4CaloResRefEta_[maxJets];
Int_t Vs4CaloResRefPart_[maxJets];


Int_t nVs5CaloRes_;
Float_t Vs5CaloResPt_[maxJets];
Float_t Vs5CaloResPhi_[maxJets];
Float_t Vs5CaloResEta_[maxJets];
Float_t Vs5CaloResTrkMax_[maxJets];
Float_t Vs5CaloResRawPt_[maxJets];
Float_t Vs5CaloResRefPt_[maxJets];
Float_t Vs5CaloResRefPhi_[maxJets];
Float_t Vs5CaloResRefEta_[maxJets];
Int_t Vs5CaloResRefPart_[maxJets];


Int_t nT2_;
Float_t T2Pt_[maxJets];
Float_t T2Phi_[maxJets];
Float_t T2Eta_[maxJets];
Int_t T2Part_[maxJets];

Int_t nT3_;
Float_t T3Pt_[maxJets];
Float_t T3Phi_[maxJets];
Float_t T3Eta_[maxJets];
Int_t T3Part_[maxJets];

Int_t nT4_;
Float_t T4Pt_[maxJets];
Float_t T4Phi_[maxJets];
Float_t T4Eta_[maxJets];
Int_t T4Part_[maxJets];

Int_t nT5_;
Float_t T5Pt_[maxJets];
Float_t T5Phi_[maxJets];
Float_t T5Eta_[maxJets];
Int_t T5Part_[maxJets];

Int_t nPu3PF_;
Float_t Pu3PFPt_[maxJets];
Float_t Pu3PFPhi_[maxJets];
Float_t Pu3PFEta_[maxJets];
Float_t Pu3PFTrkMax_[maxJets];
Float_t Pu3PFRawPt_[maxJets];
Float_t Pu3PFRefPt_[maxJets];
Float_t Pu3PFRefPhi_[maxJets];
Float_t Pu3PFRefEta_[maxJets];
Int_t Pu3PFRefPart_[maxJets];

Int_t nVs3PF_;
Float_t Vs3PFPt_[maxJets];
Float_t Vs3PFPhi_[maxJets];
Float_t Vs3PFEta_[maxJets];
Float_t Vs3PFTrkMax_[maxJets];
Float_t Vs3PFRawPt_[maxJets];
Float_t Vs3PFRefPt_[maxJets];
Float_t Vs3PFRefPhi_[maxJets];
Float_t Vs3PFRefEta_[maxJets];
Int_t Vs3PFRefPart_[maxJets];

Float_t TrkJtPt_[5];
Float_t TrkJtPhi_[5];
Float_t TrkJtEta_[5];

Int_t nLeadJtConst_;
Float_t TrkLeadJtConstPt_[maxTracks];
Float_t TrkLeadJtConstPhi_[maxTracks];
Float_t TrkLeadJtConstEta_[maxTracks];
Float_t TrkLeadJtConstCorr_[maxTracks];

Int_t nSubLeadJtConst_;
Float_t TrkSubLeadJtConstPt_[maxTracks];
Float_t TrkSubLeadJtConstPhi_[maxTracks];
Float_t TrkSubLeadJtConstEta_[maxTracks];
Float_t TrkSubLeadJtConstCorr_[maxTracks];

Int_t nThirdJtConst_;
Float_t TrkThirdJtConstPt_[maxTracks];
Float_t TrkThirdJtConstPhi_[maxTracks];
Float_t TrkThirdJtConstEta_[maxTracks];
Float_t TrkThirdJtConstCorr_[maxTracks];

Int_t nFourthJtConst_;
Float_t TrkFourthJtConstPt_[maxTracks];
Float_t TrkFourthJtConstPhi_[maxTracks];
Float_t TrkFourthJtConstEta_[maxTracks];
Float_t TrkFourthJtConstCorr_[maxTracks];

Int_t nFifthJtConst_;
Float_t TrkFifthJtConstPt_[maxTracks];
Float_t TrkFifthJtConstPhi_[maxTracks];
Float_t TrkFifthJtConstEta_[maxTracks];
Float_t TrkFifthJtConstCorr_[maxTracks];


//Gen Tree Variables

Int_t nGen_;
Float_t genPt_[maxEntrySim];
Float_t genPhi_[maxEntrySim];
Float_t genEta_[maxEntrySim];

void SetIniBranches(sampleType sType = kHIDATA, Bool_t justJt = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

  std::cout << "Branches Set" << std::endl;

  if(!justJt){
    //Track Tree Branches

    trackTreeIni_p->Branch("nTrk", &nTrk_, "nTrk/I");
  
    trackTreeIni_p->Branch("trkPt", trkPt_, "trkPt[nTrk]/F");
    trackTreeIni_p->Branch("trkPhi", trkPhi_, "trkPhi[nTrk]/F");
    trackTreeIni_p->Branch("trkEta", trkEta_, "trkEta[nTrk]/F");

    //PF Cand Tree Branches

    pfCandTreeIni_p->Branch("nPF", &nPF_, "nPF/I");
    pfCandTreeIni_p->Branch("pfPt", &pfPt_, "pfPt[nPF]/F");
    pfCandTreeIni_p->Branch("pfPhi", &pfPhi_, "pfPhi[nPF]/F");
    pfCandTreeIni_p->Branch("pfEta", &pfEta_, "pfEta[nPF]/F");
  } 

  //Jet Tree Branches

  jetTreeIni_p->Branch("runIni", &runIni_, "runIni/I");
  jetTreeIni_p->Branch("evtIni", &evtIni_, "evtIni/I");
  jetTreeIni_p->Branch("lumiIni", &lumiIni_, "lumiIni/I");

  if(hi)
    jetTreeIni_p->Branch("hiBinIni", &hiBinIni_, "hiBinIni/I");
   
  if(montecarlo)
    jetTreeIni_p->Branch("pthatIni", &pthatIni_, "pthatIni/F");

  if(hi){
    jetTreeIni_p->Branch("hiEvtPlaneIni", &hiEvtPlaneIni_, "hiEvtPlaneIni/F");
    jetTreeIni_p->Branch("psinIni", &psinIni_, "psinIni/F");
  }    

  jetTreeIni_p->Branch("nPu3Calo", &nPu3Calo_, "nPu3Calo/I");
  jetTreeIni_p->Branch("Pu3CaloPt", Pu3CaloPt_, "Pu3CaloPt[nPu3Calo]/F");
  jetTreeIni_p->Branch("Pu3CaloPhi", Pu3CaloPhi_, "Pu3CaloPhi[nPu3Calo]/F");
  jetTreeIni_p->Branch("Pu3CaloEta", Pu3CaloEta_, "Pu3CaloEta[nPu3Calo]/F");
  jetTreeIni_p->Branch("Pu3CaloTrkMax", Pu3CaloTrkMax_, "Pu3CaloTrkMax[nPu3Calo]/F");
  jetTreeIni_p->Branch("Pu3CaloRawPt", Pu3CaloRawPt_, "Pu3CaloRawPt[nPu3Calo]/F");

  jetTreeIni_p->Branch("nPu4Calo", &nPu4Calo_, "nPu4Calo/I");
  jetTreeIni_p->Branch("Pu4CaloPt", Pu4CaloPt_, "Pu4CaloPt[nPu4Calo]/F");
  jetTreeIni_p->Branch("Pu4CaloPhi", Pu4CaloPhi_, "Pu4CaloPhi[nPu4Calo]/F");
  jetTreeIni_p->Branch("Pu4CaloEta", Pu4CaloEta_, "Pu4CaloEta[nPu4Calo]/F");
  jetTreeIni_p->Branch("Pu4CaloTrkMax", Pu4CaloTrkMax_, "Pu4CaloTrkMax[nPu4Calo]/F");
  jetTreeIni_p->Branch("Pu4CaloRawPt", Pu4CaloRawPt_, "Pu4CaloRawPt[nPu4Calo]/F");

  jetTreeIni_p->Branch("nPu5Calo", &nPu5Calo_, "nPu5Calo/I");
  jetTreeIni_p->Branch("Pu5CaloPt", Pu5CaloPt_, "Pu5CaloPt[nPu5Calo]/F");
  jetTreeIni_p->Branch("Pu5CaloPhi", Pu5CaloPhi_, "Pu5CaloPhi[nPu5Calo]/F");
  jetTreeIni_p->Branch("Pu5CaloEta", Pu5CaloEta_, "Pu5CaloEta[nPu5Calo]/F");
  jetTreeIni_p->Branch("Pu5CaloTrkMax", Pu5CaloTrkMax_, "Pu5CaloTrkMax[nPu5Calo]/F");
  jetTreeIni_p->Branch("Pu5CaloRawPt", Pu5CaloRawPt_, "Pu5CaloRawPt[nPu5Calo]/F");

  if(montecarlo){
    jetTreeIni_p->Branch("Pu3CaloRefPt", Pu3CaloRefPt_, "Pu3CaloRefPt[nPu3Calo]/F");
    jetTreeIni_p->Branch("Pu3CaloRefPhi", Pu3CaloRefPhi_, "Pu3CaloRefPhi[nPu3Calo]/F");
    jetTreeIni_p->Branch("Pu3CaloRefEta", Pu3CaloRefEta_, "Pu3CaloRefEta[nPu3Calo]/F");
    jetTreeIni_p->Branch("Pu3CaloRefPart", Pu3CaloRefPart_, "Pu3CaloRefPart[nPu3Calo]/I");

    jetTreeIni_p->Branch("Pu4CaloRefPt", Pu4CaloRefPt_, "Pu4CaloRefPt[nPu4Calo]/F");
    jetTreeIni_p->Branch("Pu4CaloRefPhi", Pu4CaloRefPhi_, "Pu4CaloRefPhi[nPu4Calo]/F");
    jetTreeIni_p->Branch("Pu4CaloRefEta", Pu4CaloRefEta_, "Pu4CaloRefEta[nPu4Calo]/F");
    jetTreeIni_p->Branch("Pu4CaloRefPart", Pu4CaloRefPart_, "Pu4CaloRefPart[nPu4Calo]/I");

    jetTreeIni_p->Branch("Pu5CaloRefPt", Pu5CaloRefPt_, "Pu5CaloRefPt[nPu5Calo]/F");
    jetTreeIni_p->Branch("Pu5CaloRefPhi", Pu5CaloRefPhi_, "Pu5CaloRefPhi[nPu5Calo]/F");
    jetTreeIni_p->Branch("Pu5CaloRefEta", Pu5CaloRefEta_, "Pu5CaloRefEta[nPu5Calo]/F");
    jetTreeIni_p->Branch("Pu5CaloRefPart", Pu5CaloRefPart_, "Pu5CaloRefPart[nPu5Calo]/I");
  }    

  jetTreeIni_p->Branch("nVs2Calo", &nVs2Calo_, "nVs2Calo/I");
  jetTreeIni_p->Branch("Vs2CaloPt", Vs2CaloPt_, "Vs2CaloPt[nVs2Calo]/F");
  jetTreeIni_p->Branch("Vs2CaloPhi", Vs2CaloPhi_, "Vs2CaloPhi[nVs2Calo]/F");
  jetTreeIni_p->Branch("Vs2CaloEta", Vs2CaloEta_, "Vs2CaloEta[nVs2Calo]/F");
  jetTreeIni_p->Branch("Vs2CaloTrkMax", Vs2CaloTrkMax_, "Vs2CaloTrkMax[nVs2Calo]/F");
  jetTreeIni_p->Branch("Vs2CaloRawPt", Vs2CaloRawPt_, "Vs2CaloRawPt[nVs2Calo]/F");

  jetTreeIni_p->Branch("nVs3Calo", &nVs3Calo_, "nVs3Calo/I");
  jetTreeIni_p->Branch("Vs3CaloPt", Vs3CaloPt_, "Vs3CaloPt[nVs3Calo]/F");
  jetTreeIni_p->Branch("Vs3CaloPhi", Vs3CaloPhi_, "Vs3CaloPhi[nVs3Calo]/F");
  jetTreeIni_p->Branch("Vs3CaloEta", Vs3CaloEta_, "Vs3CaloEta[nVs3Calo]/F");
  jetTreeIni_p->Branch("Vs3CaloTrkMax", Vs3CaloTrkMax_, "Vs3CaloTrkMax[nVs3Calo]/F");
  jetTreeIni_p->Branch("Vs3CaloRawPt", Vs3CaloRawPt_, "Vs3CaloRawPt[nVs3Calo]/F");

  jetTreeIni_p->Branch("nVs4Calo", &nVs4Calo_, "nVs4Calo/I");
  jetTreeIni_p->Branch("Vs4CaloPt", Vs4CaloPt_, "Vs4CaloPt[nVs4Calo]/F");
  jetTreeIni_p->Branch("Vs4CaloPhi", Vs4CaloPhi_, "Vs4CaloPhi[nVs4Calo]/F");
  jetTreeIni_p->Branch("Vs4CaloEta", Vs4CaloEta_, "Vs4CaloEta[nVs4Calo]/F");
  jetTreeIni_p->Branch("Vs4CaloTrkMax", Vs4CaloTrkMax_, "Vs4CaloTrkMax[nVs4Calo]/F");
  jetTreeIni_p->Branch("Vs4CaloRawPt", Vs4CaloRawPt_, "Vs4CaloRawPt[nVs4Calo]/F");

  jetTreeIni_p->Branch("nVs5Calo", &nVs5Calo_, "nVs5Calo/I");
  jetTreeIni_p->Branch("Vs5CaloPt", Vs5CaloPt_, "Vs5CaloPt[nVs5Calo]/F");
  jetTreeIni_p->Branch("Vs5CaloPhi", Vs5CaloPhi_, "Vs5CaloPhi[nVs5Calo]/F");
  jetTreeIni_p->Branch("Vs5CaloEta", Vs5CaloEta_, "Vs5CaloEta[nVs5Calo]/F");
  jetTreeIni_p->Branch("Vs5CaloTrkMax", Vs5CaloTrkMax_, "Vs5CaloTrkMax[nVs5Calo]/F");
  jetTreeIni_p->Branch("Vs5CaloRawPt", Vs5CaloRawPt_, "Vs5CaloRawPt[nVs5Calo]/F");


  jetTreeIni_p->Branch("nVs2CaloFrag", &nVs2CaloFrag_, "nVs2CaloFrag/I");
  jetTreeIni_p->Branch("Vs2CaloFragPt", Vs2CaloFragPt_, "Vs2CaloFragPt[nVs2CaloFrag]/F");
  jetTreeIni_p->Branch("Vs2CaloFragPhi", Vs2CaloFragPhi_, "Vs2CaloFragPhi[nVs2CaloFrag]/F");
  jetTreeIni_p->Branch("Vs2CaloFragEta", Vs2CaloFragEta_, "Vs2CaloFragEta[nVs2CaloFrag]/F");
  jetTreeIni_p->Branch("Vs2CaloFragTrkMax", Vs2CaloFragTrkMax_, "Vs2CaloFragTrkMax[nVs2CaloFrag]/F");
  jetTreeIni_p->Branch("Vs2CaloFragRawPt", Vs2CaloFragRawPt_, "Vs2CaloFragRawPt[nVs2CaloFrag]/F");

  jetTreeIni_p->Branch("nVs3CaloFrag", &nVs3CaloFrag_, "nVs3CaloFrag/I");
  jetTreeIni_p->Branch("Vs3CaloFragPt", Vs3CaloFragPt_, "Vs3CaloFragPt[nVs3CaloFrag]/F");
  jetTreeIni_p->Branch("Vs3CaloFragPhi", Vs3CaloFragPhi_, "Vs3CaloFragPhi[nVs3CaloFrag]/F");
  jetTreeIni_p->Branch("Vs3CaloFragEta", Vs3CaloFragEta_, "Vs3CaloFragEta[nVs3CaloFrag]/F");
  jetTreeIni_p->Branch("Vs3CaloFragTrkMax", Vs3CaloFragTrkMax_, "Vs3CaloFragTrkMax[nVs3CaloFrag]/F");
  jetTreeIni_p->Branch("Vs3CaloFragRawPt", Vs3CaloFragRawPt_, "Vs3CaloFragRawPt[nVs3CaloFrag]/F");

  jetTreeIni_p->Branch("nVs4CaloFrag", &nVs4CaloFrag_, "nVs4CaloFrag/I");
  jetTreeIni_p->Branch("Vs4CaloFragPt", Vs4CaloFragPt_, "Vs4CaloFragPt[nVs4CaloFrag]/F");
  jetTreeIni_p->Branch("Vs4CaloFragPhi", Vs4CaloFragPhi_, "Vs4CaloFragPhi[nVs4CaloFrag]/F");
  jetTreeIni_p->Branch("Vs4CaloFragEta", Vs4CaloFragEta_, "Vs4CaloFragEta[nVs4CaloFrag]/F");
  jetTreeIni_p->Branch("Vs4CaloFragTrkMax", Vs4CaloFragTrkMax_, "Vs4CaloFragTrkMax[nVs4CaloFrag]/F");
  jetTreeIni_p->Branch("Vs4CaloFragRawPt", Vs4CaloFragRawPt_, "Vs4CaloFragRawPt[nVs4CaloFrag]/F");

  jetTreeIni_p->Branch("nVs5CaloFrag", &nVs5CaloFrag_, "nVs5CaloFrag/I");
  jetTreeIni_p->Branch("Vs5CaloFragPt", Vs5CaloFragPt_, "Vs5CaloFragPt[nVs5CaloFrag]/F");
  jetTreeIni_p->Branch("Vs5CaloFragPhi", Vs5CaloFragPhi_, "Vs5CaloFragPhi[nVs5CaloFrag]/F");
  jetTreeIni_p->Branch("Vs5CaloFragEta", Vs5CaloFragEta_, "Vs5CaloFragEta[nVs5CaloFrag]/F");
  jetTreeIni_p->Branch("Vs5CaloFragTrkMax", Vs5CaloFragTrkMax_, "Vs5CaloFragTrkMax[nVs5CaloFrag]/F");
  jetTreeIni_p->Branch("Vs5CaloFragRawPt", Vs5CaloFragRawPt_, "Vs5CaloFragRawPt[nVs5CaloFrag]/F");

  jetTreeIni_p->Branch("nVs2CaloRes", &nVs2CaloRes_, "nVs2CaloRes/I");
  jetTreeIni_p->Branch("Vs2CaloResPt", Vs2CaloResPt_, "Vs2CaloResPt[nVs2CaloRes]/F");
  jetTreeIni_p->Branch("Vs2CaloResPhi", Vs2CaloResPhi_, "Vs2CaloResPhi[nVs2CaloRes]/F");
  jetTreeIni_p->Branch("Vs2CaloResEta", Vs2CaloResEta_, "Vs2CaloResEta[nVs2CaloRes]/F");
  jetTreeIni_p->Branch("Vs2CaloResTrkMax", Vs2CaloResTrkMax_, "Vs2CaloResTrkMax[nVs2CaloRes]/F");
  jetTreeIni_p->Branch("Vs2CaloResRawPt", Vs2CaloResRawPt_, "Vs2CaloResRawPt[nVs2CaloRes]/F");

  jetTreeIni_p->Branch("nVs3CaloRes", &nVs3CaloRes_, "nVs3CaloRes/I");
  jetTreeIni_p->Branch("Vs3CaloResPt", Vs3CaloResPt_, "Vs3CaloResPt[nVs3CaloRes]/F");
  jetTreeIni_p->Branch("Vs3CaloResPhi", Vs3CaloResPhi_, "Vs3CaloResPhi[nVs3CaloRes]/F");
  jetTreeIni_p->Branch("Vs3CaloResEta", Vs3CaloResEta_, "Vs3CaloResEta[nVs3CaloRes]/F");
  jetTreeIni_p->Branch("Vs3CaloResTrkMax", Vs3CaloResTrkMax_, "Vs3CaloResTrkMax[nVs3CaloRes]/F");
  jetTreeIni_p->Branch("Vs3CaloResRawPt", Vs3CaloResRawPt_, "Vs3CaloResRawPt[nVs3CaloRes]/F");

  jetTreeIni_p->Branch("nVs4CaloRes", &nVs4CaloRes_, "nVs4CaloRes/I");
  jetTreeIni_p->Branch("Vs4CaloResPt", Vs4CaloResPt_, "Vs4CaloResPt[nVs4CaloRes]/F");
  jetTreeIni_p->Branch("Vs4CaloResPhi", Vs4CaloResPhi_, "Vs4CaloResPhi[nVs4CaloRes]/F");
  jetTreeIni_p->Branch("Vs4CaloResEta", Vs4CaloResEta_, "Vs4CaloResEta[nVs4CaloRes]/F");
  jetTreeIni_p->Branch("Vs4CaloResTrkMax", Vs4CaloResTrkMax_, "Vs4CaloResTrkMax[nVs4CaloRes]/F");
  jetTreeIni_p->Branch("Vs4CaloResRawPt", Vs4CaloResRawPt_, "Vs4CaloResRawPt[nVs4CaloRes]/F");

  jetTreeIni_p->Branch("nVs5CaloRes", &nVs5CaloRes_, "nVs5CaloRes/I");
  jetTreeIni_p->Branch("Vs5CaloResPt", Vs5CaloResPt_, "Vs5CaloResPt[nVs5CaloRes]/F");
  jetTreeIni_p->Branch("Vs5CaloResPhi", Vs5CaloResPhi_, "Vs5CaloResPhi[nVs5CaloRes]/F");
  jetTreeIni_p->Branch("Vs5CaloResEta", Vs5CaloResEta_, "Vs5CaloResEta[nVs5CaloRes]/F");
  jetTreeIni_p->Branch("Vs5CaloResTrkMax", Vs5CaloResTrkMax_, "Vs5CaloResTrkMax[nVs5CaloRes]/F");
  jetTreeIni_p->Branch("Vs5CaloResRawPt", Vs5CaloResRawPt_, "Vs5CaloResRawPt[nVs5CaloRes]/F");


  if(montecarlo){
    jetTreeIni_p->Branch("Vs2CaloRefPt", Vs2CaloRefPt_, "Vs2CaloRefPt[nVs2Calo]/F");
    jetTreeIni_p->Branch("Vs2CaloRefPhi", Vs2CaloRefPhi_, "Vs2CaloRefPhi[nVs2Calo]/F");
    jetTreeIni_p->Branch("Vs2CaloRefEta", Vs2CaloRefEta_, "Vs2CaloRefEta[nVs2Calo]/F");
    jetTreeIni_p->Branch("Vs2CaloRefPart", Vs2CaloRefPart_, "Vs2CaloRefPart[nVs2Calo]/I");

    jetTreeIni_p->Branch("Vs3CaloRefPt", Vs3CaloRefPt_, "Vs3CaloRefPt[nVs3Calo]/F");
    jetTreeIni_p->Branch("Vs3CaloRefPhi", Vs3CaloRefPhi_, "Vs3CaloRefPhi[nVs3Calo]/F");
    jetTreeIni_p->Branch("Vs3CaloRefEta", Vs3CaloRefEta_, "Vs3CaloRefEta[nVs3Calo]/F");
    jetTreeIni_p->Branch("Vs3CaloRefPart", Vs3CaloRefPart_, "Vs3CaloRefPart[nVs3Calo]/I");

    jetTreeIni_p->Branch("Vs4CaloRefPt", Vs4CaloRefPt_, "Vs4CaloRefPt[nVs4Calo]/F");
    jetTreeIni_p->Branch("Vs4CaloRefPhi", Vs4CaloRefPhi_, "Vs4CaloRefPhi[nVs4Calo]/F");
    jetTreeIni_p->Branch("Vs4CaloRefEta", Vs4CaloRefEta_, "Vs4CaloRefEta[nVs4Calo]/F");
    jetTreeIni_p->Branch("Vs4CaloRefPart", Vs4CaloRefPart_, "Vs4CaloRefPart[nVs4Calo]/I");

    jetTreeIni_p->Branch("Vs5CaloRefPt", Vs5CaloRefPt_, "Vs5CaloRefPt[nVs5Calo]/F");
    jetTreeIni_p->Branch("Vs5CaloRefPhi", Vs5CaloRefPhi_, "Vs5CaloRefPhi[nVs5Calo]/F");
    jetTreeIni_p->Branch("Vs5CaloRefEta", Vs5CaloRefEta_, "Vs5CaloRefEta[nVs5Calo]/F");
    jetTreeIni_p->Branch("Vs5CaloRefPart", Vs5CaloRefPart_, "Vs5CaloRefPart[nVs5Calo]/I");


    jetTreeIni_p->Branch("Vs2CaloFragRefPt", Vs2CaloFragRefPt_, "Vs2CaloFragRefPt[nVs2CaloFrag]/F");
    jetTreeIni_p->Branch("Vs2CaloFragRefPhi", Vs2CaloFragRefPhi_, "Vs2CaloFragRefPhi[nVs2CaloFrag]/F");
    jetTreeIni_p->Branch("Vs2CaloFragRefEta", Vs2CaloFragRefEta_, "Vs2CaloFragRefEta[nVs2CaloFrag]/F");
    jetTreeIni_p->Branch("Vs2CaloFragRefPart", Vs2CaloFragRefPart_, "Vs2CaloFragRefPart[nVs2CaloFrag]/I");

    jetTreeIni_p->Branch("Vs3CaloFragRefPt", Vs3CaloFragRefPt_, "Vs3CaloFragRefPt[nVs3CaloFrag]/F");
    jetTreeIni_p->Branch("Vs3CaloFragRefPhi", Vs3CaloFragRefPhi_, "Vs3CaloFragRefPhi[nVs3CaloFrag]/F");
    jetTreeIni_p->Branch("Vs3CaloFragRefEta", Vs3CaloFragRefEta_, "Vs3CaloFragRefEta[nVs3CaloFrag]/F");
    jetTreeIni_p->Branch("Vs3CaloFragRefPart", Vs3CaloFragRefPart_, "Vs3CaloFragRefPart[nVs3CaloFrag]/I");

    jetTreeIni_p->Branch("Vs4CaloFragRefPt", Vs4CaloFragRefPt_, "Vs4CaloFragRefPt[nVs4CaloFrag]/F");
    jetTreeIni_p->Branch("Vs4CaloFragRefPhi", Vs4CaloFragRefPhi_, "Vs4CaloFragRefPhi[nVs4CaloFrag]/F");
    jetTreeIni_p->Branch("Vs4CaloFragRefEta", Vs4CaloFragRefEta_, "Vs4CaloFragRefEta[nVs4CaloFrag]/F");
    jetTreeIni_p->Branch("Vs4CaloFragRefPart", Vs4CaloFragRefPart_, "Vs4CaloFragRefPart[nVs4CaloFrag]/I");

    jetTreeIni_p->Branch("Vs5CaloFragRefPt", Vs5CaloFragRefPt_, "Vs5CaloFragRefPt[nVs5CaloFrag]/F");
    jetTreeIni_p->Branch("Vs5CaloFragRefPhi", Vs5CaloFragRefPhi_, "Vs5CaloFragRefPhi[nVs5CaloFrag]/F");
    jetTreeIni_p->Branch("Vs5CaloFragRefEta", Vs5CaloFragRefEta_, "Vs5CaloFragRefEta[nVs5CaloFrag]/F");
    jetTreeIni_p->Branch("Vs5CaloFragRefPart", Vs5CaloFragRefPart_, "Vs5CaloFragRefPart[nVs5CaloFrag]/I");

    jetTreeIni_p->Branch("Vs2CaloResRefPt", Vs2CaloResRefPt_, "Vs2CaloResRefPt[nVs2CaloRes]/F");
    jetTreeIni_p->Branch("Vs2CaloResRefPhi", Vs2CaloResRefPhi_, "Vs2CaloResRefPhi[nVs2CaloRes]/F");
    jetTreeIni_p->Branch("Vs2CaloResRefEta", Vs2CaloResRefEta_, "Vs2CaloResRefEta[nVs2CaloRes]/F");
    jetTreeIni_p->Branch("Vs2CaloResRefPart", Vs2CaloResRefPart_, "Vs2CaloResRefPart[nVs2CaloRes]/I");

    jetTreeIni_p->Branch("Vs3CaloResRefPt", Vs3CaloResRefPt_, "Vs3CaloResRefPt[nVs3CaloRes]/F");
    jetTreeIni_p->Branch("Vs3CaloResRefPhi", Vs3CaloResRefPhi_, "Vs3CaloResRefPhi[nVs3CaloRes]/F");
    jetTreeIni_p->Branch("Vs3CaloResRefEta", Vs3CaloResRefEta_, "Vs3CaloResRefEta[nVs3CaloRes]/F");
    jetTreeIni_p->Branch("Vs3CaloResRefPart", Vs3CaloResRefPart_, "Vs3CaloResRefPart[nVs3CaloRes]/I");

    jetTreeIni_p->Branch("Vs4CaloResRefPt", Vs4CaloResRefPt_, "Vs4CaloResRefPt[nVs4CaloRes]/F");
    jetTreeIni_p->Branch("Vs4CaloResRefPhi", Vs4CaloResRefPhi_, "Vs4CaloResRefPhi[nVs4CaloRes]/F");
    jetTreeIni_p->Branch("Vs4CaloResRefEta", Vs4CaloResRefEta_, "Vs4CaloResRefEta[nVs4CaloRes]/F");
    jetTreeIni_p->Branch("Vs4CaloResRefPart", Vs4CaloResRefPart_, "Vs4CaloResRefPart[nVs4CaloRes]/I");

    jetTreeIni_p->Branch("Vs5CaloResRefPt", Vs5CaloResRefPt_, "Vs5CaloResRefPt[nVs5CaloRes]/F");
    jetTreeIni_p->Branch("Vs5CaloResRefPhi", Vs5CaloResRefPhi_, "Vs5CaloResRefPhi[nVs5CaloRes]/F");
    jetTreeIni_p->Branch("Vs5CaloResRefEta", Vs5CaloResRefEta_, "Vs5CaloResRefEta[nVs5CaloRes]/F");
    jetTreeIni_p->Branch("Vs5CaloResRefPart", Vs5CaloResRefPart_, "Vs5CaloResRefPart[nVs5CaloRes]/I");


    jetTreeIni_p->Branch("nT2", &nT2_, "nT2/I");
    jetTreeIni_p->Branch("T2Pt", T2Pt_, "T2Pt[nT2]/F");
    jetTreeIni_p->Branch("T2Phi", T2Phi_, "T2Phi[nT2]/F");
    jetTreeIni_p->Branch("T2Eta", T2Eta_, "T2Eta[nT2]/F");
    jetTreeIni_p->Branch("T2Part", T2Part_, "T2Part[nT2]/I");

    jetTreeIni_p->Branch("nT3", &nT3_, "nT3/I");
    jetTreeIni_p->Branch("T3Pt", T3Pt_, "T3Pt[nT3]/F");
    jetTreeIni_p->Branch("T3Phi", T3Phi_, "T3Phi[nT3]/F");
    jetTreeIni_p->Branch("T3Eta", T3Eta_, "T3Eta[nT3]/F");
    jetTreeIni_p->Branch("T3Part", T3Part_, "T3Part[nT3]/I");

    jetTreeIni_p->Branch("nT4", &nT4_, "nT4/I");
    jetTreeIni_p->Branch("T4Pt", T4Pt_, "T4Pt[nT4]/F");
    jetTreeIni_p->Branch("T4Phi", T4Phi_, "T4Phi[nT4]/F");
    jetTreeIni_p->Branch("T4Eta", T4Eta_, "T4Eta[nT4]/F");
    jetTreeIni_p->Branch("T4Part", T4Part_, "T4Part[nT4]/I");

    jetTreeIni_p->Branch("nT5", &nT5_, "nT5/I");
    jetTreeIni_p->Branch("T5Pt", T5Pt_, "T5Pt[nT5]/F");
    jetTreeIni_p->Branch("T5Phi", T5Phi_, "T5Phi[nT5]/F");
    jetTreeIni_p->Branch("T5Eta", T5Eta_, "T5Eta[nT5]/F");
    jetTreeIni_p->Branch("T5Part", T5Part_, "T5Part[nT5]/I");
  }

  jetTreeIni_p->Branch("nPu3PF", &nPu3PF_, "nPu3PF/I");
  jetTreeIni_p->Branch("Pu3PFPt", Pu3PFPt_, "Pu3PFPt[nPu3PF]/F");
  jetTreeIni_p->Branch("Pu3PFPhi", Pu3PFPhi_, "Pu3PFPhi[nPu3PF]/F");
  jetTreeIni_p->Branch("Pu3PFEta", Pu3PFEta_, "Pu3PFEta[nPu3PF]/F");
  jetTreeIni_p->Branch("Pu3PFTrkMax", Pu3PFTrkMax_, "Pu3PFTrkMax[nPu3PF]/F");
  jetTreeIni_p->Branch("Pu3PFRawPt", Pu3PFRawPt_, "Pu3PFRawPt[nPu3PF]/F");

  if(montecarlo){
    jetTreeIni_p->Branch("Pu3PFRefPt", Pu3PFRefPt_, "Pu3PFRefPt[nPu3PF]/F");
    jetTreeIni_p->Branch("Pu3PFRefPhi", Pu3PFRefPhi_, "Pu3PFRefPhi[nPu3PF]/F");
    jetTreeIni_p->Branch("Pu3PFRefEta", Pu3PFRefEta_, "Pu3PFRefEta[nPu3PF]/F");
    jetTreeIni_p->Branch("Pu3PFRefPart", Pu3PFRefPart_, "Pu3PFRefPart[nPu3PF]/I");
  }    

  jetTreeIni_p->Branch("nVs3PF", &nVs3PF_, "nVs3PF/I");
  jetTreeIni_p->Branch("Vs3PFPt", Vs3PFPt_, "Vs3PFPt[nVs3PF]/F");
  jetTreeIni_p->Branch("Vs3PFPhi", Vs3PFPhi_, "Vs3PFPhi[nVs3PF]/F");
  jetTreeIni_p->Branch("Vs3PFEta", Vs3PFEta_, "Vs3PFEta[nVs3PF]/F");
  jetTreeIni_p->Branch("Vs3PFTrkMax", Vs3PFTrkMax_, "Vs3PFTrkMax[nVs3PF]/F");
  jetTreeIni_p->Branch("Vs3PFRawPt", Vs3PFRawPt_, "Vs3PFRawPt[nVs3PF]/F");

  if(montecarlo){
    jetTreeIni_p->Branch("Vs3PFRefPt", Vs3PFRefPt_, "Vs3PFRefPt[nVs3PF]/F");
    jetTreeIni_p->Branch("Vs3PFRefPhi", Vs3PFRefPhi_, "Vs3PFRefPhi[nVs3PF]/F");
    jetTreeIni_p->Branch("Vs3PFRefEta", Vs3PFRefEta_, "Vs3PFRefEta[nVs3PF]/F");
    jetTreeIni_p->Branch("Vs3PFRefPart", Vs3PFRefPart_, "Vs3PFRefPart[nVs3PF]/I");
  }    

  if(justJt){
    jetTreeIni_p->Branch("TrkJtPt", TrkJtPt_, "TrkJtPt[5]/F");
    jetTreeIni_p->Branch("TrkJtPhi", TrkJtPhi_, "TrkJtPhi[5]/F");
    jetTreeIni_p->Branch("TrkJtEta", TrkJtEta_, "TrkJtEta[5]/F");

    jetTreeIni_p->Branch("nLeadJtConst", &nLeadJtConst_, "nLeadJtConst/I");
    jetTreeIni_p->Branch("TrkLeadJtConstPt", TrkLeadJtConstPt_, "TrkLeadJtConstPt[nLeadJtConst]/F");
    jetTreeIni_p->Branch("TrkLeadJtConstPhi", TrkLeadJtConstPhi_, "TrkLeadJtConstPhi[nLeadJtConst]/F");
    jetTreeIni_p->Branch("TrkLeadJtConstEta", TrkLeadJtConstEta_, "TrkLeadJtConstEta[nLeadJtConst]/F");
    jetTreeIni_p->Branch("TrkLeadJtConstCorr", TrkLeadJtConstCorr_, "TrkLeadJtConstCorr[nLeadJtConst]/F");

    jetTreeIni_p->Branch("nSubLeadJtConst", &nSubLeadJtConst_, "nSubLeadJtConst/I");
    jetTreeIni_p->Branch("TrkSubLeadJtConstPt", TrkSubLeadJtConstPt_, "TrkSubLeadJtConstPt[nSubLeadJtConst]/F");
    jetTreeIni_p->Branch("TrkSubLeadJtConstPhi", TrkSubLeadJtConstPhi_, "TrkSubLeadJtConstPhi[nSubLeadJtConst]/F");
    jetTreeIni_p->Branch("TrkSubLeadJtConstEta", TrkSubLeadJtConstEta_, "TrkSubLeadJtConstEta[nSubLeadJtConst]/F");
    jetTreeIni_p->Branch("TrkSubLeadJtConstCorr", TrkSubLeadJtConstCorr_, "TrkSubLeadJtConstCorr[nSubLeadJtConst]/F");

    jetTreeIni_p->Branch("nThirdJtConst", &nThirdJtConst_, "nThirdJtConst/I");
    jetTreeIni_p->Branch("TrkThirdJtConstPt", TrkThirdJtConstPt_, "TrkThirdJtConstPt[nThirdJtConst]/F");
    jetTreeIni_p->Branch("TrkThirdJtConstPhi", TrkThirdJtConstPhi_, "TrkThirdJtConstPhi[nThirdJtConst]/F");
    jetTreeIni_p->Branch("TrkThirdJtConstEta", TrkThirdJtConstEta_, "TrkThirdJtConstEta[nThirdJtConst]/F");
    jetTreeIni_p->Branch("TrkThirdJtConstCorr", TrkThirdJtConstCorr_, "TrkThirdJtConstCorr[nThirdJtConst]/F");

    jetTreeIni_p->Branch("nFourthJtConst", &nFourthJtConst_, "nFourthJtConst/I");
    jetTreeIni_p->Branch("TrkFourthJtConstPt", TrkFourthJtConstPt_, "TrkFourthJtConstPt[nFourthJtConst]/F");
    jetTreeIni_p->Branch("TrkFourthJtConstPhi", TrkFourthJtConstPhi_, "TrkFourthJtConstPhi[nFourthJtConst]/F");
    jetTreeIni_p->Branch("TrkFourthJtConstEta", TrkFourthJtConstEta_, "TrkFourthJtConstEta[nFourthJtConst]/F");
    jetTreeIni_p->Branch("TrkFourthJtConstCorr", TrkFourthJtConstCorr_, "TrkFourthJtConstCorr[nFourthJtConst]/F");

    jetTreeIni_p->Branch("nFifthJtConst", &nFifthJtConst_, "nFifthJtConst/I");
    jetTreeIni_p->Branch("TrkFifthJtConstPt", TrkFifthJtConstPt_, "TrkFifthJtConstPt[nFifthJtConst]/F");
    jetTreeIni_p->Branch("TrkFifthJtConstPhi", TrkFifthJtConstPhi_, "TrkFifthJtConstPhi[nFifthJtConst]/F");
    jetTreeIni_p->Branch("TrkFifthJtConstEta", TrkFifthJtConstEta_, "TrkFifthJtConstEta[nFifthJtConst]/F");
    jetTreeIni_p->Branch("TrkFifthJtConstCorr", TrkFifthJtConstCorr_, "TrkFifthJtConstCorr[nFifthJtConst]/F");
  }

  //Gen Tree Branches

  if(montecarlo && !justJt){
    genTreeIni_p->Branch("nGen", &nGen_, "nGen/I");
    
    genTreeIni_p->Branch("genPt", genPt_, "genPt[nGen]/F");
    genTreeIni_p->Branch("genPhi", genPhi_, "genPhi[nGen]/F");
    genTreeIni_p->Branch("genEta", genEta_, "genEta[nGen]/F");
    
  }    
}


void GetIniBranches(sampleType sType = kHIDATA, Bool_t justJt = false)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

  std::cout << "Get Branches" << std::endl;

  if(!justJt){
    //Track Tree Branches

    trackTreeIni_p->SetBranchAddress("nTrk", &nTrk_);
    trackTreeIni_p->SetBranchAddress("trkPt", trkPt_);
    trackTreeIni_p->SetBranchAddress("trkPhi", trkPhi_);
    trackTreeIni_p->SetBranchAddress("trkEta", trkEta_);

    pfCandTreeIni_p->SetBranchAddress("nPF", &nPF_);
    pfCandTreeIni_p->SetBranchAddress("pfPt", pfPt_);
    pfCandTreeIni_p->SetBranchAddress("pfPhi", pfPhi_);
    pfCandTreeIni_p->SetBranchAddress("pfEta", pfEta_);    
  }

  //Jet Tree Branches

  jetTreeIni_p->SetBranchAddress("runIni", &runIni_);
  jetTreeIni_p->SetBranchAddress("evtIni", &evtIni_);
  jetTreeIni_p->SetBranchAddress("lumiIni", &lumiIni_);

  if(hi)
    jetTreeIni_p->SetBranchAddress("hiBinIni", &hiBinIni_);

  if(montecarlo)
    jetTreeIni_p->SetBranchAddress("pthatIni", &pthatIni_);

  if(hi){
    jetTreeIni_p->SetBranchAddress("hiEvtPlaneIni", &hiEvtPlaneIni_);
    jetTreeIni_p->SetBranchAddress("psinIni", &psinIni_);
  }  

  jetTreeIni_p->SetBranchAddress("nPu3Calo", &nPu3Calo_);
  jetTreeIni_p->SetBranchAddress("Pu3CaloPt", Pu3CaloPt_);
  jetTreeIni_p->SetBranchAddress("Pu3CaloPhi", Pu3CaloPhi_);
  jetTreeIni_p->SetBranchAddress("Pu3CaloEta", Pu3CaloEta_);
  jetTreeIni_p->SetBranchAddress("Pu3CaloTrkMax", Pu3CaloTrkMax_);
  jetTreeIni_p->SetBranchAddress("Pu3CaloRawPt", Pu3CaloRawPt_);

  jetTreeIni_p->SetBranchAddress("nPu4Calo", &nPu4Calo_);
  jetTreeIni_p->SetBranchAddress("Pu4CaloPt", Pu4CaloPt_);
  jetTreeIni_p->SetBranchAddress("Pu4CaloPhi", Pu4CaloPhi_);
  jetTreeIni_p->SetBranchAddress("Pu4CaloEta", Pu4CaloEta_);
  jetTreeIni_p->SetBranchAddress("Pu4CaloTrkMax", Pu4CaloTrkMax_);
  jetTreeIni_p->SetBranchAddress("Pu4CaloRawPt", Pu4CaloRawPt_);

  jetTreeIni_p->SetBranchAddress("nPu5Calo", &nPu5Calo_);
  jetTreeIni_p->SetBranchAddress("Pu5CaloPt", Pu5CaloPt_);
  jetTreeIni_p->SetBranchAddress("Pu5CaloPhi", Pu5CaloPhi_);
  jetTreeIni_p->SetBranchAddress("Pu5CaloEta", Pu5CaloEta_);
  jetTreeIni_p->SetBranchAddress("Pu5CaloTrkMax", Pu5CaloTrkMax_);
  jetTreeIni_p->SetBranchAddress("Pu5CaloRawPt", Pu5CaloRawPt_);

  if(montecarlo){
    jetTreeIni_p->SetBranchAddress("Pu3CaloRefPt", Pu3CaloRefPt_);
    jetTreeIni_p->SetBranchAddress("Pu3CaloRefPhi", Pu3CaloRefPhi_);
    jetTreeIni_p->SetBranchAddress("Pu3CaloRefEta", Pu3CaloRefEta_);
    jetTreeIni_p->SetBranchAddress("Pu3CaloRefPart", Pu3CaloRefPart_);

    jetTreeIni_p->SetBranchAddress("Pu4CaloRefPt", Pu4CaloRefPt_);
    jetTreeIni_p->SetBranchAddress("Pu4CaloRefPhi", Pu4CaloRefPhi_);
    jetTreeIni_p->SetBranchAddress("Pu4CaloRefEta", Pu4CaloRefEta_);
    jetTreeIni_p->SetBranchAddress("Pu4CaloRefPart", Pu4CaloRefPart_);

    jetTreeIni_p->SetBranchAddress("Pu5CaloRefPt", Pu5CaloRefPt_);
    jetTreeIni_p->SetBranchAddress("Pu5CaloRefPhi", Pu5CaloRefPhi_);
    jetTreeIni_p->SetBranchAddress("Pu5CaloRefEta", Pu5CaloRefEta_);
    jetTreeIni_p->SetBranchAddress("Pu5CaloRefPart", Pu5CaloRefPart_);
  }

  jetTreeIni_p->SetBranchAddress("nVs2Calo", &nVs2Calo_);
  jetTreeIni_p->SetBranchAddress("Vs2CaloPt", Vs2CaloPt_);
  jetTreeIni_p->SetBranchAddress("Vs2CaloPhi", Vs2CaloPhi_);
  jetTreeIni_p->SetBranchAddress("Vs2CaloEta", Vs2CaloEta_);
  jetTreeIni_p->SetBranchAddress("Vs2CaloTrkMax", Vs2CaloTrkMax_);

  jetTreeIni_p->SetBranchAddress("nVs3Calo", &nVs3Calo_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloPt", Vs3CaloPt_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloPhi", Vs3CaloPhi_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloEta", Vs3CaloEta_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloTrkMax", Vs3CaloTrkMax_);

  jetTreeIni_p->SetBranchAddress("nVs4Calo", &nVs4Calo_);
  jetTreeIni_p->SetBranchAddress("Vs4CaloPt", Vs4CaloPt_);
  jetTreeIni_p->SetBranchAddress("Vs4CaloPhi", Vs4CaloPhi_);
  jetTreeIni_p->SetBranchAddress("Vs4CaloEta", Vs4CaloEta_);
  jetTreeIni_p->SetBranchAddress("Vs4CaloTrkMax", Vs4CaloTrkMax_);

  jetTreeIni_p->SetBranchAddress("nVs5Calo", &nVs5Calo_);
  jetTreeIni_p->SetBranchAddress("Vs5CaloPt", Vs5CaloPt_);
  jetTreeIni_p->SetBranchAddress("Vs5CaloPhi", Vs5CaloPhi_);
  jetTreeIni_p->SetBranchAddress("Vs5CaloEta", Vs5CaloEta_);
  jetTreeIni_p->SetBranchAddress("Vs5CaloTrkMax", Vs5CaloTrkMax_);


  jetTreeIni_p->SetBranchAddress("nVs2CaloFrag", &nVs2CaloFrag_);
  jetTreeIni_p->SetBranchAddress("Vs2CaloFragPt", Vs2CaloFragPt_);
  jetTreeIni_p->SetBranchAddress("Vs2CaloFragPhi", Vs2CaloFragPhi_);
  jetTreeIni_p->SetBranchAddress("Vs2CaloFragEta", Vs2CaloFragEta_);
  jetTreeIni_p->SetBranchAddress("Vs2CaloFragTrkMax", Vs2CaloFragTrkMax_);

  jetTreeIni_p->SetBranchAddress("nVs3CaloFrag", &nVs3CaloFrag_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloFragPt", Vs3CaloFragPt_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloFragPhi", Vs3CaloFragPhi_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloFragEta", Vs3CaloFragEta_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloFragTrkMax", Vs3CaloFragTrkMax_);

  jetTreeIni_p->SetBranchAddress("nVs4CaloFrag", &nVs4CaloFrag_);
  jetTreeIni_p->SetBranchAddress("Vs4CaloFragPt", Vs4CaloFragPt_);
  jetTreeIni_p->SetBranchAddress("Vs4CaloFragPhi", Vs4CaloFragPhi_);
  jetTreeIni_p->SetBranchAddress("Vs4CaloFragEta", Vs4CaloFragEta_);
  jetTreeIni_p->SetBranchAddress("Vs4CaloFragTrkMax", Vs4CaloFragTrkMax_);

  jetTreeIni_p->SetBranchAddress("nVs5CaloFrag", &nVs5CaloFrag_);
  jetTreeIni_p->SetBranchAddress("Vs5CaloFragPt", Vs5CaloFragPt_);
  jetTreeIni_p->SetBranchAddress("Vs5CaloFragPhi", Vs5CaloFragPhi_);
  jetTreeIni_p->SetBranchAddress("Vs5CaloFragEta", Vs5CaloFragEta_);
  jetTreeIni_p->SetBranchAddress("Vs5CaloFragTrkMax", Vs5CaloFragTrkMax_);

  jetTreeIni_p->SetBranchAddress("nVs2CaloRes", &nVs2CaloRes_);
  jetTreeIni_p->SetBranchAddress("Vs2CaloResPt", Vs2CaloResPt_);
  jetTreeIni_p->SetBranchAddress("Vs2CaloResPhi", Vs2CaloResPhi_);
  jetTreeIni_p->SetBranchAddress("Vs2CaloResEta", Vs2CaloResEta_);
  jetTreeIni_p->SetBranchAddress("Vs2CaloResTrkMax", Vs2CaloResTrkMax_);

  jetTreeIni_p->SetBranchAddress("nVs3CaloRes", &nVs3CaloRes_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloResPt", Vs3CaloResPt_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloResPhi", Vs3CaloResPhi_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloResEta", Vs3CaloResEta_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloResTrkMax", Vs3CaloResTrkMax_);

  jetTreeIni_p->SetBranchAddress("nVs4CaloRes", &nVs4CaloRes_);
  jetTreeIni_p->SetBranchAddress("Vs4CaloResPt", Vs4CaloResPt_);
  jetTreeIni_p->SetBranchAddress("Vs4CaloResPhi", Vs4CaloResPhi_);
  jetTreeIni_p->SetBranchAddress("Vs4CaloResEta", Vs4CaloResEta_);
  jetTreeIni_p->SetBranchAddress("Vs4CaloResTrkMax", Vs4CaloResTrkMax_);

  jetTreeIni_p->SetBranchAddress("nVs5CaloRes", &nVs5CaloRes_);
  jetTreeIni_p->SetBranchAddress("Vs5CaloResPt", Vs5CaloResPt_);
  jetTreeIni_p->SetBranchAddress("Vs5CaloResPhi", Vs5CaloResPhi_);
  jetTreeIni_p->SetBranchAddress("Vs5CaloResEta", Vs5CaloResEta_);
  jetTreeIni_p->SetBranchAddress("Vs5CaloResTrkMax", Vs5CaloResTrkMax_);

  if(montecarlo){
    jetTreeIni_p->SetBranchAddress("Vs2CaloRefPt", Vs2CaloRefPt_);
    jetTreeIni_p->SetBranchAddress("Vs2CaloRefPhi", Vs2CaloRefPhi_);
    jetTreeIni_p->SetBranchAddress("Vs2CaloRefEta", Vs2CaloRefEta_);
    jetTreeIni_p->SetBranchAddress("Vs2CaloRefPart", Vs2CaloRefPart_);

    jetTreeIni_p->SetBranchAddress("Vs3CaloRefPt", Vs3CaloRefPt_);
    jetTreeIni_p->SetBranchAddress("Vs3CaloRefPhi", Vs3CaloRefPhi_);
    jetTreeIni_p->SetBranchAddress("Vs3CaloRefEta", Vs3CaloRefEta_);
    jetTreeIni_p->SetBranchAddress("Vs3CaloRefPart", Vs3CaloRefPart_);

    jetTreeIni_p->SetBranchAddress("Vs4CaloRefPt", Vs4CaloRefPt_);
    jetTreeIni_p->SetBranchAddress("Vs4CaloRefPhi", Vs4CaloRefPhi_);
    jetTreeIni_p->SetBranchAddress("Vs4CaloRefEta", Vs4CaloRefEta_);
    jetTreeIni_p->SetBranchAddress("Vs4CaloRefPart", Vs4CaloRefPart_);

    jetTreeIni_p->SetBranchAddress("Vs5CaloRefPt", Vs5CaloRefPt_);
    jetTreeIni_p->SetBranchAddress("Vs5CaloRefPhi", Vs5CaloRefPhi_);
    jetTreeIni_p->SetBranchAddress("Vs5CaloRefEta", Vs5CaloRefEta_);
    jetTreeIni_p->SetBranchAddress("Vs5CaloRefPart", Vs5CaloRefPart_);

    jetTreeIni_p->SetBranchAddress("Vs2CaloFragRefPt", Vs2CaloFragRefPt_);
    jetTreeIni_p->SetBranchAddress("Vs2CaloFragRefPhi", Vs2CaloFragRefPhi_);
    jetTreeIni_p->SetBranchAddress("Vs2CaloFragRefEta", Vs2CaloFragRefEta_);
    jetTreeIni_p->SetBranchAddress("Vs2CaloFragRefPart", Vs2CaloFragRefPart_);

    jetTreeIni_p->SetBranchAddress("Vs3CaloFragRefPt", Vs3CaloFragRefPt_);
    jetTreeIni_p->SetBranchAddress("Vs3CaloFragRefPhi", Vs3CaloFragRefPhi_);
    jetTreeIni_p->SetBranchAddress("Vs3CaloFragRefEta", Vs3CaloFragRefEta_);
    jetTreeIni_p->SetBranchAddress("Vs3CaloFragRefPart", Vs3CaloFragRefPart_);

    jetTreeIni_p->SetBranchAddress("Vs4CaloFragRefPt", Vs4CaloFragRefPt_);
    jetTreeIni_p->SetBranchAddress("Vs4CaloFragRefPhi", Vs4CaloFragRefPhi_);
    jetTreeIni_p->SetBranchAddress("Vs4CaloFragRefEta", Vs4CaloFragRefEta_);
    jetTreeIni_p->SetBranchAddress("Vs4CaloFragRefPart", Vs4CaloFragRefPart_);

    jetTreeIni_p->SetBranchAddress("Vs5CaloFragRefPt", Vs5CaloFragRefPt_);
    jetTreeIni_p->SetBranchAddress("Vs5CaloFragRefPhi", Vs5CaloFragRefPhi_);
    jetTreeIni_p->SetBranchAddress("Vs5CaloFragRefEta", Vs5CaloFragRefEta_);
    jetTreeIni_p->SetBranchAddress("Vs5CaloFragRefPart", Vs5CaloFragRefPart_);

    jetTreeIni_p->SetBranchAddress("Vs2CaloResRefPt", Vs2CaloResRefPt_);
    jetTreeIni_p->SetBranchAddress("Vs2CaloResRefPhi", Vs2CaloResRefPhi_);
    jetTreeIni_p->SetBranchAddress("Vs2CaloResRefEta", Vs2CaloResRefEta_);
    jetTreeIni_p->SetBranchAddress("Vs2CaloResRefPart", Vs2CaloResRefPart_);

    jetTreeIni_p->SetBranchAddress("Vs3CaloResRefPt", Vs3CaloResRefPt_);
    jetTreeIni_p->SetBranchAddress("Vs3CaloResRefPhi", Vs3CaloResRefPhi_);
    jetTreeIni_p->SetBranchAddress("Vs3CaloResRefEta", Vs3CaloResRefEta_);
    jetTreeIni_p->SetBranchAddress("Vs3CaloResRefPart", Vs3CaloResRefPart_);

    jetTreeIni_p->SetBranchAddress("Vs4CaloResRefPt", Vs4CaloResRefPt_);
    jetTreeIni_p->SetBranchAddress("Vs4CaloResRefPhi", Vs4CaloResRefPhi_);
    jetTreeIni_p->SetBranchAddress("Vs4CaloResRefEta", Vs4CaloResRefEta_);
    jetTreeIni_p->SetBranchAddress("Vs4CaloResRefPart", Vs4CaloResRefPart_);

    jetTreeIni_p->SetBranchAddress("Vs5CaloResRefPt", Vs5CaloResRefPt_);
    jetTreeIni_p->SetBranchAddress("Vs5CaloResRefPhi", Vs5CaloResRefPhi_);
    jetTreeIni_p->SetBranchAddress("Vs5CaloResRefEta", Vs5CaloResRefEta_);
    jetTreeIni_p->SetBranchAddress("Vs5CaloResRefPart", Vs5CaloResRefPart_);


    jetTreeIni_p->SetBranchAddress("nT2", &nT2_);
    jetTreeIni_p->SetBranchAddress("T2Pt", T2Pt_);
    jetTreeIni_p->SetBranchAddress("T2Phi", T2Phi_);
    jetTreeIni_p->SetBranchAddress("T2Eta", T2Eta_);
    jetTreeIni_p->SetBranchAddress("T2Part", T2Part_);

    jetTreeIni_p->SetBranchAddress("nT3", &nT3_);
    jetTreeIni_p->SetBranchAddress("T3Pt", T3Pt_);
    jetTreeIni_p->SetBranchAddress("T3Phi", T3Phi_);
    jetTreeIni_p->SetBranchAddress("T3Eta", T3Eta_);
    jetTreeIni_p->SetBranchAddress("T3Part", T3Part_);

    jetTreeIni_p->SetBranchAddress("nT4", &nT4_);
    jetTreeIni_p->SetBranchAddress("T4Pt", T4Pt_);
    jetTreeIni_p->SetBranchAddress("T4Phi", T4Phi_);
    jetTreeIni_p->SetBranchAddress("T4Eta", T4Eta_);
    jetTreeIni_p->SetBranchAddress("T4Part", T4Part_);

    jetTreeIni_p->SetBranchAddress("nT5", &nT5_);
    jetTreeIni_p->SetBranchAddress("T5Pt", T5Pt_);
    jetTreeIni_p->SetBranchAddress("T5Phi", T5Phi_);
    jetTreeIni_p->SetBranchAddress("T5Eta", T5Eta_);
    jetTreeIni_p->SetBranchAddress("T5Part", T5Part_);
  }


  jetTreeIni_p->SetBranchAddress("nPu3PF", &nPu3PF_);
  jetTreeIni_p->SetBranchAddress("Pu3PFPt", Pu3PFPt_);
  jetTreeIni_p->SetBranchAddress("Pu3PFPhi", Pu3PFPhi_);
  jetTreeIni_p->SetBranchAddress("Pu3PFEta", Pu3PFEta_);
  jetTreeIni_p->SetBranchAddress("Pu3PFTrkMax", Pu3PFTrkMax_);
  jetTreeIni_p->SetBranchAddress("Pu3PFRawPt", Pu3PFRawPt_);

  if(montecarlo){
    jetTreeIni_p->SetBranchAddress("Pu3PFRefPt", Pu3PFRefPt_);
    jetTreeIni_p->SetBranchAddress("Pu3PFRefPhi", Pu3PFRefPhi_);
    jetTreeIni_p->SetBranchAddress("Pu3PFRefEta", Pu3PFRefEta_);
    jetTreeIni_p->SetBranchAddress("Pu3PFRefPart", Pu3PFRefPart_);
  }

  jetTreeIni_p->SetBranchAddress("nVs3PF", &nVs3PF_);
  jetTreeIni_p->SetBranchAddress("Vs3PFPt", Vs3PFPt_);
  jetTreeIni_p->SetBranchAddress("Vs3PFPhi", Vs3PFPhi_);
  jetTreeIni_p->SetBranchAddress("Vs3PFEta", Vs3PFEta_);
  jetTreeIni_p->SetBranchAddress("Vs3PFTrkMax", Vs3PFTrkMax_);
  jetTreeIni_p->SetBranchAddress("Vs3PFRawPt", Vs3PFRawPt_);

  if(montecarlo){
    jetTreeIni_p->SetBranchAddress("Vs3PFRefPt", Vs3PFRefPt_);
    jetTreeIni_p->SetBranchAddress("Vs3PFRefPhi", Vs3PFRefPhi_);
    jetTreeIni_p->SetBranchAddress("Vs3PFRefEta", Vs3PFRefEta_);
    jetTreeIni_p->SetBranchAddress("Vs3PFRefPart", Vs3PFRefPart_);
  }

  if(justJt){
    jetTreeIni_p->SetBranchAddress("TrkJtPt", TrkJtPt_);
    jetTreeIni_p->SetBranchAddress("TrkJtPhi", TrkJtPhi_);
    jetTreeIni_p->SetBranchAddress("TrkJtEta", TrkJtEta_);

    jetTreeIni_p->SetBranchAddress("nLeadJtConst", &nLeadJtConst_);
    jetTreeIni_p->SetBranchAddress("TrkLeadJtConstPt", TrkLeadJtConstPt_);
    jetTreeIni_p->SetBranchAddress("TrkLeadJtConstPhi", TrkLeadJtConstPhi_);
    jetTreeIni_p->SetBranchAddress("TrkLeadJtConstEta", TrkLeadJtConstEta_);
    jetTreeIni_p->SetBranchAddress("TrkLeadJtConstCorr", TrkLeadJtConstCorr_);

    jetTreeIni_p->SetBranchAddress("nSubLeadJtConst", &nSubLeadJtConst_);
    jetTreeIni_p->SetBranchAddress("TrkSubLeadJtConstPt", TrkSubLeadJtConstPt_);
    jetTreeIni_p->SetBranchAddress("TrkSubLeadJtConstPhi", TrkSubLeadJtConstPhi_);
    jetTreeIni_p->SetBranchAddress("TrkSubLeadJtConstEta", TrkSubLeadJtConstEta_);
    jetTreeIni_p->SetBranchAddress("TrkSubLeadJtConstCorr", TrkSubLeadJtConstCorr_);

    jetTreeIni_p->SetBranchAddress("nThirdJtConst", &nThirdJtConst_);
    jetTreeIni_p->SetBranchAddress("TrkThirdJtConstPt", TrkThirdJtConstPt_);
    jetTreeIni_p->SetBranchAddress("TrkThirdJtConstPhi", TrkThirdJtConstPhi_);
    jetTreeIni_p->SetBranchAddress("TrkThirdJtConstEta", TrkThirdJtConstEta_);
    jetTreeIni_p->SetBranchAddress("TrkThirdJtConstCorr", TrkThirdJtConstCorr_);

    jetTreeIni_p->SetBranchAddress("nFourthJtConst", &nFourthJtConst_);
    jetTreeIni_p->SetBranchAddress("TrkFourthJtConstPt", TrkFourthJtConstPt_);
    jetTreeIni_p->SetBranchAddress("TrkFourthJtConstPhi", TrkFourthJtConstPhi_);
    jetTreeIni_p->SetBranchAddress("TrkFourthJtConstEta", TrkFourthJtConstEta_);
    jetTreeIni_p->SetBranchAddress("TrkFourthJtConstCorr", TrkFourthJtConstCorr_);

    jetTreeIni_p->SetBranchAddress("nFifthJtConst", &nFifthJtConst_);
    jetTreeIni_p->SetBranchAddress("TrkFifthJtConstPt", TrkFifthJtConstPt_);
    jetTreeIni_p->SetBranchAddress("TrkFifthJtConstPhi", TrkFifthJtConstPhi_);
    jetTreeIni_p->SetBranchAddress("TrkFifthJtConstEta", TrkFifthJtConstEta_);
    jetTreeIni_p->SetBranchAddress("TrkFifthJtConstCorr", TrkFifthJtConstCorr_);
  }

  //Gen Tree Branches

  if(montecarlo && !justJt){
    genTreeIni_p->SetBranchAddress("nGen", &nGen_);
    genTreeIni_p->SetBranchAddress("genPt", genPt_);
    genTreeIni_p->SetBranchAddress("genPhi", genPhi_);
    genTreeIni_p->SetBranchAddress("genEta", genEta_);
  }
}


void InitDiJetIniSkim(sampleType sType = kHIDATA, Bool_t justJt = false)
{
  std::cout << "Init DiJet IniSkim" << std::endl;

  if(!justJt) trackTreeIni_p = new TTree("trackTreeIni", "trackTreeIni");
  if(!justJt) pfCandTreeIni_p = new TTree("pfCandTreeIni", "pfCandTreeIni");
  jetTreeIni_p = new TTree("jetTreeIni", "jetTreeIni");
  if(isMonteCarlo(sType) && !justJt) genTreeIni_p = new TTree("genTreeIni", "genTreeIni");

  SetIniBranches(sType, justJt);
}


void CleanupDiJetIniSkim()
{
  if(trackTreeIni_p != 0) delete trackTreeIni_p;
  if(pfCandTreeIni_p != 0) delete pfCandTreeIni_p;
  if(jetTreeIni_p != 0) delete jetTreeIni_p;
  if(genTreeIni_p != 0) delete genTreeIni_p;
}


void GetDiJetIniSkim(TFile* iniFile_p, sampleType sType = kHIDATA, Bool_t justJt = false)
{
  std::cout << "Get DiJet IniSkim" << std::endl;

  if(!justJt) trackTreeIni_p = (TTree*)iniFile_p->Get("trackTreeIni");
  if(!justJt) pfCandTreeIni_p = (TTree*)iniFile_p->Get("pfCandTreeIni");
  jetTreeIni_p = (TTree*)iniFile_p->Get("jetTreeIni");
  if(isMonteCarlo(sType) && !justJt) genTreeIni_p = (TTree*)iniFile_p->Get("genTreeIni");

  GetIniBranches(sType, justJt);
}


void InitTrkJts(Bool_t justJt = false)
{
  if(justJt){
    for(Int_t jtIter = 0; jtIter < 5; jtIter++){
      TrkJtPt_[jtIter] = -999;
      TrkJtPhi_[jtIter] = -999;
      TrkJtEta_[jtIter] = -999;
    }
  }
}

#endif




