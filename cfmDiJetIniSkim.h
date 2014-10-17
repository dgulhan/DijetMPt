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

enum AlgoType_PbPb{
  Pu3Calo, //0
  Pu4Calo, //1
  Pu5Calo, //2
  Vs2Calo, //3
  Vs3Calo, //4
  Vs4Calo, //5
  Vs5Calo, //6
  T,       //7
  PuPF,    //8
  VsPF     //9
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
Float_t trkCorr_[maxTracks];

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

Int_t nT3_;
Float_t T3Pt_[maxJets];
Float_t T3Phi_[maxJets];
Float_t T3Eta_[maxJets];
Int_t T3Part_[maxJets];

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
    trackTreeIni_p->Branch("trkCorr", trkCorr_, "trkCorr[nTrk]/F");

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

    jetTreeIni_p->Branch("nT3", &nT3_, "nT3/I");
    jetTreeIni_p->Branch("T3Pt", T3Pt_, "T3Pt[nT3]/F");
    jetTreeIni_p->Branch("T3Phi", T3Phi_, "T3Phi[nT3]/F");
    jetTreeIni_p->Branch("T3Eta", T3Eta_, "T3Eta[nT3]/F");
    jetTreeIni_p->Branch("T3Part", T3Part_, "T3Part[nT3]/I");
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
    trackTreeIni_p->SetBranchAddress("trkCorr", trkCorr_);

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

    jetTreeIni_p->SetBranchAddress("nT3", &nT3_);
    jetTreeIni_p->SetBranchAddress("T3Pt", T3Pt_);
    jetTreeIni_p->SetBranchAddress("T3Phi", T3Phi_);
    jetTreeIni_p->SetBranchAddress("T3Eta", T3Eta_);
    jetTreeIni_p->SetBranchAddress("T3Part", T3Part_);
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




