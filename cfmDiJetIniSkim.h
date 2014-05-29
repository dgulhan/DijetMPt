//=============================================                                 
// Author: Chris McGinn                                                         
//                                                                              
// DiJet Initial Skim Class (MC)                                                    
//                                                                              
//=============================================  
#ifndef cfmDiJetSkim_h
#define cfmDiJetSkim_h

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"

enum sampleType {
  kHIDATA, //0
  kHIMC,   //1
  kPPDATA, //2
  kPPMC,   //3
  kPADATA, //4 
  kPAMC    //5
};


enum AlgoType_PbPb {
  PuCalo,
  VsCalo,
  T
};


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
TTree* jetTreeIni_p;
TTree* genTreeIni_p;

//Track Tree Variables

const int MAXTRKS = 20000; //From SetupTrackTree.h

Int_t nTrk_;
Float_t trkPt_[MAXTRKS];
Float_t trkPhi_[MAXTRKS];
Float_t trkEta_[MAXTRKS];

//Jet Tree Variables

const int MAXJETS = 504; //From SetupJetTree.h

Int_t run_;
Int_t evt_;
Int_t lumi_;
Int_t hiBin_;

Float_t pthat_;

Int_t nPu3Calo_;
Float_t Pu3CaloPt_[MAXJETS];
Float_t Pu3CaloPhi_[MAXJETS];
Float_t Pu3CaloEta_[MAXJETS];

Int_t nVs3Calo_;
Float_t Vs3CaloPt_[MAXJETS];
Float_t Vs3CaloPhi_[MAXJETS];
Float_t Vs3CaloEta_[MAXJETS];

Int_t nT3_;
Float_t T3Pt_[MAXJETS];
Float_t T3Phi_[MAXJETS];
Float_t T3Eta_[MAXJETS];

//Gen Tree Variables

const int MAXGEN = 50000; //From SetupGenParticleTree.h

Int_t nGen_;
Float_t genPt_[MAXGEN];
Float_t genPhi_[MAXGEN];
Float_t genEta_[MAXGEN];

void SetIniBranches(Bool_t montecarlo = false)
{
  //Track Tree Branches

  std::cout << "Branches Set" << std::endl;
  
  trackTreeIni_p->Branch("nTrk", &nTrk_, "nTrk/I");
  
  trackTreeIni_p->Branch("trkPt", &trkPt_, "trkPt[nTrk]/F");
  trackTreeIni_p->Branch("trkPhi", &trkPhi_, "trkPhi[nTrk]/F");
  trackTreeIni_p->Branch("trkEta", &trkEta_, "trkEta[nTrk]/F");
  
  //Jet Tree Branches

  jetTreeIni_p->Branch("run", &run_, "run/I");
  jetTreeIni_p->Branch("evt", &evt_, "evt/I");
  jetTreeIni_p->Branch("lumi", &lumi_, "lumi/I");
  jetTreeIni_p->Branch("hiBin", &hiBin_, "hiBin/I");

   
  if(montecarlo)
    jetTreeIni_p->Branch("pthat", &pthat_, "pthat/F");

  
  jetTreeIni_p->Branch("nPu3Calo", &nPu3Calo_, "nPu3Calo/I");
  jetTreeIni_p->Branch("Pu3CaloPt", &Pu3CaloPt_, "Pu3CaloPt[nPu3Calo]/F");
  jetTreeIni_p->Branch("Pu3CaloPhi", &Pu3CaloPhi_, "Pu3CaloPhi[nPu3Calo]/F");
  jetTreeIni_p->Branch("Pu3CaloEta", &Pu3CaloEta_, "Pu3CaloEta[nPu3Calo]/F");

  jetTreeIni_p->Branch("nVs3Calo", &nVs3Calo_, "nVs3Calo/I");
  jetTreeIni_p->Branch("Vs3CaloPt", &Vs3CaloPt_, "Vs3CaloPt[nVs3Calo]/F");
  jetTreeIni_p->Branch("Vs3CaloPhi", &Vs3CaloPhi_, "Vs3CaloPhi[nVs3Calo]/F");
  jetTreeIni_p->Branch("Vs3CaloEta", &Vs3CaloEta_, "Vs3CaloEta[nVs3Calo]/F");

  jetTreeIni_p->Branch("nT3", &nT3_, "nT3/I");
  jetTreeIni_p->Branch("T3Pt", &T3Pt_, "T3Pt[nT3]/F");
  jetTreeIni_p->Branch("T3Phi", &T3Phi_, "T3Phi[nT3]/F");
  jetTreeIni_p->Branch("T3Eta", &T3Eta_, "T3Eta[nT3]/F");
  

  //Gen Tree Branches

  if(montecarlo){
    genTreeIni_p->Branch("nGen", &nGen_, "nGen/I");
    
    genTreeIni_p->Branch("genPt", &genPt_, "genPt[nGen]/F");
    genTreeIni_p->Branch("genPhi", &genPhi_, "genPhi[nGen]/F");
    genTreeIni_p->Branch("genEta", &genEta_, "genEta[nGen]/F");
    
  }    
}


void getIniBranches(Bool_t montecarlo = false)
{
  //Track Tree Branches

  std::cout << "Get Branches" << std::endl;

  trackTreeIni_p->SetBranchAddress("nTrk", &nTrk_);
  trackTreeIni_p->SetBranchAddress("trkPt", &trkPt_);
  trackTreeIni_p->SetBranchAddress("trkPhi", &trkPhi_);
  trackTreeIni_p->SetBranchAddress("trkEta", &trkEta_);

  //Jet Tree Branches

  jetTreeIni_p->SetBranchAddress("run", &run_);
  jetTreeIni_p->SetBranchAddress("evt", &evt_);
  jetTreeIni_p->SetBranchAddress("lumi", &lumi_);
  jetTreeIni_p->SetBranchAddress("hiBin", &hiBin_);

  if(montecarlo)
    jetTreeIni_p->SetBranchAddress("pthat", &pthat_);

  jetTreeIni_p->SetBranchAddress("nPu3Calo", &nPu3Calo_);
  jetTreeIni_p->SetBranchAddress("Pu3CaloPt", &Pu3CaloPt_);
  jetTreeIni_p->SetBranchAddress("Pu3CaloPhi", &Pu3CaloPhi_);
  jetTreeIni_p->SetBranchAddress("Pu3CaloEta", &Pu3CaloEta_);

  jetTreeIni_p->SetBranchAddress("nVs3Calo", &nVs3Calo_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloPt", &Vs3CaloPt_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloPhi", &Vs3CaloPhi_);
  jetTreeIni_p->SetBranchAddress("Vs3CaloEta", &Vs3CaloEta_);

  jetTreeIni_p->SetBranchAddress("nT3", &nT3_);
  jetTreeIni_p->SetBranchAddress("T3Pt", &T3Pt_);
  jetTreeIni_p->SetBranchAddress("T3Phi", &T3Phi_);
  jetTreeIni_p->SetBranchAddress("T3Eta", &T3Eta_);

  //Gen Tree Branches

  if(montecarlo){
    genTreeIni_p->Branch("nGen", &nGen_);
    genTreeIni_p->Branch("genPt", &genPt_);
    genTreeIni_p->Branch("genPhi", &genPhi_);
    genTreeIni_p->Branch("genEta", &genEta_);
  }
}


void InitDiJetIniSkim(Bool_t montecarlo = false)
{
  std::cout << "Init DiJet IniSkim" << std::endl;

  trackTreeIni_p = new TTree("trackTreeIni", "trackTreeIni");
  jetTreeIni_p = new TTree("jetTreeIni", "jetTreeIni");

  if(montecarlo)
    genTreeIni_p = new TTree("genTreeIni", "genTreeIni");

  SetIniBranches(montecarlo);
}


void GetDiJetIniSkim(TFile* iniFile, Bool_t montecarlo = false)
{
  std::cout << "Get DiJet IniSkim" << std::endl;

  trackTreeIni_p = (TTree*)inFile_p->Get("trackTreeIni");
  jetTreeIni_p = (TTree*)inFile_p->Get("jetTreeIni");

  if(montecarlo)
    genTreeIni_p = (TTree*)inFile_p->Get("genTreeIni");

  GetIniBranches(montecarlo);
}


#endif
