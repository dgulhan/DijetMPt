//=============================================                                 
// Author: Chris McGinn                                                         
//                                                                              
// DiJet Mean Tree Class (MC)                                                    
//                                                                              
//=============================================  
#ifndef cfmDiJetMeanTree_h
#define cfmDiJetMeanTree_h

#include "/net/hisrv0001/home/cfmcginn/DijetMPt/CMSSW_5_3_20/src/DijetAnalysisSkim/cfmDiJetAnaSkim.h"
#include <string>

TTree* trackTreeMean_p = 0;
TTree* genTreeMean_p = 0;

const Int_t nAjBins = 4;
const Float_t ajBins[nAjBins] = {0.11, 0.22, 0.33, 1.00};
const Float_t rBins[nRBins] = {0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00};
const Int_t nAjCentBins = 4;
const Float_t ajCentBins[nAjCentBins] = {20, 60, 100, 200};
const Int_t nRCentBins = 2;
const Float_t rCentBins[nRCentBins] = {60, 200};

//Track Tree Variables

Int_t rProjA_MULT_[2*nSumAlg][nPtBins][nAjCentBins][nAjBins];
Float_t rProjA_MEAN_[2*nSumAlg][nPtBins][nAjCentBins][nAjBins];
Float_t rProjA_SIG_[2*nSumAlg][nPtBins][nAjCentBins][nAjBins];
Float_t rProjA_WEIGHT_[2*nSumAlg][nPtBins][nAjCentBins][nAjBins];
Float_t rProjA_WEIGHT2_[2*nSumAlg][nPtBins][nAjCentBins][nAjBins];

Int_t rProjAR_MULT_[2*nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t rProjAR_MEAN_[2*nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t rProjAR_SIG_[2*nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t rProjAR_WEIGHT_[2*nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t rProjAR_WEIGHT2_[2*nSumAlg][nPtBins][nRCentBins][nRBins];

Int_t rProjARD_MULT_[2*nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t rProjARD_MEAN_[2*nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t rProjARD_SIG_[2*nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t rProjARD_WEIGHT_[2*nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t rProjARD_WEIGHT2_[2*nSumAlg][nPtBins][nRCentBins][nRBins];

Int_t rProjARU_MULT_[2*nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t rProjARU_MEAN_[2*nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t rProjARU_SIG_[2*nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t rProjARU_WEIGHT_[2*nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t rProjARU_WEIGHT2_[2*nSumAlg][nPtBins][nRCentBins][nRBins];

//Gen Tree Variables

Int_t gProjA_MULT_[nSumAlg][nPtBins][nAjCentBins][nAjBins];
Float_t gProjA_MEAN_[nSumAlg][nPtBins][nAjCentBins][nAjBins];
Float_t gProjA_SIG_[nSumAlg][nPtBins][nAjCentBins][nAjBins];
Float_t gProjA_WEIGHT_[nSumAlg][nPtBins][nAjCentBins][nAjBins];
Float_t gProjA_WEIGHT2_[nSumAlg][nPtBins][nAjCentBins][nAjBins];

Int_t gProjAR_MULT_[nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t gProjAR_MEAN_[nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t gProjAR_SIG_[nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t gProjAR_WEIGHT_[nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t gProjAR_WEIGHT2_[nSumAlg][nPtBins][nRCentBins][nRBins];

Int_t gProjARD_MULT_[nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t gProjARD_MEAN_[nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t gProjARD_SIG_[nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t gProjARD_WEIGHT_[nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t gProjARD_WEIGHT2_[nSumAlg][nPtBins][nRCentBins][nRBins];

Int_t gProjARU_MULT_[nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t gProjARU_MEAN_[nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t gProjARU_SIG_[nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t gProjARU_WEIGHT_[nSumAlg][nPtBins][nRCentBins][nRBins];
Float_t gProjARU_WEIGHT2_[nSumAlg][nPtBins][nRCentBins][nRBins];


//Cut var
const Float_t ajBounds[nAjBins] = {0.11, 0.22, 0.33, 1.00};

void SetMeanBranches(sampleType sType = kHIDATA)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  std::cout << "Branches Set" << std::endl;

  //Track Tree Branches
  //Tracks proj. onto Alg, ordered according to enum above, All, Cone, and NotCone

  trackTreeMean_p->Branch("rProjA_MULT", rProjA_MULT_, Form("rProjA_MULT[%d][%d][%d][%d]/I", 2*nSumAlg, nPtBins, nAjCentBins, nAjBins));
  trackTreeMean_p->Branch("rProjA_MEAN", rProjA_MEAN_, Form("rProjA_MEAN[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nAjCentBins, nAjBins));
  trackTreeMean_p->Branch("rProjA_SIG", rProjA_SIG_, Form("rProjA_SIG[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nAjCentBins, nAjBins));
  trackTreeMean_p->Branch("rProjA_WEIGHT", rProjA_WEIGHT_, Form("rProjA_WEIGHT[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nAjCentBins, nAjBins));
  trackTreeMean_p->Branch("rProjA_WEIGHT2", rProjA_WEIGHT2_, Form("rProjA_WEIGHT2[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nAjCentBins, nAjBins));


  trackTreeMean_p->Branch("rProjAR_MULT", rProjAR_MULT_, Form("rProjAR_MULT[%d][%d][%d][%d]/I", 2*nSumAlg, nPtBins, nRCentBins, nRBins));
  trackTreeMean_p->Branch("rProjAR_MEAN", rProjAR_MEAN_, Form("rProjAR_MEAN[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRCentBins, nRBins));
  trackTreeMean_p->Branch("rProjAR_SIG", rProjAR_SIG_, Form("rProjAR_SIG[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRCentBins, nRBins));
  trackTreeMean_p->Branch("rProjAR_WEIGHT", rProjAR_WEIGHT_, Form("rProjAR_WEIGHT[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRCentBins, nRBins));
  trackTreeMean_p->Branch("rProjAR_WEIGHT2", rProjAR_WEIGHT2_, Form("rProjAR_WEIGHT2[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRCentBins, nRBins));

  trackTreeMean_p->Branch("rProjARD_MULT", rProjARD_MULT_, Form("rProjARD_MULT[%d][%d][%d][%d]/I", 2*nSumAlg, nPtBins, nRCentBins, nRBins));
  trackTreeMean_p->Branch("rProjARD_MEAN", rProjARD_MEAN_, Form("rProjARD_MEAN[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRCentBins, nRBins));
  trackTreeMean_p->Branch("rProjARD_SIG", rProjARD_SIG_, Form("rProjARD_SIG[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRCentBins, nRBins));
  trackTreeMean_p->Branch("rProjARD_WEIGHT", rProjARD_WEIGHT_, Form("rProjARD_WEIGHT[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRCentBins, nRBins));
  trackTreeMean_p->Branch("rProjARD_WEIGHT2", rProjARD_WEIGHT2_, Form("rProjARD_WEIGHT2[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRCentBins, nRBins));

  trackTreeMean_p->Branch("rProjARU_MULT", rProjARU_MULT_, Form("rProjARU_MULT[%d][%d][%d][%d]/I", 2*nSumAlg, nPtBins, nRCentBins, nRBins));
  trackTreeMean_p->Branch("rProjARU_MEAN", rProjARU_MEAN_, Form("rProjARU_MEAN[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRCentBins, nRBins));
  trackTreeMean_p->Branch("rProjARU_SIG", rProjARU_SIG_, Form("rProjARU_SIG[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRCentBins, nRBins));
  trackTreeMean_p->Branch("rProjARU_WEIGHT", rProjARU_WEIGHT_, Form("rProjARU_WEIGHT[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRCentBins, nRBins));
  trackTreeMean_p->Branch("rProjARU_WEIGHT2", rProjARU_WEIGHT2_, Form("rProjARU_WEIGHT2[%d][%d][%d][%d]/F", 2*nSumAlg, nPtBins, nRCentBins, nRBins));



  //Gen Tree Branches
  //Gen. proj. onto jetAlg, array ordered according to enum
  if(montecarlo){
    genTreeMean_p->Branch("gProjA_MULT", gProjA_MULT_, Form("gProjA_MULT[%d][%d][%d][%d]I", nSumAlg, nPtBins, nAjCentBins, nAjBins));
    genTreeMean_p->Branch("gProjA_MEAN", gProjA_MEAN_, Form("gProjA_MEAN[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nAjCentBins, nAjBins));
    genTreeMean_p->Branch("gProjA_SIG", gProjA_SIG_, Form("gProjA_SIG[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nAjCentBins, nAjBins));
    genTreeMean_p->Branch("gProjA_WEIGHT", gProjA_WEIGHT_, Form("gProjA_WEIGHT[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nAjCentBins, nAjBins));
    genTreeMean_p->Branch("gProjA_WEIGHT2", gProjA_WEIGHT2_, Form("gProjA_WEIGHT2[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nAjCentBins, nAjBins));

    genTreeMean_p->Branch("gProjAR_MULT", gProjAR_MULT_, Form("gProjAR_MULT[%d][%d][%d][%d]I", nSumAlg, nPtBins, nRCentBins, nRBins));
    genTreeMean_p->Branch("gProjAR_MEAN", gProjAR_MEAN_, Form("gProjAR_MEAN[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nRCentBins, nRBins));
    genTreeMean_p->Branch("gProjAR_SIG", gProjAR_SIG_, Form("gProjAR_SIG[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nRCentBins, nRBins));
    genTreeMean_p->Branch("gProjAR_WEIGHT", gProjAR_WEIGHT_, Form("gProjAR_WEIGHT[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nRCentBins, nRBins));
    genTreeMean_p->Branch("gProjAR_WEIGHT2", gProjAR_WEIGHT2_, Form("gProjAR_WEIGHT2[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nRCentBins, nRBins));

    genTreeMean_p->Branch("gProjARD_MULT", gProjARD_MULT_, Form("gProjARD_MULT[%d][%d][%d][%d]I", nSumAlg, nPtBins, nRCentBins, nRBins));
    genTreeMean_p->Branch("gProjARD_MEAN", gProjARD_MEAN_, Form("gProjARD_MEAN[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nRCentBins, nRBins));
    genTreeMean_p->Branch("gProjARD_SIG", gProjARD_SIG_, Form("gProjARD_SIG[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nRCentBins, nRBins));
    genTreeMean_p->Branch("gProjARD_WEIGHT", gProjARD_WEIGHT_, Form("gProjARD_WEIGHT[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nRCentBins, nRBins));
    genTreeMean_p->Branch("gProjARD_WEIGHT2", gProjARD_WEIGHT2_, Form("gProjARD_WEIGHT2[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nRCentBins, nRBins));

    genTreeMean_p->Branch("gProjARU_MULT", gProjARU_MULT_, Form("gProjARU_MULT[%d][%d][%d][%d]I", nSumAlg, nPtBins, nRCentBins, nRBins));
    genTreeMean_p->Branch("gProjARU_MEAN", gProjARU_MEAN_, Form("gProjARU_MEAN[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nRCentBins, nRBins));
    genTreeMean_p->Branch("gProjARU_SIG", gProjARU_SIG_, Form("gProjARU_SIG[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nRCentBins, nRBins));
    genTreeMean_p->Branch("gProjARU_WEIGHT", gProjARU_WEIGHT_, Form("gProjARU_WEIGHT[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nRCentBins, nRBins));
    genTreeMean_p->Branch("gProjARU_WEIGHT2", gProjARU_WEIGHT2_, Form("gProjARU_WEIGHT2[%d][%d][%d][%d]/F", nSumAlg, nPtBins, nRCentBins, nRBins));
  }

  return;
}


void GetMeanBranches(sampleType sType = kHIDATA)
{
  Bool_t montecarlo = isMonteCarlo(sType);
  std::cout << "Get Branches" << std::endl;
  
  //Track Tree Branches

  trackTreeMean_p->SetBranchAddress("rProjA_MULT", rProjA_MULT_);
  trackTreeMean_p->SetBranchAddress("rProjA_MEAN", rProjA_MEAN_);
  trackTreeMean_p->SetBranchAddress("rProjA_SIG", rProjA_SIG_);
  trackTreeMean_p->SetBranchAddress("rProjA_WEIGHT", rProjA_WEIGHT_);
  trackTreeMean_p->SetBranchAddress("rProjA_WEIGHT2", rProjA_WEIGHT2_);

  trackTreeMean_p->SetBranchAddress("rProjAR_MULT", rProjAR_MULT_);
  trackTreeMean_p->SetBranchAddress("rProjAR_MEAN", rProjAR_MEAN_);
  trackTreeMean_p->SetBranchAddress("rProjAR_SIG", rProjAR_SIG_);
  trackTreeMean_p->SetBranchAddress("rProjAR_WEIGHT", rProjAR_WEIGHT_);
  trackTreeMean_p->SetBranchAddress("rProjAR_WEIGHT2", rProjAR_WEIGHT2_);

  trackTreeMean_p->SetBranchAddress("rProjARD_MULT", rProjARD_MULT_);
  trackTreeMean_p->SetBranchAddress("rProjARD_MEAN", rProjARD_MEAN_);
  trackTreeMean_p->SetBranchAddress("rProjARD_SIG", rProjARD_SIG_);
  trackTreeMean_p->SetBranchAddress("rProjARD_WEIGHT", rProjARD_WEIGHT_);
  trackTreeMean_p->SetBranchAddress("rProjARD_WEIGHT2", rProjARD_WEIGHT2_);

  trackTreeMean_p->SetBranchAddress("rProjARU_MULT", rProjARU_MULT_);
  trackTreeMean_p->SetBranchAddress("rProjARU_MEAN", rProjARU_MEAN_);
  trackTreeMean_p->SetBranchAddress("rProjARU_SIG", rProjARU_SIG_);
  trackTreeMean_p->SetBranchAddress("rProjARU_WEIGHT", rProjARU_WEIGHT_);
  trackTreeMean_p->SetBranchAddress("rProjARU_WEIGHT2", rProjARU_WEIGHT2_);
  
  //Gen Tree Variables
  if(montecarlo){
    genTreeMean_p->SetBranchAddress("gProjA_MULT", gProjA_MULT_);
    genTreeMean_p->SetBranchAddress("gProjA_MEAN", gProjA_MEAN_);
    genTreeMean_p->SetBranchAddress("gProjA_SIG", gProjA_SIG_);
    genTreeMean_p->SetBranchAddress("gProjA_WEIGHT", gProjA_WEIGHT_);
    genTreeMean_p->SetBranchAddress("gProjA_WEIGHT2", gProjA_WEIGHT2_);

    genTreeMean_p->SetBranchAddress("gProjAR_MULT", gProjAR_MULT_);
    genTreeMean_p->SetBranchAddress("gProjAR_MEAN", gProjAR_MEAN_);
    genTreeMean_p->SetBranchAddress("gProjAR_SIG", gProjAR_SIG_);
    genTreeMean_p->SetBranchAddress("gProjAR_WEIGHT", gProjAR_WEIGHT_);
    genTreeMean_p->SetBranchAddress("gProjAR_WEIGHT2", gProjAR_WEIGHT2_);

    genTreeMean_p->SetBranchAddress("gProjARD_MULT", gProjARD_MULT_);
    genTreeMean_p->SetBranchAddress("gProjARD_MEAN", gProjARD_MEAN_);
    genTreeMean_p->SetBranchAddress("gProjARD_SIG", gProjARD_SIG_);
    genTreeMean_p->SetBranchAddress("gProjARD_WEIGHT", gProjARD_WEIGHT_);
    genTreeMean_p->SetBranchAddress("gProjARD_WEIGHT2", gProjARD_WEIGHT2_);

    genTreeMean_p->SetBranchAddress("gProjARU_MULT", gProjARU_MULT_);
    genTreeMean_p->SetBranchAddress("gProjARU_MEAN", gProjARU_MEAN_);
    genTreeMean_p->SetBranchAddress("gProjARU_SIG", gProjARU_SIG_);
    genTreeMean_p->SetBranchAddress("gProjARU_WEIGHT", gProjARU_WEIGHT_);
    genTreeMean_p->SetBranchAddress("gProjARU_WEIGHT2", gProjARU_WEIGHT2_);
  }    

  return;
}


void InitDiJetMeanTree(sampleType sType = kHIDATA)
{
  std::cout << "Init DiJet MeanTree" << std::endl;

  trackTreeMean_p = new TTree("trackTreeMean", "trackTreeMean");
  if(isMonteCarlo(sType)) genTreeMean_p = new TTree("genTreeMean", "genTreeMean");

  SetMeanBranches(sType);

  return;
}


void CleanupDiJetMeanTree()
{
  if(trackTreeMean_p != 0) delete trackTreeMean_p;
  if(genTreeMean_p != 0) delete genTreeMean_p;

  return;
}


void GetDiJetMeanTree(TFile* meanFile_p, sampleType sType = kHIDATA){
  std::cout << "Get DiJet MeanTree" << std::endl;

  trackTreeMean_p = (TTree*)meanFile_p->Get("trackTreeMean");
  if(isMonteCarlo(sType)) genTreeMean_p = (TTree*)meanFile_p->Get("genTreeMean");

  GetMeanBranches(sType);

  return;
}


#endif
