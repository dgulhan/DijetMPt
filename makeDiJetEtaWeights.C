//=============================================
// Author: Chris McGinn                                                                                
//
// DiJet Eta Weight Maker
//
//=============================================

#include <string>
#include <iostream>

#include "TTree.h"
#include "TDatime.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCut.h"
#include "TMath.h"

const Float_t leadJtPt = 120.0;
const Float_t subleadJtPt = 50.0;
const Float_t jtEta = 1.6;
const Float_t jtDelPhi = 5*TMath::Pi()/6;

const std::string algString[7] = {"akPu3Calo", "akPu4Calo", "akPu5Calo", "akVs2Calo", "akVs3Calo", "akVs4Calo", "akVs5Calo"};

void makeDijetEtaWeights(TTree* inTree_p)
{
  Bool_t eventSet_[10];
  Float_t AlgJtPt_[10][4];
  Float_t AlgJtEta_[10][4];
  Float_t AlgJtDelPhi12_[10];

  inTree_p->SetBranchStatus("*", 0);

  inTree_p->SetBranchStatus("eventSet", 1);
  inTree_p->SetBranchStatus("AlgJtPt", 1);
  inTree_p->SetBranchStatus("AlgJtEta", 1);
  inTree_p->SetBranchStatus("AlgJtDelPhi12", 1);

  inTree_p->SetBranchAddress("eventSet", eventSet_);
  inTree_p->SetBranchAddress("AlgJtPt", AlgJtPt_);
  inTree_p->SetBranchAddress("AlgJtEta", AlgJtEta_);
  inTree_p->SetBranchAddress("AlgJtDelPhi12", AlgJtDelPhi12_);

  Int_t nEntries = inTree_p->GetEntries();
  Float_t eps = .01;
  Int_t count = 0;

  TFile* outFile_p;

  while(true){
    std::cout << count << std::endl;

    TH1F* leadEta_p = new TH1F("leadEta_h", "leadEta_h", 100, -jtEta, jtEta);
    TH1F* subleadEta_p = new TH1F("subleadEta_h", "subleadEta_h", 100, -jtEta, jtEta);
    TH1F* divHist_p = new TH1F(Form("divHist_h"), Form("divHist_h"), 100, -jtEta, jtEta);

    TH1F* leadWeight_p;
    TH1F* subleadWeight_p;
    outFile_p = new TFile("etaWeightFile.root", "UPDATE");
    if(count > 0) subleadWeight_p = (TH1F*)outFile_p->Get(Form("subleadEtaWeightHist_h"));
    if(count > 1) leadWeight_p = (TH1F*)outFile_p->Get(Form("leadEtaWeightHist_h"));

    std::cout << "A" << std::endl;

    for(Int_t iter = 0; iter < nEntries; iter++){
      inTree_p->GetEntry(iter);
      
      if(!eventSet_[4]) continue;
      if(AlgJtPt_[4][0] < leadJtPt || AlgJtPt_[4][1] < subleadJtPt) continue;
      if(TMath::Abs(AlgJtEta_[4][0]) > jtEta || TMath::Abs(AlgJtEta_[4][1]) > jtEta) continue;
      if(AlgJtDelPhi12_[4] < jtDelPhi) continue;
      
      Float_t fillWeight = 1.0;
      if(count > 0) fillWeight *= subleadWeight_p->GetBinContent(subleadWeight_p->FindBin(AlgJtEta_[4][1]));
      if(count > 1) fillWeight *= leadWeight_p->GetBinContent(leadWeight_p->FindBin(AlgJtEta_[4][0]));

      leadEta_p->Fill(AlgJtEta_[4][0], fillWeight);
      subleadEta_p->Fill(AlgJtEta_[4][1], fillWeight);
    }
						       
    leadEta_p->SetMinimum(0.0001);
    subleadEta_p->SetMinimum(0.0001);
    leadEta_p->Write("", TObject::kOverwrite);
    subleadEta_p->Write("", TObject::kOverwrite);
    
    std::cout << "B" << std::endl;
    
    if(count%2 == 0){
      divHist_p->Divide(leadEta_p, subleadEta_p);
      if(TMath::Abs(divHist_p->GetMaximum() - 1) < eps && TMath::Abs(divHist_p->GetMinimum() - 1) < eps) break;
      if(count/2 != 0) divHist_p->Multiply(subleadWeight_p);
      divHist_p->Write(Form("subleadEtaWeightHist_h"), TObject::kOverwrite);
    }
    else{
      divHist_p->Divide(subleadEta_p, leadEta_p);
      if(TMath::Abs(divHist_p->GetMaximum() - 1) < eps && TMath::Abs(divHist_p->GetMinimum() - 1) < eps) break;     
      if(count/2 != 0) divHist_p->Multiply(leadWeight_p);
      divHist_p->Write(Form("leadEtaWeightHist_h"), TObject::kOverwrite);
    }    
    std::cout << "D" << std::endl;


    outFile_p->Close();
    delete outFile_p;

    delete divHist_p;
    delete leadEta_p;
    delete subleadEta_p;

    std::cout << "E" << std::endl;

    count++;
  }    

  return;
}

void runMakeDijetEtaWeights(const std::string inName)
{
  TH1::SetDefaultSumw2();
  TFile* inFile_p = new TFile(inName.c_str(), "READ");
  TTree* getTree_p = (TTree*)inFile_p->Get("jetTreeAna");

  makeDijetEtaWeights(getTree_p);

  inFile_p->Close();
  delete inFile_p;

  return;
}
