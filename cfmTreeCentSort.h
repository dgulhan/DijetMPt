#include "TTree.h"
#include "TMath.h"
#include <vector>

std::vector<Int_t>* centEntryVect_p[200];
Int_t tempHiBin_ = 0;
Float_t tempVz_ = 0;
Int_t pcollisionEventSelection_ = 0;

void InitCentEntryVect()
{
  for(Int_t iter = 0; iter < 200; iter++){
    centEntryVect_p[iter] = new std::vector<Int_t>;
  }
  return;
}


void FillCentEntryVect(TTree* inTree_p)
{
  Int_t nentries = inTree_p->GetEntries();

  inTree_p->SetBranchStatus("*", 0);
  inTree_p->SetBranchStatus("hiBin", 1);
  inTree_p->SetBranchStatus("vz", 1);
  inTree_p->SetBranchStatus("pcollisionEventSelection", 1);

  inTree_p->SetBranchAddress("hiBin", &tempHiBin_);
  inTree_p->SetBranchAddress("vz", &tempVz_);
  inTree_p->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection_);

  for(Int_t iter = 0; iter < nentries; iter++){
    inTree_p->GetEntry(iter);
    if(!pcollisionEventSelection_) continue;
    if(TMath::Abs(tempVz_) > 15.) continue;

    centEntryVect_p[tempHiBin_]->push_back(iter);
  }
  return;
}


void RunFillCentEntryVect(TTree* inTree_p)
{
  InitCentEntryVect();
  FillCentEntryVect(inTree_p);
  return;
}


void CleanCentEntryVect()
{
  for(Int_t iter = 0; iter < 200; iter++){
    centEntryVect_p[iter]->clear();
    delete centEntryVect_p[iter];
  }
  return;
}


