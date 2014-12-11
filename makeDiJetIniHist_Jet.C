#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"

#include "TCanvas.h"
#include "TLine.h"

#include "commonSetup.h"

#include <iostream>
#include <fstream>
#include <string>

#include "commonUtility.h"
#include "TMath.h"
#include "TLatex.h"

#include <vector>
#include <algorithm>
#include "cfmVectFunc.h"

TChain* getChain_p[3] = {0, 0, 0};

const Float_t fracR3 = 0.940546;
const Float_t jtPtCut = 40.0;

const Float_t ajLowBound = (120.0-50.0)/(120.0+50.0);

Int_t pcollisionEventSelection_;
Float_t vz_;
Float_t pthat_;
Int_t hiBin_;

Int_t nref_;
Float_t jtpt_[maxJets];
Float_t jteta_[maxJets];
Float_t rawpt_[maxJets];
Float_t refpt_[maxJets];
Float_t refeta_[maxJets];

Int_t ngen_;
Float_t genpt_[maxJets];
Float_t geneta_[maxJets];
Int_t genmatchindex_[maxJets];

void BookChain(Bool_t isPbPb, Bool_t isMonteCarlo)
{
  getChain_p[0]->SetBranchStatus("*", 0);

  if(isPbPb) getChain_p[0]->SetBranchStatus("pcollisionEventSelection", 1);

  getChain_p[0]->SetBranchStatus("vz", 1);
  if(isMonteCarlo) getChain_p[0]->SetBranchStatus("pthat", 1);
  if(isPbPb) getChain_p[0]->SetBranchStatus("hiBin", 1);

  getChain_p[0]->SetBranchStatus("nref", 1);
  getChain_p[0]->SetBranchStatus("jtpt", 1);
  getChain_p[0]->SetBranchStatus("jteta", 1);
  getChain_p[0]->SetBranchStatus("rawpt", 1);

  if(isMonteCarlo){
    getChain_p[0]->SetBranchStatus("refpt", 1);
    getChain_p[0]->SetBranchStatus("refeta", 1);

    getChain_p[0]->SetBranchStatus("ngen", 1);
    getChain_p[0]->SetBranchStatus("genpt", 1);
    getChain_p[0]->SetBranchStatus("geneta", 1);
    getChain_p[0]->SetBranchStatus("genmatchindex", 1);
  }

  if(isPbPb) getChain_p[0]->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection_);

  getChain_p[0]->SetBranchAddress("vz", &vz_);
  if(isMonteCarlo) getChain_p[0]->SetBranchAddress("pthat", &pthat_);
  if(isPbPb) getChain_p[0]->SetBranchAddress("hiBin", &hiBin_);

  getChain_p[0]->SetBranchAddress("nref", &nref_);
  getChain_p[0]->SetBranchAddress("jtpt", &jtpt_);
  getChain_p[0]->SetBranchAddress("jteta", &jteta_);
  getChain_p[0]->SetBranchAddress("rawpt", &rawpt_);

  if(isMonteCarlo){
    getChain_p[0]->SetBranchAddress("refpt", &refpt_);
    getChain_p[0]->SetBranchAddress("refeta", &refeta_);

    getChain_p[0]->SetBranchAddress("ngen", &ngen_);
    getChain_p[0]->SetBranchAddress("genpt", &genpt_);
    getChain_p[0]->SetBranchAddress("geneta", &geneta_);
    getChain_p[0]->SetBranchAddress("genmatchindex", &genmatchindex_);
  }

  return;
}


void GetChain(std::vector<std::string> inList, const std::string alg, Bool_t isPbPb, Bool_t isMonteCarlo)
{
  getChain_p[0] = new TChain(Form("%sJetAnalyzer/t", alg.c_str()));
  getChain_p[1] = new TChain("skimanalysis/HltTree");
  getChain_p[2] = new TChain("hiEvtAnalyzer/HiTree");

  for(Int_t iter = 0; iter < (Int_t)(inList.size()); iter++){
    getChain_p[0]->Add(inList[iter].c_str());
    getChain_p[1]->Add(inList[iter].c_str());
    getChain_p[2]->Add(inList[iter].c_str());
  }

  getChain_p[0]->AddFriend(getChain_p[2]);
  getChain_p[0]->AddFriend(getChain_p[1]);

  BookChain(isPbPb, isMonteCarlo);

  return;
}


void CleanChain()
{
  if(getChain_p[0] != 0){
    delete getChain_p[0];
    getChain_p[0] = 0;
  }

  if(getChain_p[1] != 0){
    delete getChain_p[1];
    getChain_p[1] = 0;
  }

  if(getChain_p[2] != 0){
    delete getChain_p[2];
    getChain_p[2] = 0;
  }
  return;
}


Float_t getHatWeight(Float_t inHat, Bool_t isPbPb, Bool_t isMonteCarlo)
{
  if(!isMonteCarlo) return 1;

  if(isPbPb){
    return 1;
  }
  else{
    return 1;
  }

  std::cout << inHat << std::endl;
  std::cout << "No weight assigned; check for error." << std::endl;

  return 0;
}


void getLogBins(const Float_t lower, const Float_t higher, const Int_t nBins, Float_t bins[])
{
  Float_t logBins[nBins+1];
  bins[0] = lower;
  bins[nBins] = higher;

  logBins[0] = TMath::Log10(lower);
  logBins[nBins] = TMath::Log10(higher);

  Float_t interval = (logBins[nBins] - logBins[0])/nBins;

  for(Int_t iter = 1; iter < nBins; iter++){
    logBins[iter] = logBins[0] + iter*interval;
    bins[iter] = TMath::Power(10, logBins[iter]);
  }

  return;
}


void getEtaBins(const Float_t lower, const Float_t higher, const Int_t nBins, Float_t bins[])
{
  Float_t interval = (higher - lower)/nBins;

  for(Int_t iter = 0; iter < nBins+1; iter++){
    bins[iter] = lower + interval*iter;
  }
  bins[0] = -1.99999999999;
  bins[nBins] = 1.99999999999;

  return;
}


Int_t posSearch(Float_t val, Int_t arrSize, Float_t arr[])
{
  Int_t pos = arrSize/2;
  Int_t low = 0;
  Int_t high = arrSize-1;

  Int_t iter = 0;

  if(val > arr[high]) return high - 1;

  while(val < arr[pos] || val >= arr[pos+1]){
    iter++;

    if(val < arr[pos]) high = pos;
    else low = pos;
    pos = (high+low)/2;
  }

  return pos;
}


void CleanTH1(TH1F* cleanHist_p)
{
  delete cleanHist_p;
  cleanHist_p = 0;
}


void makeDiJetIniHist(std::vector<std::string> inList, const std::string outName, const std::string alg = "akVs3Calo", Bool_t isPbPb = false, Bool_t isMonteCarlo = false)
{
  TH1::SetDefaultSumw2();

  std::cout << "Init: " << alg << std::endl;

  GetChain(inList, alg, isPbPb, isMonteCarlo);

  const Int_t nPtBins = 25;
  const Float_t ptLower = jtPtCut;
  const Float_t ptHigher = 700.;

  const Int_t nEtaBins = 25;
  const Float_t etaLower = -2.0;
  const Float_t etaHigher = 2.0;

  //Edit Here

  Float_t ptBins[nPtBins+1];
  Float_t etaBins[nEtaBins+1];

  getLogBins(ptLower, ptHigher, nPtBins, ptBins);
  getEtaBins(etaLower, etaHigher, nEtaBins, etaBins);

  std::vector<Float_t>* jtPt_p = new std::vector<Float_t>;

  TH1F* jetPt_p = new TH1F(Form("%s_jetPt_h", alg.c_str()), Form("%s_jetPt_h", alg.c_str()), nPtBins, ptBins);
  TH1F* jetEta_p = new TH1F(Form("%s_jetEta_h", alg.c_str()), Form("%s_jetEta_h", alg.c_str()), nEtaBins, etaBins);;

  Int_t nentries = getChain_p[0]->GetEntries();

  for(Int_t jEntry = 0; jEntry < nentries; jEntry++){
    getChain_p[0]->GetEntry(jEntry);

    if(jEntry%200000 == 0) std::cout << jEntry << std::endl;

    if(!pcollisionEventSelection_ && isPbPb) continue;

    if(TMath::Abs(vz_) > 15) continue;

    //    Float_t hatWeight = getHatWeight(pthat_, isPbPb, isMonteCarlo);

    for(Int_t jtIter = 0; jtIter < nref_; jtIter++){
      if(TMath::Abs(jteta_[jtIter]) >= 2.0) continue;

      if(jtpt_[jtIter] < jtPtCut) break;

      //      jetPt_p->Fill(jtpt_[jtIter], hatWeight);
      //      jetEta_p->Fill(jteta_[jtIter], hatWeight);

      jtPt_p->push_back(jtpt_[jtIter]);

      break;
    }
  }

  std::sort(jtPt_p->begin(), jtPt_p->end());
  Int_t cutPos = -1;
  Int_t fracPos = -1;

  Int_t nJts = jtPt_p->size();

  for(Int_t jtIter = 0; jtIter < nJts; jtIter++){
    if(jtPt_p->at(jtIter) > 120.00){
      cutPos = jtIter;
      break;
   }
  }

  for(Int_t jtIter = 0; jtIter < nJts; jtIter++){
    if(((Float_t)jtIter)/((Float_t)nJts) > fracR3){
      fracPos = jtIter;
      break;
   }
  }

  std::cout << "CutPos, Percentile cut, " << alg << ": " << cutPos << ", " << ((Float_t)cutPos)/((Float_t)nJts) << std::endl;

  Float_t subPt = jtPt_p->at(fracPos)*(1-ajLowBound)/(1+ajLowBound);

  std::cout << "FracPos, ptCut1, ptCut2, frac, " << alg << ": " << fracPos << ", " << jtPt_p->at(fracPos) << ", " << subPt << ", " << ((Float_t)fracPos)/((Float_t)nJts) << std::endl;

  TFile* out = new TFile(outName.c_str(), "UPDATE");
  std::cout << outName << std::endl;

  jetPt_p->Write("", TObject::kOverwrite);
  jetEta_p->Write("", TObject::kOverwrite);

  out->Close();
  delete out;

  CleanTH1(jetPt_p);
  CleanTH1(jetEta_p);

  CleanChain();

  jtPt_p->clear();
  delete jtPt_p;

  return;
}


int runMakeDiJetIniHist(std::string fList = "", const char* outFileName = "raw_rawOverGen", Bool_t isPbPb = false, Bool_t isMonteCarlo = false)
{
  TH1::SetDefaultSumw2();

  std::string buffer;
  std::vector<std::string> listOfFiles;
  int nLines = 0;
  ifstream inFile(fList.data());

  std::cout << fList << std::endl;
  std::cout << inFile.is_open() << std::endl;

  if(!inFile.is_open()){
    std::cout << "Error opening file. Exiting." << std::endl;
    return 1;
  }
  else{
    while(true){
      inFile >> buffer;
      if(inFile.eof()) break;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }

  const std::string PFCaloString[2] = {"PF", "Calo"};

  for(Int_t numIter = 1; numIter < 5; numIter++){
    for(Int_t pfCaloIter = 1; pfCaloIter < 2; pfCaloIter++){      
      makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), Form("akVs%d%s", numIter+1, PFCaloString[pfCaloIter].c_str()), isPbPb, isMonteCarlo);
    }    
  }

  return(0);
}


int main(int argc, char* argv[])
{
  if(argc != 5){
    std::cout << "Usage: runMakeDiJetIniHist <inputList> <outFileName> <isPbPb> <isMonteCarlo>" << std::endl;
    return 1;
  }

  int rStatus = -1;

  rStatus = runMakeDiJetIniHist(argv[1], argv[2], (Bool_t)(atoi(argv[3])), (Bool_t)(atoi(argv[4])));

  return rStatus;
}
