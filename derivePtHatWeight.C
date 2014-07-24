#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"

void derivePtHatWeights(std::string fList = "")
{
  std::string buffer;
  std::vector<std::string> listOfFiles;
  int nLines = 0;
  ifstream inFile(fList.data());

  std::cout << fList << std::endl;
  std::cout << inFile.is_open() << std::endl;

  if(!inFile.is_open()){
    std::cout << "Error opening file. Exiting." <<std::endl;
    return;
  }
  else{
    while(!inFile.eof()){
      inFile >> buffer;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }

  std::cout << "FileList Loaded" << std::endl;

  TChain* ptHatChain_p = new TChain("akVs3CaloJetAnalyzer/t");

  for(Int_t iter = 0; iter < (Int_t)(listOfFiles.size()); iter++){
    ptHatChain_p->Add(listOfFiles[iter].c_str());
  }

  Float_t ptHat_ = 0;

  ptHatChain_p->SetBranchStatus("*", 0);
  ptHatChain_p->SetBranchStatus("pthat", 1);
  ptHatChain_p->SetBranchAddress("pthat", &ptHat_);

  Int_t nEntries = ptHatChain_p->GetEntries();
  std::cout << nEntries << std::endl;

  Int_t hatEntries[5] = {0, 0, 0, 0, 0};
  Int_t ptHatCuts[6] = {30, 50, 80, 120, 170, 100000};
  Float_t crossSections[6] = {.01075, .001025, .00009865, .00001129, .000001448, 0.000000000};
  Float_t hatWeight[5] = {0, 0, 0, 0, 0};

  for(Int_t evtIter = 0; evtIter < nEntries; evtIter++){
    ptHatChain_p->GetEntry(evtIter);

    for(Int_t hatIter = 0; hatIter < 5; hatIter++){
      if(ptHat_ > ptHatCuts[hatIter] && ptHat_ < ptHatCuts[hatIter+1]){
	hatEntries[hatIter]++;
	break;
      }
    }
  }

  for(Int_t hatIter = 0; hatIter < 5; hatIter++){
    std::cout << hatIter << std::endl;
    std::cout << "  hatEntries: " << hatEntries[hatIter] << std::endl;
    hatWeight[hatIter] = (crossSections[hatIter] - crossSections[hatIter+1])/hatEntries[hatIter];
    std::cout << "  hatWeight: " << hatWeight[hatIter] << std::endl;
  }

  return;
}
