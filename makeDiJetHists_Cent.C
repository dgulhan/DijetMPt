//=============================================                                                
// Author: Chris McGinn                                                                        
//                                                                                             
// DiJet Histogram Maker, Jet Properties                                                              
//                                                                                             
//=============================================     

#include "TTree.h"
#include "TDatime.h"
#include "TFile.h"
#include "diJetFileTag.h"

TFile* inFile_p = 0;
TFile* outFile_p = 0;

void makeCentHist(TTree* anaTree_p, const char* outName, Int_t setNum, sampleType sType = kHIDATA)
{
  inFile_p->cd();

  const char* title = Form("hiBin_%s_%s", algType[setNum], fileTag);

  TH1F* centHist_p;
  
  TString name = Form("%s_h(200, -.5, 199.5)", title);
  
  SetCuts(setNum, sType);

  if(sType == kHIDATA)  
    anaTree_p->Project(name, "hiBin", setCut && pthat);
  else
    anaTree_p->Project(name, "hiBin", Form("pthatWeight")*(setCut && pthat));  

  centHist_p = (TH1F*)inFile_p->Get(Form("%s_h", title));

  centHist_p->Scale(1./centHist_p->Integral());
  centHist_p->GetXaxis()->SetTitleOffset(.75);
  centHist_p->SetXTitle(title);

  outFile_p = new TFile(outName, "UPDATE");
  centHist_p->Write(Form("%s_h", title));
  outFile_p->Close();

  delete outFile_p;
  return;
}


void makeRatHist(const char* outName, Int_t setNum, sampleType sType = kHIDATA, const char* fileTag2 = "")
{
  TH1F* numHist_p;
  TH1F* denomHist_p;

  outFile_p = new TFile(outName, "READ");

  if(sType == kHIDATA){
    numHist_p = (TH1F*)outFile_p->Get(Form("hiBin_%s_%s_h", algType[setNum], fileTag));
    denomHist_p = (TH1F*)outFile_p->Get(Form("hiBin_%s_%s_h", algType[setNum], fileTag2));
  }
  else{
    numHist_p = (TH1F*)outFile_p->Get(Form("hiBin_%s_%s_h", algType[setNum], fileTag2));
    denomHist_p = (TH1F*)outFile_p->Get(Form("hiBin_%s_%s_h", algType[setNum], fileTag));
  }

  numHist_p->Divide(denomHist_p);
  numHist_p->SetXTitle(Form("hiBin_%s_DataOverMC_h", algType[setNum]));
  numHist_p->GetXaxis()->SetTitleOffset(.75);

  TFile* centFile_p = new TFile("centHist_eventSet.root", "UPDATE");
  numHist_p->Write(Form("hiBin_%s_DataOverMC_h", algType[setNum]));
  centFile_p->Close();
  delete centFile_p;

  outFile_p->Close();
  delete outFile_p;

  return;
}


void makeDiJetHists_Cent(const char* inName, const char* outName, sampleType sType = kHIDATA, const char* fileTag2 = "")
{
  TH1::SetDefaultSumw2();

  setFileTag(inName);

  inFile_p = new TFile(inName, "READ");
  TTree* anaTree_p = (TTree*)inFile_p->Get("jetTreeAna");

  Int_t jetAlgMax = 5;

  for(Int_t algIter = 0; algIter < jetAlgMax; algIter++){
    std::cout << "Algorithm: " << algType[algIter] << std::endl;

    if(algIter != 2){
      makeCentHist(anaTree_p, outName, algIter, sType);
      
      if(strcmp(fileTag2, "") != 0){
	makeRatHist(outName, algIter, sType, fileTag2);
      }
    }
  }

  inFile_p->Close();

  delete inFile_p;
  
  return;
}
