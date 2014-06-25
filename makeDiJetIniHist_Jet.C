#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCut.h"
#include "TMath.h"

#include "TCanvas.h"
#include "TLine.h"

#include "/net/hisrv0001/home/cfmcginn/emDiJet/CMSSW_5_3_12_patch3/tempHIFA/HiForestAnalysis/commonSetup.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "commonUtility.h"

TTree* getTree_p = 0;

Int_t pcollisionEventSelection_;
Float_t vz_;
Float_t pthat_;
Int_t hiBin_;

Int_t nref_;
Float_t jtpt_[maxJets];
Float_t jteta_[maxJets];
Float_t rawpt_[maxJets];
Float_t refpt_[maxJets];

void Book_Tree()
{
  getTree_p->SetBranchStatus("*", 0);
  getTree_p->SetBranchStatus("pcollisionEventSelection", 1);
  getTree_p->SetBranchStatus("vz", 1);
  getTree_p->SetBranchStatus("pthat", 1);
  getTree_p->SetBranchStatus("hiBin", 1);

  getTree_p->SetBranchStatus("nref", 1);
  getTree_p->SetBranchStatus("jtpt", 1);
  getTree_p->SetBranchStatus("jteta", 1);
  getTree_p->SetBranchStatus("rawpt", 1);
  getTree_p->SetBranchStatus("refpt", 1);

  getTree_p->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection_);
  getTree_p->SetBranchAddress("vz", &vz_);
  getTree_p->SetBranchAddress("pthat", &pthat_);
  getTree_p->SetBranchAddress("hiBin", &hiBin_);

  getTree_p->SetBranchAddress("nref", &nref_);
  getTree_p->SetBranchAddress("jtpt", &jtpt_);
  getTree_p->SetBranchAddress("jteta", &jteta_);
  getTree_p->SetBranchAddress("rawpt", &rawpt_);
  getTree_p->SetBranchAddress("refpt", &refpt_);
}


Float_t getMean(std::vector<Float_t>* inVect_p)
{
  if(inVect_p->size() == 0){
    std::cout << "Passed empty vector. Return 0." << std::endl;
    return 0;
  }

  Float_t sum = 0;
  for(Int_t jtIter = 0; jtIter < (Int_t)inVect_p->size(); jtIter++){
    sum += inVect_p->at(jtIter);
  }

  return sum/((Float_t)inVect_p->size());
}


Float_t getError(std::vector<Float_t>* inVect_p, Float_t mean)
{
  if(inVect_p->size() < 2){
    std::cout << "Passed vector of size 1 or less. Return 0." << std::endl;
    return 0;
  }

  Float_t error = 0;
  for(Int_t jtIter = 0; jtIter < (Int_t)inVect_p->size(); jtIter++){
    error += TMath::Power((inVect_p->at(jtIter) - mean), 2);
  }

  error = TMath::Sqrt(error/((Float_t)inVect_p->size() - 1.0));
  return error/TMath::Sqrt((Float_t)inVect_p->size());
}


void makeDiJetIniHist(std::vector<std::string> inList, Int_t num, const char* outName, const char* alg = "akVs3Calo")
{
  TH1::SetDefaultSumw2();

  std::cout << "Init: " << alg << std::endl;

  TFile* f = new TFile(inList[num].c_str(), "READ");
  getTree_p = (TTree*)f->Get(Form("%sJetAnalyzer/t", alg));
  getTree_p->AddFriend("skimanalysis/HltTree");
  getTree_p->AddFriend("hiEvtAnalyzer/HiTree");

  Book_Tree();

  const Int_t nBins = 100;
  Float_t bins[nBins+1];

  std::vector<Float_t>* mean_rawOverGen_p[nBins];
  std::vector<Float_t>* mean_jetOverGen_p[nBins];

  bins[0] = 0;
  for(Int_t iter = 0; iter < nBins; iter++){
    bins[iter+1] = bins[iter] + 7;
    mean_rawOverGen_p[iter] = new std::vector<Float_t>;
    mean_jetOverGen_p[iter] = new std::vector<Float_t>;
  }

  TH1F* rawOverGenHist_p = new TH1F(Form("%s_rawOverGen_h", alg), Form("%s_rawOverGen_h", alg), nBins, bins);
  TH1F* jetOverGenHist_p = new TH1F(Form("%s_jetOverGen_h", alg), Form("%s_jetOverGen_h", alg), nBins, bins);

  handsomeTH1(rawOverGenHist_p);
  rawOverGenHist_p->SetXTitle("p_{T}^{raw}");
  rawOverGenHist_p->SetYTitle("p_{T}^{raw}/p_{T}^{gen}");

  handsomeTH1(jetOverGenHist_p);
  jetOverGenHist_p->SetXTitle("p_{T}^{gen}");
  jetOverGenHist_p->SetYTitle("p_{T}^{jet}/p_{T}^{gen}");

  for(Int_t jEntry = 0; jEntry <  getTree_p->GetEntries(); jEntry++){
    getTree_p->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(!pcollisionEventSelection_)
      continue;

    if(TMath::Abs(vz_) > 15)
      continue;

    for(Int_t jtIter = 0; jtIter < nref_; jtIter++){
      if(TMath::Abs(jteta_[jtIter]) > 2.0 || refpt_[jtIter]  < 0)
	continue;

      for(Int_t rawIter = 0; rawIter < nBins; rawIter++){
	if(rawpt_[jtIter] > bins[rawIter] && rawpt_[jtIter] < bins[rawIter + 1])
	  mean_rawOverGen_p[rawIter]->push_back(rawpt_[jtIter]/refpt_[jtIter]);

	if(refpt_[jtIter] > bins[rawIter] && refpt_[jtIter] < bins[rawIter + 1])
	  mean_jetOverGen_p[rawIter]->push_back(jtpt_[jtIter]/refpt_[jtIter]);
      }
    }
  }

  for(Int_t rawIter = 0; rawIter < nBins; rawIter++){
    if(mean_rawOverGen_p[rawIter]->size() != 0){
      Float_t mean = getMean(mean_rawOverGen_p[rawIter]);
      rawOverGenHist_p->SetBinContent(rawIter+1, mean);
      rawOverGenHist_p->SetBinError(rawIter+1, getError(mean_rawOverGen_p[rawIter], mean));
    }

    if(mean_jetOverGen_p[rawIter]->size() != 0){
      Float_t mean = getMean(mean_jetOverGen_p[rawIter]);
      jetOverGenHist_p->SetBinContent(rawIter+1, mean);
      jetOverGenHist_p->SetBinError(rawIter+1, getError(mean_jetOverGen_p[rawIter], mean));
    }
  }

  TFile* out = new TFile(outName, "UPDATE");
  std::cout << outName << std::endl;
  rawOverGenHist_p->Write();
  jetOverGenHist_p->Write();
  out->Close();
  delete out;

  delete rawOverGenHist_p;
  rawOverGenHist_p = 0;

  delete jetOverGenHist_p;
  jetOverGenHist_p = 0;

  for(Int_t rawIter = 0; rawIter < nBins; rawIter++){
    mean_rawOverGen_p[rawIter]->clear();
    delete mean_rawOverGen_p[rawIter];
    mean_rawOverGen_p[rawIter] = 0;

    mean_jetOverGen_p[rawIter]->clear();
    delete mean_jetOverGen_p[rawIter];
    mean_jetOverGen_p[rawIter] = 0;
  }

  delete f;
}


void makeDiJetIniHist_Cent(std::vector<std::string> inList, Int_t num, const char* outName, const char* alg = "akVs3Calo")
{
  TH1::SetDefaultSumw2();

  std::cout << "Init: " << alg << std::endl;

  TFile* f = new TFile(inList[num].c_str(), "READ");
  getTree_p = (TTree*)f->Get(Form("%sJetAnalyzer/t", alg));
  getTree_p->AddFriend("skimanalysis/HltTree");
  getTree_p->AddFriend("hiEvtAnalyzer/HiTree");

  Book_Tree();

  const Int_t nBins = 100;
  std::vector<Float_t>* mean_cent_p[nBins];
  for(Int_t iter = 0; iter < nBins; iter++){
    mean_cent_p[iter] = new std::vector<Float_t>;
  }

  TH1F* cent_p = new TH1F(Form("%s_cent_h", alg), Form("%s_cent_h", alg), nBins, -.5, 199.5);
  handsomeTH1(cent_p);
  cent_p->SetXTitle("hiBin");
  cent_p->SetYTitle("p_{T}^{raw}/p_{T}^{gen}");

  for(Int_t jEntry = 0; jEntry < getTree_p->GetEntries(); jEntry++){
    getTree_p->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(!pcollisionEventSelection_)
      continue;

    if(TMath::Abs(vz_) > 15)
      continue;

    Int_t binPos = cent_p->FindBin(hiBin_);

    for(Int_t jtIter = 0; jtIter < nref_; jtIter++){
      if(TMath::Abs(jteta_[jtIter]) > 2.0 || refpt_[jtIter]  < 0)
        continue;

      mean_cent_p[binPos]->push_back(rawpt_[jtIter]/refpt_[jtIter]);
    }
  }

  for(Int_t binIter = 0; binIter < nBins; binIter++){
    if(mean_cent_p[binIter]->size() != 0){
      Float_t mean = getMean(mean_cent_p[binIter]);
      cent_p->SetBinContent(binIter + 1, mean);
      cent_p->SetBinError(binIter + 1, getError(mean_cent_p[binIter], mean));
    }
  }

  TFile* out = new TFile(outName, "UPDATE");
  std::cout << outName << std::endl;
  cent_p->Write();
  out->Close();
  delete out;
  out = 0;

  for(Int_t centIter = 0; centIter < nBins; centIter++){
    mean_cent_p[centIter]->clear();
    delete mean_cent_p[centIter];
    mean_cent_p[centIter] = 0;
  }

  delete f;
}


void makeDiJetIniHistRatVsPu(const char* histFileName, const char* PFCalo = "PF")
{
  TH1::SetDefaultSumw2();

  TFile* f = new TFile(histFileName, "UPDATE");

  TH1F* getRawVs_p = (TH1F*)f->Get(Form("akVs3%s_rawOverGen_h", PFCalo));
  TH1F* getRawPu_p = (TH1F*)f->Get(Form("akPu3%s_rawOverGen_h", PFCalo));
  TH1F* getJetVs_p = (TH1F*)f->Get(Form("akVs3%s_jetOverGen_h", PFCalo));
  TH1F* getJetPu_p = (TH1F*)f->Get(Form("akPu3%s_jetOverGen_h", PFCalo));

  getRawVs_p->Divide(getRawPu_p);
  getJetVs_p->Divide(getJetPu_p);

  getRawVs_p->SetYTitle("Vs/Pu p_{T}^{rat}");
  getRawVs_p->Write(Form("ak3%s_rawVsOverPu_h", PFCalo));
  getJetVs_p->SetYTitle("Vs/Pu p_{T}^{rat}");
  getJetVs_p->Write(Form("ak3%s_jetVsOverPu_h", PFCalo));

  getRawPu_p = 0;
  getRawVs_p = 0;
  getJetPu_p = 0;
  getJetVs_p = 0;
  f->Close();
  delete f;
}


void makeDiJetIniHistRatPFCalo(const char* histFileName, const char* VsPu = "Vs")
{
  TH1::SetDefaultSumw2();

  TFile* f = new TFile(histFileName, "UPDATE");

  TH1F* getRawPF_p = (TH1F*)f->Get(Form("ak%s3PF_rawOverGen_h", VsPu));
  TH1F* getRawCalo_p = (TH1F*)f->Get(Form("ak%s3Calo_rawOverGen_h", VsPu));
  TH1F* getJetPF_p = (TH1F*)f->Get(Form("ak%s3PF_jetOverGen_h", VsPu));
  TH1F* getJetCalo_p = (TH1F*)f->Get(Form("ak%s3Calo_jetOverGen_h", VsPu));

  getRawPF_p->Divide(getRawCalo_p);
  getJetPF_p->Divide(getJetCalo_p);

  getRawPF_p->SetYTitle("PF/Calo p_{T}^{rat}");
  getRawPF_p->Write(Form("ak3%s_rawPFOverCalo_h", VsPu));
  getJetPF_p->SetYTitle("PF/Calo p_{T}^{rat}");
  getJetPF_p->Write(Form("ak3%s_jetPFOverCalo_h", VsPu));

  getRawCalo_p = 0;
  getRawPF_p = 0;
  getJetCalo_p = 0;
  getJetPF_p = 0;
  f->Close();
  delete f;
}


void drawLine(){
  TLine* zeroLine_p = new TLine(0.0, 1.0, 700.0, 1.0);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw("SAME");
}


void plotDiJetIniHistRatVsPu(const char* histFileName, const char* PFCalo = "PF")
{
  TH1::SetDefaultSumw2();

  TFile *f = new TFile(histFileName, "READ");
  TH1F* getHist_p = (TH1F*)f->Get(Form("ak3%s_rawVsOverPu_h", PFCalo));
  getHist_p->DrawCopy();
  drawLine();
}


void plotDiJetIniHistRatPFCalo(const char* histFileName, const char* VsPu = "Vs")
{
  TH1::SetDefaultSumw2();

  TFile *f = new TFile(histFileName, "READ");
  TH1F* getHist_p = (TH1F*)f->Get(Form("ak3%s_PFOverCalo_h", VsPu));
  getHist_p->DrawCopy();
  drawLine();
}


void plotDiJetIniHistCent(const char* histFileName, const char* VsPu = "Vs", const char* PFCalo = "Calo")
{
  TH1::SetDefaultSumw2();

  TFile *f = new TFile(histFileName, "READ");

  TH1F* get3Hist_p = (TH1F*)f->Get(Form("ak%s3%s_cent_h", VsPu, PFCalo));
  TH1F* get4Hist_p = (TH1F*)f->Get(Form("ak%s4%s_cent_h", VsPu, PFCalo));
  TH1F* get5Hist_p = (TH1F*)f->Get(Form("ak%s5%s_cent_h", VsPu, PFCalo));


  get3Hist_p->SetMaximum(1.0);
  get3Hist_p->SetMinimum(0.6);
  get3Hist_p->DrawCopy();

  get4Hist_p->SetMarkerColor(kRed);
  get4Hist_p->SetLineColor(kRed);
  get4Hist_p->DrawCopy("SAME");

  get5Hist_p->SetMarkerColor(kBlue);
  get5Hist_p->SetLineColor(kBlue);
  get5Hist_p->DrawCopy("SAME");

  drawLine();

  TLegend* leg = new TLegend(.6, .4, .8, .6);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
  leg->SetTextSizePixels(28);
  leg->SetBorderSize(0);			

  leg->AddEntry(get3Hist_p, "R = 0.3");
  leg->AddEntry(get4Hist_p, "R = 0.4");
  leg->AddEntry(get5Hist_p, "R = 0.5");

  leg->Draw("SAME");
}


int runMakeDiJetIniHist(std::string fList = "", const char* outFileName = "raw_rawOverGen", Int_t num = 0)
{
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
    while(!inFile.eof()){
      inFile >> buffer;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }

  makeDiJetIniHist(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akVs3PF");
  makeDiJetIniHist(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akPu3PF");
  makeDiJetIniHist(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akVs3Calo");
  makeDiJetIniHist(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akPu3Calo");

  makeDiJetIniHist_Cent(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akVs3PF");
  makeDiJetIniHist_Cent(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akPu3PF");
  makeDiJetIniHist_Cent(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akVs3Calo");
  makeDiJetIniHist_Cent(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akPu3Calo");

  makeDiJetIniHist_Cent(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akVs4PF");
  makeDiJetIniHist_Cent(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akPu4PF");
  makeDiJetIniHist_Cent(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akVs4Calo");
  makeDiJetIniHist_Cent(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akPu4Calo");

  makeDiJetIniHist_Cent(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akVs5PF");
  makeDiJetIniHist_Cent(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akPu5PF");
  makeDiJetIniHist_Cent(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akVs5Calo");
  makeDiJetIniHist_Cent(listOfFiles, num, Form("%s_%d.root", outFileName, num), "akPu5Calo");

  makeDiJetIniHistRatVsPu(Form("%s_%d.root", outFileName, num), "PF");
  makeDiJetIniHistRatVsPu(Form("%s_%d.root", outFileName, num), "Calo");

  makeDiJetIniHistRatPFCalo(Form("%s_%d.root", outFileName, num), "Vs");
  makeDiJetIniHistRatPFCalo(Form("%s_%d.root", outFileName, num), "Pu");

  return(0);
}


int main(int argc, char* argv[])
{
  if(argc != 4){
    std::cout << "Usage: runMakeDiJetIniHist <inputList> <outFileName> <#>" << std::endl;
    return 1;
  }

  int rStatus = -1;

  rStatus = runMakeDiJetIniHist(argv[1], argv[2], atoi(argv[3]));

  return rStatus;
}
