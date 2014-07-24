#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"

#include "TCanvas.h"
#include "TLine.h"

#include "/net/hisrv0001/home/cfmcginn/emDiJet/CMSSW_5_3_12_patch3/tempHIFA/HiForestAnalysis/commonSetup.h"
#include "cfmVectFunc.h"

#include <iostream>
#include <fstream>

#include "commonUtility.h"

TChain* getChain_p[3] = {0, 0, 0};

Float_t ptHatCuts[6] = {30, 50, 80, 120, 170, 1000000};
Float_t ptHatWeights[5] = {.556347, .030167, .00539882, .00546899, .0000263122};

Int_t pcollisionEventSelection_;
Float_t vz_;
Float_t pthat_;
Int_t hiBin_;

Int_t nref_;
Float_t jtpt_[maxJets];
Float_t jteta_[maxJets];
Float_t rawpt_[maxJets];
Float_t refpt_[maxJets];

Float_t jtPtCut = 0.0;

void BookChain(Bool_t isPbPb)
{
  getChain_p[0]->SetBranchStatus("*", 0);

  if(isPbPb) getChain_p[0]->SetBranchStatus("pcollisionEventSelection", 1);

  getChain_p[0]->SetBranchStatus("vz", 1);
  getChain_p[0]->SetBranchStatus("pthat", 1);
  getChain_p[0]->SetBranchStatus("hiBin", 1);

  getChain_p[0]->SetBranchStatus("nref", 1);
  getChain_p[0]->SetBranchStatus("jtpt", 1);
  getChain_p[0]->SetBranchStatus("jteta", 1);
  getChain_p[0]->SetBranchStatus("rawpt", 1);
  getChain_p[0]->SetBranchStatus("refpt", 1);

  if(isPbPb) getChain_p[0]->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection_);

  getChain_p[0]->SetBranchAddress("vz", &vz_);
  getChain_p[0]->SetBranchAddress("pthat", &pthat_);
  getChain_p[0]->SetBranchAddress("hiBin", &hiBin_);

  getChain_p[0]->SetBranchAddress("nref", &nref_);
  getChain_p[0]->SetBranchAddress("jtpt", &jtpt_);
  getChain_p[0]->SetBranchAddress("jteta", &jteta_);
  getChain_p[0]->SetBranchAddress("rawpt", &rawpt_);
  getChain_p[0]->SetBranchAddress("refpt", &refpt_);

  return;
}


void GetChain(std::vector<std::string> inList, const char* alg, Bool_t isPbPb)
{
  getChain_p[0] = new TChain(Form("%sJetAnalyzer/t", alg));
  getChain_p[1] = new TChain("skimanalysis/HltTree");
  getChain_p[2] = new TChain("hiEvtAnalyzer/HiTree");

  for(Int_t iter = 0; iter < (Int_t)(inList.size()); iter++){
    getChain_p[0]->Add(inList[iter].c_str());
    getChain_p[1]->Add(inList[iter].c_str());
    getChain_p[2]->Add(inList[iter].c_str());
  }
  getChain_p[0]->AddFriend(getChain_p[1]);
  getChain_p[0]->AddFriend(getChain_p[2]);

  BookChain(isPbPb);

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


void setHatWeights()
{
  Float_t denom = 0;
  for(Int_t iter = 0; iter < 5; iter++){
    denom += ptHatWeights[iter];
  }
  for(Int_t iter = 0; iter < 5; iter++){
    ptHatWeights[iter] = ptHatWeights[iter]/denom;
  }
  return;
}


Float_t getHatWeight(Float_t inHat)
{
  for(Int_t iter = 0; iter < 5; iter++){
    if(inHat > ptHatCuts[iter] && inHat < ptHatCuts[iter+1])
      return ptHatWeights[iter];
  }

  std::cout << "No weight assigned; check for error." << std::endl;
  return 0;
}


void makeDiJetIniHist(std::vector<std::string> inList, const char* outName, const char* alg = "akVs3Calo", Bool_t isPbPb = false)
{
  TH1::SetDefaultSumw2();

  std::cout << "Init: " << alg << std::endl;

  GetChain(inList, alg, isPbPb);

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
  rawOverGenHist_p->GetXaxis()->SetTitleOffset(.75);
  rawOverGenHist_p->SetYTitle("p_{T}^{raw}/p_{T}^{gen}");

  handsomeTH1(jetOverGenHist_p);
  jetOverGenHist_p->SetXTitle("p_{T}^{gen}");
  jetOverGenHist_p->GetXaxis()->SetTitleOffset(.75);
  jetOverGenHist_p->SetYTitle("p_{T}^{jet}/p_{T}^{gen}");
  jetOverGenHist_p->SetMaximum(1.10);
  jetOverGenHist_p->SetMinimum(0.90);

  for(Int_t jEntry = 0; jEntry <  getChain_p[0]->GetEntries(); jEntry++){
    getChain_p[0]->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(!pcollisionEventSelection_ && isPbPb) continue;

    if(TMath::Abs(vz_) > 15) continue;

    for(Int_t jtIter = 0; jtIter < nref_; jtIter++){
      if(TMath::Abs(jteta_[jtIter]) > 2.0 || refpt_[jtIter] < jtPtCut)	continue;

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

  CleanChain();

  return;
}




void makeDiJetIniHist_Eta(std::vector<std::string> inList, const char* outName, const char* alg = "akVs3Calo", Bool_t isPbPb = false)
{
  TH1::SetDefaultSumw2();

  std::cout << "Init: " << alg << std::endl;

  GetChain(inList, alg, isPbPb);

  const Int_t nBins = 25;
  Float_t bins[nBins+1];

  std::vector<Float_t>* mean_jetOverGen_p[nBins];

  bins[0] = -2.0;
  for(Int_t iter = 0; iter < nBins; iter++){
    bins[iter+1] = bins[iter] + .16;
    mean_jetOverGen_p[iter] = new std::vector<Float_t>;
  }

  TH1F* jetOverGenHist_p = new TH1F(Form("%s_eta_h", alg), Form("%s_eta_h", alg), nBins, bins);

  handsomeTH1(jetOverGenHist_p);
  jetOverGenHist_p->SetXTitle("#eta");
  jetOverGenHist_p->GetXaxis()->SetTitleOffset(.75);
  jetOverGenHist_p->SetYTitle("p_{T}^{jet}/p_{T}^{gen}");

  for(Int_t jEntry = 0; jEntry <  getChain_p[0]->GetEntries(); jEntry++){
    getChain_p[0]->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(!pcollisionEventSelection_ && isPbPb)
      continue;

    if(TMath::Abs(vz_) > 15)
      continue;

    for(Int_t jtIter = 0; jtIter < nref_; jtIter++){
      if(TMath::Abs(jteta_[jtIter]) > 2.0 || refpt_[jtIter]  < jtPtCut)
	continue;

      for(Int_t rawIter = 0; rawIter < nBins; rawIter++){
	if(jteta_[jtIter] > bins[rawIter] && jteta_[jtIter] < bins[rawIter + 1])
	  mean_jetOverGen_p[rawIter]->push_back(jtpt_[jtIter]/refpt_[jtIter]);
      }
    }
  }

  for(Int_t rawIter = 0; rawIter < nBins; rawIter++){
    if(mean_jetOverGen_p[rawIter]->size() != 0){
      Float_t mean = getMean(mean_jetOverGen_p[rawIter]);
      jetOverGenHist_p->SetBinContent(rawIter+1, mean);
      jetOverGenHist_p->SetBinError(rawIter+1, getError(mean_jetOverGen_p[rawIter], mean));
    }
  }

  TFile* out = new TFile(outName, "UPDATE");
  std::cout << outName << std::endl;
  jetOverGenHist_p->Write();
  out->Close();
  delete out;

  delete jetOverGenHist_p;
  jetOverGenHist_p = 0;

  for(Int_t rawIter = 0; rawIter < nBins; rawIter++){
    mean_jetOverGen_p[rawIter]->clear();
    delete mean_jetOverGen_p[rawIter];
    mean_jetOverGen_p[rawIter] = 0;
  }

  CleanChain();

  return;
}





void makeDiJetIniHist_Cent(std::vector<std::string> inList, const char* outName, const char* alg = "akVs3Calo", Bool_t isPbPb = false)
{
  TH1::SetDefaultSumw2();

  std::cout << "Init: " << alg << std::endl;

  GetChain(inList, alg, isPbPb);

  const Int_t nBins = 100;
  std::vector<Float_t>* mean_cent_p[nBins];
  for(Int_t iter = 0; iter < nBins; iter++){
    mean_cent_p[iter] = new std::vector<Float_t>;
  }

  TH1F* cent_p = new TH1F(Form("%s_cent_h", alg), Form("%s_cent_h", alg), nBins, -.5, 199.5);
  handsomeTH1(cent_p);
  cent_p->SetXTitle("hiBin");
  cent_p->GetXaxis()->SetTitleOffset(.75);
  cent_p->SetYTitle("p_{T}^{raw}/p_{T}^{gen}");

  for(Int_t jEntry = 0; jEntry < getChain_p[0]->GetEntries(); jEntry++){
    getChain_p[0]->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(!pcollisionEventSelection_ && isPbPb)
      continue;

    if(TMath::Abs(vz_) > 15)
      continue;

    Int_t binPos = cent_p->FindBin(hiBin_);

    for(Int_t jtIter = 0; jtIter < nref_; jtIter++){
      if(TMath::Abs(jteta_[jtIter]) > 2.0 || refpt_[jtIter]  < jtPtCut)
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

  CleanChain();

  return;
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

  return;
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

  return;
}


void drawLine(){
  TLine* zeroLine_p = new TLine(0.0, 1.0, 700.0, 1.0);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw("SAME");

  return;
}


void plotDiJetIniHistRatVsPu(const char* histFileName, const char* PFCalo = "PF")
{
  TH1::SetDefaultSumw2();

  TFile *f = new TFile(histFileName, "READ");
  TH1F* getHist_p = (TH1F*)f->Get(Form("ak3%s_rawVsOverPu_h", PFCalo));
  getHist_p->DrawCopy();
  drawLine();

  return;
}


void plotDiJetIniHistRatPFCalo(const char* histFileName, const char* VsPu = "Vs")
{
  TH1::SetDefaultSumw2();

  TFile *f = new TFile(histFileName, "READ");
  TH1F* getHist_p = (TH1F*)f->Get(Form("ak3%s_PFOverCalo_h", VsPu));
  getHist_p->DrawCopy();
  drawLine();

  return;
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
  get3Hist_p->SetXTitle("hiBin");
  get3Hist_p->GetXaxis()->SetTitleOffset(.75);
  get3Hist_p->SetYTitle("p^{raw}_{T}/p^{gen}_{T}");
  TCanvas* plotCanv_p = new TCanvas(Form("eqCheck_ak%s#%s", VsPu, PFCalo), Form("eqCheck_ak%s%s", VsPu, PFCalo), 1);
  get3Hist_p->DrawCopy();

  get4Hist_p->SetMarkerColor(kRed);
  get4Hist_p->SetLineColor(kRed);
  get4Hist_p->DrawCopy("SAME");

  get5Hist_p->SetMarkerColor(kBlue);
  get5Hist_p->SetLineColor(kBlue);
  get5Hist_p->DrawCopy("SAME");

  drawLine();

  TLegend* leg = new TLegend(.6, .6, .8, .8);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
  leg->SetTextSizePixels(28);
  leg->SetBorderSize(0);			

  leg->AddEntry(get3Hist_p, "R = 0.3");
  leg->AddEntry(get4Hist_p, "R = 0.4");
  leg->AddEntry(get5Hist_p, "R = 0.5");

  leg->Draw("SAME");

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  label_p->DrawLatex(.6, .82, Form("ak%s#%s", VsPu, PFCalo));  

  claverCanvasSaving(plotCanv_p, Form("pdfDir/eqCheck_ak%s%s", VsPu, PFCalo), "pdf");

  return;
}


int runMakeDiJetIniHist(std::string fList = "", const char* outFileName = "raw_rawOverGen", Bool_t isPbPb = false)
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

  setHatWeights();

  for(Int_t iter = 0; iter < 5; iter++){
    std::cout << ptHatWeights[iter] << std::endl;
  }

  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs3PF", isPbPb);
  //  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu3PF", isPbPb);
  //  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs3Calo", isPbPb);
  //  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu3Calo", isPbPb);

  /*
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs4PF", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu4PF", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs4Calo", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu4Calo", isPbPb);

  makeDiJetIniHist_Eta(listOfFiles, Form("%s.root", outFileName), "akVs4PF", isPbPb);
  makeDiJetIniHist_Eta(listOfFiles, Form("%s.root", outFileName), "akPu4PF", isPbPb);
  makeDiJetIniHist_Eta(listOfFiles, Form("%s.root", outFileName), "akVs4Calo", isPbPb);
  makeDiJetIniHist_Eta(listOfFiles, Form("%s.root", outFileName), "akPu4Calo", isPbPb);

  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs5PF", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu5PF", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs5Calo", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu5Calo", isPbPb);
  */
  /*
  makeDiJetIniHist_Cent(listOfFiles, Form("%s.root", outFileName), "akVs3PF", isPbPb);
  makeDiJetIniHist_Cent(listOfFiles, Form("%s.root", outFileName), "akPu3PF", isPbPb);
  makeDiJetIniHist_Cent(listOfFiles, Form("%s.root", outFileName), "akVs3Calo", isPbPb);
  makeDiJetIniHist_Cent(listOfFiles, Form("%s.root", outFileName), "akPu3Calo", isPbPb);
  */
  /*
  makeDiJetIniHist_Cent(listOfFiles, Form("%s.root", outFileName), "akVs4PF", isPbPb);
  makeDiJetIniHist_Cent(listOfFiles, Form("%s.root", outFileName), "akPu4PF", isPbPb);
  makeDiJetIniHist_Cent(listOfFiles, Form("%s.root", outFileName), "akVs4Calo", isPbPb);
  makeDiJetIniHist_Cent(listOfFiles, Form("%s.root", outFileName), "akPu4Calo", isPbPb);

  makeDiJetIniHist_Cent(listOfFiles, Form("%s.root", outFileName), "akVs5PF", isPbPb);
  makeDiJetIniHist_Cent(listOfFiles, Form("%s.root", outFileName), "akPu5PF", isPbPb);
  makeDiJetIniHist_Cent(listOfFiles, Form("%s.root", outFileName), "akVs5Calo", isPbPb);
  makeDiJetIniHist_Cent(listOfFiles, Form("%s.root", outFileName), "akPu5Calo", isPbPb);
  */

  /*
  makeDiJetIniHistRatVsPu(Form("%s.root", outFileName), "PF");
  makeDiJetIniHistRatVsPu(Form("%s.root", outFileName), "Calo");

  makeDiJetIniHistRatPFCalo(Form("%s.root", outFileName), "Vs");
  makeDiJetIniHistRatPFCalo(Form("%s.root", outFileName), "Pu");
  */

  return(0);
}


int main(int argc, char* argv[])
{
  if(argc != 4){
    std::cout << "Usage: runMakeDiJetIniHist <inputList> <outFileName> <isPbPb>" << std::endl;
    return 1;
  }

  int rStatus = -1;

  rStatus = runMakeDiJetIniHist(argv[1], argv[2], (Bool_t)(atoi(argv[3])));

  return rStatus;
}
