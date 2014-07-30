#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"

#include "TCanvas.h"
#include "TLine.h"

#include "/net/hisrv0001/home/cfmcginn/emDiJet/CMSSW_5_3_12_patch3/tempHIFA/HiForestAnalysis/commonSetup.h"

#include <iostream>
#include <fstream>
#include <string>

#include "commonUtility.h"
#include "TMath.h"
#include "TLatex.h"
TChain* getChain_p[3] = {0, 0, 0};

Int_t ptHatCuts_PYTH[6] = {30, 50, 80, 120, 170, 1000000};
Float_t ptHatWeights_PYTH[5] = {.556347, .057268, .00590566, .00597628, .0000266327};
Int_t ptHatCuts_PYTHHYD[9] = {15, 30, 50, 80, 120, 220, 280, 370, 1000000};
Float_t ptHatWeights_PYTHHYD[8] = {.611066, .0399951, .00243874, .000241009, .0000273228, .00000147976, .000000618337, .000000250369};

Int_t centCutArray[9] = {0, 10, 20, 40, 60, 80, 100, 140, 200};

Int_t pcollisionEventSelection_;
Float_t vz_;
Float_t pthat_;
Int_t hiBin_;

Int_t nref_;
Float_t jtpt_[maxJets];
Float_t jteta_[maxJets];
Float_t rawpt_[maxJets];
Float_t refpt_[maxJets];

Float_t jtPtCut = 20.0;

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
  Float_t denom_PbPb = 0;
  for(Int_t iter = 0; iter < 5; iter++){
    denom += ptHatWeights_PYTH[iter];
  }
  for(Int_t iter = 0; iter < 5; iter++){
    ptHatWeights_PYTH[iter] = ptHatWeights_PYTH[iter]/denom;
  }

  for(Int_t iter = 0; iter < 5; iter++){
    denom_PbPb += ptHatWeights_PYTHHYD[iter];
  }
  for(Int_t iter = 0; iter < 8; iter++){
    ptHatWeights_PYTHHYD[iter] = ptHatWeights_PYTHHYD[iter]/denom_PbPb;
  }

  return;
}


Float_t getHatWeight(Float_t inHat, Bool_t isPbPb)
{
  if(isPbPb){
    for(Int_t iter = 0; iter < 8; iter++){
      if(inHat > ptHatCuts_PYTHHYD[iter] && inHat < ptHatCuts_PYTHHYD[iter+1])
	return ptHatWeights_PYTHHYD[iter];
    }
  }
  else{
    for(Int_t iter = 0; iter < 5; iter++){
      if(inHat > ptHatCuts_PYTH[iter] && inHat < ptHatCuts_PYTH[iter+1])
	return ptHatWeights_PYTH[iter];
    }
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


std::string getCentString(Bool_t isPbPb, Int_t centLow, Int_t centHi)
{
  if(isPbPb) return Form("%d%d", (Int_t)(centLow*.5), (Int_t)((centHi)*.5));
  else return "PP";
}


Int_t getCentPos(Int_t inHiBin)
{
  for(Int_t centIter = 0; centIter < 9; centIter++){
    if(inHiBin >= centCutArray[centIter] && inHiBin < centCutArray[centIter+1])
      return centIter;
  }
  std::cout << "Uh oh, no known centrality, getCentPos, return -1" << std::endl;
  return -1;
}


void makeDiJetIniHist(std::vector<std::string> inList, const std::string outName, const char* alg = "akVs3Calo", Bool_t isPbPb = false)
{
  TH1::SetDefaultSumw2();

  std::cout << "Init: " << alg << std::endl;

  GetChain(inList, alg, isPbPb);

  const Int_t nBins = 50;
  const Float_t lower = 20.;
  const Float_t higher = 700.;

  Int_t nCentBins;

  if(isPbPb) nCentBins = 8;
  else nCentBins = 1;

  std::string centString[nCentBins];

  if(isPbPb){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      centString[centIter] = getCentString(isPbPb, centCutArray[centIter], centCutArray[centIter+1]);
    }
  }
  else centString[0] = getCentString(isPbPb, 0, 0);

  //Edit Here

  Float_t bins[nBins+1];

  getLogBins(lower, higher, nBins, bins);

  TH1F* pthatNoWeight_p = new TH1F("pthatNoWeight_h", "pthatNoWeight_h", 100, 30, 300);
  TH1F* pthatWeight_p = new TH1F("pthatWeight_h", "pthatWeight_h", 100, 30, 300);
  TH1F* mean_jetOverGen_p[nCentBins][nBins];

  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    for(Int_t iter = 0; iter < nBins; iter++){
      mean_jetOverGen_p[centIter][iter] = new TH1F(Form("tempJetHist_%s_%d", centString[centIter].c_str(), iter), Form("tempJetHist_%s_%d", centString[centIter].c_str(), iter), 2000, 0, 2);
    }
  }

  TH1F* jetOverGenHist_p[nCentBins];
  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    jetOverGenHist_p[centIter] = new TH1F(Form("%s_jetOverGen_%s_h", alg, centString[centIter].c_str()), Form("%s_jetOverGen_%s_h", alg, centString[centIter].c_str()), nBins, bins);
    handsomeTH1(jetOverGenHist_p[centIter]);
    jetOverGenHist_p[centIter]->SetXTitle("p_{T}^{gen}");
    jetOverGenHist_p[centIter]->GetXaxis()->SetTitleOffset(.75);
    jetOverGenHist_p[centIter]->SetYTitle("p_{T}^{jet}/p_{T}^{gen}");
    jetOverGenHist_p[centIter]->SetMaximum(1.04999);
    jetOverGenHist_p[centIter]->SetMinimum(0.95001);
  }

  for(Int_t jEntry = 0; jEntry <  getChain_p[0]->GetEntries(); jEntry++){
    getChain_p[0]->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(!pcollisionEventSelection_ && isPbPb) continue;

    if(TMath::Abs(vz_) > 15) continue;

    if(pthat_ < 30) continue;

    Float_t hatWeight = getHatWeight(pthat_, isPbPb);
    Int_t centPos = 0;
    if(isPbPb) centPos = getCentPos(hiBin_);

    pthatNoWeight_p->Fill(pthat_);
    pthatWeight_p->Fill(pthat_, hatWeight);

    for(Int_t jtIter = 0; jtIter < nref_; jtIter++){
      if(TMath::Abs(jteta_[jtIter]) > 2.0 || refpt_[jtIter] < jtPtCut)	continue;

      for(Int_t rawIter = 0; rawIter < nBins; rawIter++){
	if(refpt_[jtIter] > bins[rawIter] && refpt_[jtIter] < bins[rawIter + 1]){
	  mean_jetOverGen_p[centPos][rawIter]->Fill(jtpt_[jtIter]/refpt_[jtIter], hatWeight);
	  break;
	}
      }
    }
  }


  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    for(Int_t rawIter = 0; rawIter < nBins; rawIter++){
      if(mean_jetOverGen_p[centIter][rawIter]->GetEntries() != 0){
	mean_jetOverGen_p[centIter][rawIter]->Fit("gaus", "Q WL");
	jetOverGenHist_p[centIter]->SetBinContent(rawIter+1, mean_jetOverGen_p[centIter][rawIter]->GetFunction("gaus")->GetParameter(1));
	
	if(mean_jetOverGen_p[centIter][rawIter]->GetEntries() != 1){
	  jetOverGenHist_p[centIter]->SetBinError(rawIter+1, mean_jetOverGen_p[centIter][rawIter]->GetFunction("gaus")->GetParError(1));
	  
	  if(jetOverGenHist_p[centIter]->GetBinError(rawIter+1) > .05){
	    std::cout << rawIter << ", " << bins[rawIter] << ", " << bins[rawIter+1] << ", " << jetOverGenHist_p[centIter]->GetBinContent(rawIter + 1) << ", " << jetOverGenHist_p[centIter]->GetBinError(rawIter+1) << std::endl;
	  }
	}
      }    
    }
  }


  TFile* out = new TFile(outName.c_str(), "UPDATE");
  std::cout << outName << std::endl;
  if(!strcmp(alg, "akVs3PF")){
    pthatWeight_p->Scale(1./pthatWeight_p->Integral());
    handsomeTH1(pthatWeight_p);
    pthatWeight_p->GetXaxis()->SetTitleOffset(1.15);
    pthatWeight_p->SetXTitle("p_{T}^{#hat}");
    pthatWeight_p->GetXaxis()->SetTitleOffset(1.15);
    pthatWeight_p->SetYTitle("Event Fraction");
    pthatWeight_p->Write();

    pthatNoWeight_p->Scale(1./pthatNoWeight_p->Integral());
    handsomeTH1(pthatNoWeight_p);
    pthatNoWeight_p->GetXaxis()->SetTitleOffset(1.15);
    pthatNoWeight_p->SetXTitle("p_{T}^{#hat}");
    pthatNoWeight_p->GetXaxis()->SetTitleOffset(1.15);
    pthatNoWeight_p->SetYTitle("Event Fraction");
    pthatNoWeight_p->Write();
  }
  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    jetOverGenHist_p[centIter]->Write();
  }
  out->Close();
  delete out;

  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    delete jetOverGenHist_p[centIter];
    jetOverGenHist_p[centIter] = 0;
  }

  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    for(Int_t rawIter = 0; rawIter < nBins; rawIter++){
      delete mean_jetOverGen_p[centIter][rawIter];
      mean_jetOverGen_p[centIter][rawIter] = 0;   
    }
  }

  delete pthatWeight_p;
  pthatWeight_p = 0;

  delete pthatNoWeight_p;
  pthatNoWeight_p = 0;

  CleanChain();

  return;
}



/*
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

*/
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

  TLine* fiftyLine_p = new TLine(50.0, 0.95, 50.0, 1.05);
  fiftyLine_p->SetLineColor(1);
  fiftyLine_p->SetLineStyle(2);
  fiftyLine_p->Draw("SAME");

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


void plotInitHist(TH1F* inHist_p)
{
  inHist_p->GetXaxis()->SetTitleSize(0.07);
  inHist_p->GetXaxis()->SetTitleOffset(1.1);
  inHist_p->GetXaxis()->SetLabelSize(0.06);
  inHist_p->GetYaxis()->SetTitleSize(0.07);
  inHist_p->GetYaxis()->SetTitleOffset(1.2);
  inHist_p->GetYaxis()->SetLabelSize(0.055);
}


void plotDiJetIniHist_PYTH(const std::string histFileName, const std::string VsPu, const std::string PFCalo)
{
  TH1::SetDefaultSumw2();
  TFile *f = new TFile(histFileName.c_str(), "UPDATE");

  TH1F* getHist3_p = (TH1F*)f->Get(Form("ak%s3%s_jetOverGen_PP_h", VsPu.c_str(), PFCalo.c_str()));
  TH1F* getHist4_p = (TH1F*)f->Get(Form("ak%s4%s_jetOverGen_PP_h", VsPu.c_str(), PFCalo.c_str()));
  TH1F* getHist5_p = (TH1F*)f->Get(Form("ak%s5%s_jetOverGen_PP_h", VsPu.c_str(), PFCalo.c_str()));

  TCanvas* plotCanv_p = new TCanvas(Form("ak%s%s_jetOverGen_PP_c", VsPu.c_str(), PFCalo.c_str()), Form("ak%s%s_jetOverGen_PP_c", VsPu.c_str(), PFCalo.c_str()), 3*300, 1*350);
  plotCanv_p->Divide(3, 1, 0.0, 0.0);

  plotInitHist(getHist3_p);
  plotInitHist(getHist4_p);
  plotInitHist(getHist5_p);

  plotCanv_p->cd(1);
  gPad->SetLogx();
  getHist3_p->Draw();
  drawLine();

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(23);

  label_p->DrawLatex(.5, .9, Form("ak%s3%s", VsPu.c_str(), PFCalo.c_str()));
  label_p->DrawLatex(.5, .80, Form("p_{T}^{gen} > 20 GeV/c"));


  plotCanv_p->cd(2);
  gPad->SetLogx();
  getHist4_p->Draw();
  drawLine();
  label_p->DrawLatex(.5, .9, Form("ak%s4%s", VsPu.c_str(), PFCalo.c_str()));
  label_p->DrawLatex(.5, .80, Form("|#eta| < 2.0"));

  plotCanv_p->cd(3);
  gPad->SetLogx();
  getHist5_p->Draw();
  drawLine();
  label_p->DrawLatex(.6, .9, Form("ak%s5%s", VsPu.c_str(), PFCalo.c_str()));
  label_p->DrawLatex(.6, .80, "PYTHIA");

  plotCanv_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(plotCanv_p, Form("pdfDir/ak%s%s_jetOverGen", VsPu.c_str(), PFCalo.c_str()), "pdf");

  delete label_p;
  delete plotCanv_p;
  f->Close();
  delete f;
}


void drawCentHist(TH1F* inHist_p)
{
  gPad->SetLogx();
  inHist_p->Draw();
  drawLine();
  return;
}


void plotDiJetIniHist_PYTHHYD(const std::string histFileName, const std::string alg)
{
  TH1::SetDefaultSumw2();
  TFile *f = new TFile(histFileName.c_str(), "UPDATE");

  TH1F* getHist_p[8];

  for(Int_t iter = 0; iter < 8; iter++){
    getHist_p[iter] = (TH1F*)f->Get(Form("%s_jetOverGen_%d%d_h", alg.c_str(), centCutArray[iter]/2, centCutArray[iter+1]/2));
  }

  TCanvas* plotCanv_p = new TCanvas(Form("%s_jetOverGen_PbPb_c", alg.c_str()), Form("%s_jetOverGen_PbPb_c", alg.c_str()), 4*300, 2*350);
  plotCanv_p->Divide(4, 2, 0.0, 0.0);

  for(Int_t iter = 0; iter < 8; iter++){
    plotInitHist(getHist_p[iter]);
  }

  plotCanv_p->cd(1);
  drawCentHist(getHist_p[7]);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(23);

  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[7]/2, centCutArray[8]/2));
  label_p->DrawLatex(.5, .70, Form("p_{T}^{gen} > 20 GeV/c"));

  plotCanv_p->cd(2);
  drawCentHist(getHist_p[6]);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[6]/2, centCutArray[7]/2));
  label_p->DrawLatex(.5, .70, Form("|#eta| < 2.0"));

  plotCanv_p->cd(3);
  drawCentHist(getHist_p[5]);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[5]/2, centCutArray[6]/2));
  label_p->DrawLatex(.5, .70, "PYT+HYD");

  plotCanv_p->cd(4);
  drawCentHist(getHist_p[4]);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[4]/2, centCutArray[5]/2));

  plotCanv_p->cd(5);
  drawCentHist(getHist_p[3]);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[3]/2, centCutArray[4]/2));

  plotCanv_p->cd(6);
  drawCentHist(getHist_p[2]);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[2]/2, centCutArray[3]/2));

  plotCanv_p->cd(7);
  drawCentHist(getHist_p[1]);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[1]/2, centCutArray[2]/2));

  plotCanv_p->cd(8);
  drawCentHist(getHist_p[0]);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[0]/2, centCutArray[1]/2));


  plotCanv_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(plotCanv_p, Form("pdfDir/%s_jetOverGen", alg.c_str()), "pdf");

  delete label_p;
  delete plotCanv_p;
  f->Close();
  delete f;

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


void runPlotDiJetIniHist_PYTH(const std::string histFileName)
{
  plotDiJetIniHist_PYTH(histFileName, "", "PF");
  plotDiJetIniHist_PYTH(histFileName, "", "Calo");

  plotDiJetIniHist_PYTH(histFileName, "Vs", "PF");
  plotDiJetIniHist_PYTH(histFileName, "Vs", "Calo");

  plotDiJetIniHist_PYTH(histFileName, "Pu", "PF");
  plotDiJetIniHist_PYTH(histFileName, "Pu", "Calo");

  return;
}


void runPlotDiJetIniHist_PYTHHYD(const std::string histFileName)
{
  plotDiJetIniHist_PYTHHYD(histFileName, "akVs3PF");
  plotDiJetIniHist_PYTHHYD(histFileName, "akVs3Calo");

  plotDiJetIniHist_PYTHHYD(histFileName, "akPu3PF");
  plotDiJetIniHist_PYTHHYD(histFileName, "akPu3Calo");

  plotDiJetIniHist_PYTHHYD(histFileName, "akVs4PF");
  plotDiJetIniHist_PYTHHYD(histFileName, "akVs4Calo");

  plotDiJetIniHist_PYTHHYD(histFileName, "akPu4PF");
  plotDiJetIniHist_PYTHHYD(histFileName, "akPu4Calo");

  plotDiJetIniHist_PYTHHYD(histFileName, "akVs5PF");
  plotDiJetIniHist_PYTHHYD(histFileName, "akVs5Calo");

  plotDiJetIniHist_PYTHHYD(histFileName, "akPu5PF");
  plotDiJetIniHist_PYTHHYD(histFileName, "akPu5Calo");

  return;
}


int runMakeDiJetIniHist(std::string fList = "", const char* outFileName = "raw_rawOverGen", Bool_t isPbPb = false)
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
    while(!inFile.eof()){
      inFile >> buffer;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }

  if(!isPbPb){
    for(Int_t iter = 0; iter < 5; iter++){
      std::cout << ptHatWeights_PYTH[iter] << std::endl;
    }
  }
  else{
    for(Int_t iter = 0; iter < 4; iter++){
      std::cout << ptHatWeights_PYTHHYD[iter] << std::endl;
    }
  }

  if(!isPbPb){
    makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "ak3PF", isPbPb);
    makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "ak3Calo", isPbPb);
  }

  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs3PF", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu3PF", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs3Calo", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu3Calo", isPbPb);

  if(!isPbPb){
    makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "ak4PF", isPbPb);
    makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "ak4Calo", isPbPb);
  }

  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs4PF", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu4PF", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs4Calo", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu4Calo", isPbPb);

  if(!isPbPb){
    makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "ak5PF", isPbPb);
    makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "ak5Calo", isPbPb);
  }

  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs5PF", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu5PF", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs5Calo", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu5Calo", isPbPb);

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
