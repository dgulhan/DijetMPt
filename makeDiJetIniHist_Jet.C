#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TCut.h"
#include "TMath.h"

#include "TCanvas.h"
#include "TLine.h"

#include <iostream>

void makeDiJetIniHist(const char* inName, const char* outName, const char* alg = "akVs3Calo")
{
  TH1::SetDefaultSumw2();

  TFile* f = new TFile(inName, "READ");
  TTree* getTree_p = (TTree*)f->Get(Form("%sJetAnalyzer/t", alg));
  getTree_p->AddFriend("skimanalysis/HltTree");
  getTree_p->AddFriend("hiEvtAnalyzer/HiTree");

  const Int_t nBins = 100;
  Float_t bins[nBins+1];
  bins[0] = 0;
  for(Int_t iter = 1; iter < nBins+1; iter++){
    bins[iter] = bins[iter - 1] + 7;
  }

  TH1F* rawOverGenHist_p = new TH1F(Form("%s_rawOverGen_h", alg), Form("%s_rawOverGen_h", alg), nBins, bins);
  TH1F* getHist_p = 0;
  TCut overallCut = "pcollisionEventSelection && TMath::Abs(jteta) < 2.0 && TMath::Abs(vz) < 15.0";

  for(Int_t iter = 0; iter < nBins; iter++){
    TCut rawCut = Form("rawpt > %f && rawpt < %f", bins[iter], bins[iter+1]);
    getTree_p->Project(Form("getHist%d_h", iter), "rawpt/genpt", overallCut && rawCut);
    getHist_p = (TH1F*)f->Get(Form("getHist%d_h", iter));

    rawOverGenHist_p->SetBinContent(iter+1, getHist_p->GetMean());
    rawOverGenHist_p->SetBinError(iter+1, getHist_p->GetMeanError());

    getHist_p = 0;
  }

  rawOverGenHist_p->SetXTitle("p_{T}^{raw}");
  rawOverGenHist_p->GetXaxis()->SetTitleOffset(1.0);
  rawOverGenHist_p->GetXaxis()->CenterTitle();
  rawOverGenHist_p->SetYTitle("p_{T}^{raw}/p_{T}^{gen}");
  rawOverGenHist_p->GetYaxis()->CenterTitle();

  TFile* out = new TFile(outName, "UPDATE");
  rawOverGenHist_p->Write();
  out->Close();
  delete out;

  delete rawOverGenHist_p;
  rawOverGenHist_p = 0;
  getTree_p = 0;

  f->Close();
  delete f;
}


void makeDiJetIniHistRatVsPu(const char* histFileName, const char* PFCalo = "PF")
{
  TH1::SetDefaultSumw2();

  TFile* f = new TFile(histFileName, "UPDATE");
  TH1F* getVs_p = (TH1F*)f->Get(Form("akVs3%s_rawOverGen_h", PFCalo));
  TH1F* getPu_p = (TH1F*)f->Get(Form("akPu3%s_rawOverGen_h", PFCalo));

  getVs_p->Divide(getPu_p);

  getVs_p->SetYTitle("Vs/Pu p_{T}^{rat}");
  getVs_p->Write(Form("ak3%s_VsOverPu_h", PFCalo));

  getPu_p = 0;
  getVs_p = 0;
  f->Close();
  delete f;
}


void makeDiJetIniHistRatPFCalo(const char* histFileName, const char* VsPu = "Vs")
{
  TH1::SetDefaultSumw2();

  TFile* f = new TFile(histFileName, "UPDATE");
  TH1F* getPF_p = (TH1F*)f->Get(Form("ak%s3PF_rawOverGen_h", VsPu));
  TH1F* getCalo_p = (TH1F*)f->Get(Form("ak%s3Calo_rawOverGen_h", VsPu));

  getPF_p->Divide(getCalo_p);

  getPF_p->SetYTitle("PF/Calo p_{T}^{rat}");
  getPF_p->Write(Form("ak3%s_PFOverCalo_h", VsPu));

  getCalo_p = 0;
  getPF_p = 0;
  f->Close();
  delete f;
}


void plotDiJetIniHistRatVsPu(const char* histFileName, const char* PFCalo = "PF")
{
  TFile *f = new TFile(histFileName, "READ");
  TH1F* getHist_p = (TH1F*)f->Get(Form("ak3%s_VsOverPu_h", PFCalo));
  TCanvas* plotCanv_p = new TCanvas("plotCanv", "plotCanv", 1);
  getHist_p->DrawCopy();
  TLine* zeroLine_p = new TLine(0.0, 1.0, 700.0, 1.0);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw("SAME");
}


void plotDiJetIniHistRatPFCalo(const char* histFileName, const char* VsPu = "Vs")
{
  TFile *f = new TFile(histFileName, "READ");
  TH1F* getHist_p = (TH1F*)f->Get(Form("ak3%s_PFOverCalo_h", VsPu));
  TCanvas* plotCanv_p = new TCanvas("plotCanv", "plotCanv", 1);
  getHist_p->DrawCopy();
  TLine* zeroLine_p = new TLine(0.0, 1.0, 700.0, 1.0);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(2);
  zeroLine_p->Draw("SAME");
}


void runMakeDiJetIniHist(const char* inName, const char* outName = "raw_rawOverGen_hat80.root")
{
  TFile* f = new TFile(outName, "RECREATE");
  f->Close();
  delete f;
  f = 0;

  makeDiJetIniHist(inName, outName, "akVs3PF");
  makeDiJetIniHist(inName, outName, "akPu3PF");
  makeDiJetIniHist(inName, outName, "akVs3Calo");
  makeDiJetIniHist(inName, outName, "akPu3Calo");

  makeDiJetIniHistRatPFCalo(outName, "Vs");
  makeDiJetIniHistRatPFCalo(outName, "Pu");
  
  makeDiJetIniHistRatVsPu(outName, "PF");
  makeDiJetIniHistRatVsPu(outName, "Calo");
}
