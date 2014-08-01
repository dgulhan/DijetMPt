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

#include <vector>
#include "cfmVectFunc.h"

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
Float_t refeta_[maxJets];

Int_t ngen_;
Float_t genpt_[maxJets];
Float_t geneta_[maxJets];
Int_t genmatchindex_[maxJets];

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
  getChain_p[0]->SetBranchStatus("refeta", 1);

  getChain_p[0]->SetBranchStatus("ngen", 1);
  getChain_p[0]->SetBranchStatus("genpt", 1);
  getChain_p[0]->SetBranchStatus("geneta", 1);
  getChain_p[0]->SetBranchStatus("genmatchindex", 1);

  if(isPbPb) getChain_p[0]->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection_);

  getChain_p[0]->SetBranchAddress("vz", &vz_);
  getChain_p[0]->SetBranchAddress("pthat", &pthat_);
  getChain_p[0]->SetBranchAddress("hiBin", &hiBin_);

  getChain_p[0]->SetBranchAddress("nref", &nref_);
  getChain_p[0]->SetBranchAddress("jtpt", &jtpt_);
  getChain_p[0]->SetBranchAddress("jteta", &jteta_);
  getChain_p[0]->SetBranchAddress("rawpt", &rawpt_);
  getChain_p[0]->SetBranchAddress("refpt", &refpt_);
  getChain_p[0]->SetBranchAddress("refeta", &refeta_);

  getChain_p[0]->SetBranchAddress("ngen", &ngen_);
  getChain_p[0]->SetBranchAddress("genpt", &genpt_);
  getChain_p[0]->SetBranchAddress("geneta", &geneta_);
  getChain_p[0]->SetBranchAddress("genmatchindex", &genmatchindex_);

  return;
}


void GetChain(std::vector<std::string> inList, const std::string alg, Bool_t isPbPb)
{
  getChain_p[0] = new TChain(Form("%sJetAnalyzer/t", alg.c_str()));
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

    if(iter > 100){
      std::cout << iter << ", " << val << ", " << pos << ", " << arr[pos] << std::endl;
    }
  }

  return pos;
}


void makeDiJetIniHist(std::vector<std::string> inList, const std::string outName, const std::string alg = "akVs3Calo", Bool_t isPbPb = false)
{
  TH1::SetDefaultSumw2();

  std::cout << "Init: " << alg << std::endl;

  GetChain(inList, alg, isPbPb);

  const Int_t nPtBins = 25;
  const Float_t ptLower = 20.;
  const Float_t ptHigher = 700.;

  const Int_t nEtaBins = 25;
  const Float_t etaLower = -2.0;
  const Float_t etaHigher = 2.0;

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

  Float_t ptBins[nPtBins+1];
  Float_t etaBins[nEtaBins+1];

  getLogBins(ptLower, ptHigher, nPtBins, ptBins);
  getEtaBins(etaLower, etaHigher, nEtaBins, etaBins);

  TH1F* pthatNoWeight_p = new TH1F("pthatNoWeight_h", "pthatNoWeight_h", 50, 0, 500);
  TH1F* pthatWeight_p = new TH1F("pthatWeight_h", "pthatWeight_h", 50, 0, 500);
  TH1F* temp_jetOverGen_Pt_p[nCentBins][nPtBins];
  TH1F* temp_jetOverGen_Eta_p[nCentBins][nEtaBins];
  TH1F* temp_genEff_Pt_p[nCentBins][nPtBins];
  TH1F* temp_genEff_Eta_p[nCentBins][nPtBins];

  Float_t effBinning[3] = {-.5, .5, 1.5};

  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    for(Int_t iter = 0; iter < nPtBins; iter++){
      temp_jetOverGen_Pt_p[centIter][iter] = new TH1F(Form("%s_tempJetHist_Pt_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), Form("%s_tempJetHist_Pt_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), 200, 0, 2);
      temp_jetOverGen_Eta_p[centIter][iter] = new TH1F(Form("%s_tempJetHist_Eta_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), Form("%s_tempJetHist_Eta_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), 200, 0, 2);

      temp_genEff_Pt_p[centIter][iter] = new TH1F(Form("%s_tempGenEffHist_Pt_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), Form("%s_tempGenEffHist_Pt_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), 2, effBinning);
      temp_genEff_Eta_p[centIter][iter] = new TH1F(Form("%s_tempGenEffHist_Eta_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), Form("%s_tempGenEffHist_Eta_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), 2, effBinning);

    }
  }

  TH1F* jetOverGenMean_Pt_p[nCentBins];
  TH1F* jetOverGenMean_Eta_p[nCentBins];
  TH1F* jetOverGenRes_Pt_p[nCentBins];
  TH1F* jetOverGenRes_Eta_p[nCentBins];
  TH1F* genEff_Pt_p[nCentBins];
  TH1F* genEff_Eta_p[nCentBins];

  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    jetOverGenMean_Pt_p[centIter] = new TH1F(Form("%s_jetOverGenMean_Pt_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_jetOverGenMean_Pt_%s_h", alg.c_str(), centString[centIter].c_str()), nPtBins, ptBins);
    handsomeTH1(jetOverGenMean_Pt_p[centIter]);
    jetOverGenMean_Pt_p[centIter]->SetXTitle("p_{T}^{gen}");
    jetOverGenMean_Pt_p[centIter]->GetXaxis()->SetTitleOffset(.75);
    jetOverGenMean_Pt_p[centIter]->SetYTitle("<p_{T}^{jet}/p_{T}^{gen}>");

    jetOverGenMean_Eta_p[centIter] = new TH1F(Form("%s_jetOverGenMean_Eta_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_jetOverGenMean_Eta_%s_h", alg.c_str(), centString[centIter].c_str()), nEtaBins, etaBins);
    handsomeTH1(jetOverGenMean_Eta_p[centIter]);
    jetOverGenMean_Eta_p[centIter]->SetXTitle("#eta_{gen}");
    jetOverGenMean_Eta_p[centIter]->GetXaxis()->SetTitleOffset(.75);
    jetOverGenMean_Eta_p[centIter]->SetYTitle("<p_{T}^{jet}/p_{T}^{gen}>");

    jetOverGenRes_Pt_p[centIter] = new TH1F(Form("%s_jetOverGenRes_Pt_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_jetOverGenRes_Pt_%s_h", alg.c_str(), centString[centIter].c_str()), nPtBins, ptBins);
    handsomeTH1(jetOverGenRes_Pt_p[centIter]);
    jetOverGenRes_Pt_p[centIter]->SetXTitle("p_{T}^{gen}");
    jetOverGenRes_Pt_p[centIter]->GetXaxis()->SetTitleOffset(.75);
    jetOverGenRes_Pt_p[centIter]->SetYTitle("#sigma_{reco/gen}");

    jetOverGenRes_Eta_p[centIter] = new TH1F(Form("%s_jetOverGenRes_Eta_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_jetOverGenRes_Eta_%s_h", alg.c_str(), centString[centIter].c_str()), nEtaBins, etaBins);
    handsomeTH1(jetOverGenRes_Eta_p[centIter]);
    jetOverGenRes_Eta_p[centIter]->SetXTitle("#eta_{gen}");
    jetOverGenRes_Eta_p[centIter]->GetXaxis()->SetTitleOffset(.75);
    jetOverGenRes_Eta_p[centIter]->SetYTitle("#sigma_{reco/gen}");

    genEff_Pt_p[centIter] = new TH1F(Form("%s_genEff_Pt_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_genEff_Pt_%s_h", alg.c_str(), centString[centIter].c_str()), nPtBins, ptBins);
    handsomeTH1(genEff_Pt_p[centIter]);
    genEff_Pt_p[centIter]->SetXTitle("p_{T}^{gen}");
    genEff_Pt_p[centIter]->GetXaxis()->SetTitleOffset(.75);
    genEff_Pt_p[centIter]->SetYTitle("Efficiency");

    genEff_Eta_p[centIter] = new TH1F(Form("%s_genEff_Eta_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_genEff_Eta_%s_h", alg.c_str(), centString[centIter].c_str()), nEtaBins, etaBins);
    handsomeTH1(genEff_Eta_p[centIter]);
    genEff_Eta_p[centIter]->SetXTitle("#eta_{gen}");
    genEff_Eta_p[centIter]->GetXaxis()->SetTitleOffset(.75);
    genEff_Eta_p[centIter]->SetYTitle("Efficiency");
  }

  for(Int_t jEntry = 0; jEntry <  getChain_p[0]->GetEntries(); jEntry++){
    getChain_p[0]->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(!pcollisionEventSelection_ && isPbPb) continue;

    if(TMath::Abs(vz_) > 15) continue;

    Float_t hatWeight = getHatWeight(pthat_, isPbPb);
    Int_t centPos = 0;
    if(isPbPb) centPos = getCentPos(hiBin_);

    pthatNoWeight_p->Fill(pthat_);
    pthatWeight_p->Fill(pthat_, hatWeight);

    for(Int_t jtIter = 0; jtIter < nref_; jtIter++){
      if(TMath::Abs(refeta_[jtIter]) >= 2 || refpt_[jtIter] < jtPtCut)	continue;

      if(refpt_[jtIter] > pthat_ + 30) continue;

      Int_t ptPos = posSearch(refpt_[jtIter], nPtBins+1, ptBins);
      temp_jetOverGen_Pt_p[centPos][ptPos]->Fill(jtpt_[jtIter]/refpt_[jtIter], hatWeight);

      if(refpt_[jtIter] > 50){
	Int_t etaPos = posSearch(refeta_[jtIter], nEtaBins+1, etaBins);
	temp_jetOverGen_Eta_p[centPos][etaPos]->Fill(jtpt_[jtIter]/refpt_[jtIter], hatWeight);
      }
    }

    for(Int_t genIter = 0; genIter < ngen_; genIter++){
      if(genpt_[genIter] < jtPtCut) break;

      if(TMath::Abs(geneta_[genIter]) >= 2) continue;

      Int_t ptPos = posSearch(genpt_[genIter], nPtBins+1, ptBins);
      if(genmatchindex_[genIter] >= 0) temp_genEff_Pt_p[centPos][ptPos]->Fill(1.0, hatWeight);
      else temp_genEff_Pt_p[centPos][ptPos]->Fill(0.0, hatWeight);
      
      if(genpt_[genIter] > 50){
	Int_t etaPos = posSearch(geneta_[genIter], nEtaBins+1, etaBins);
	if(genmatchindex_[genIter] >= 0) temp_genEff_Eta_p[centPos][etaPos]->Fill(1.0, hatWeight);
	else temp_genEff_Eta_p[centPos][etaPos]->Fill(0.0, hatWeight);
      }
    }
  }


  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    for(Int_t rawIter = 0; rawIter < nPtBins; rawIter++){

      if(temp_jetOverGen_Pt_p[centIter][rawIter]->GetEntries() != 0){
	temp_jetOverGen_Pt_p[centIter][rawIter]->Fit("gaus", "Q WL");
	jetOverGenMean_Pt_p[centIter]->SetBinContent(rawIter+1, temp_jetOverGen_Pt_p[centIter][rawIter]->GetFunction("gaus")->GetParameter(1));
	jetOverGenRes_Pt_p[centIter]->SetBinContent(rawIter+1, temp_jetOverGen_Pt_p[centIter][rawIter]->GetFunction("gaus")->GetParameter(2));
	
	if(temp_jetOverGen_Pt_p[centIter][rawIter]->GetEntries() != 1){
	  jetOverGenMean_Pt_p[centIter]->SetBinError(rawIter+1, temp_jetOverGen_Pt_p[centIter][rawIter]->GetFunction("gaus")->GetParError(1));
	  jetOverGenRes_Pt_p[centIter]->SetBinError(rawIter+1, temp_jetOverGen_Pt_p[centIter][rawIter]->GetFunction("gaus")->GetParError(2));
	  
	  if(jetOverGenMean_Pt_p[centIter]->GetBinError(rawIter+1) > .05){
	    std::cout << rawIter << ", " << ptBins[rawIter] << ", " << ptBins[rawIter+1] << ", " << jetOverGenMean_Pt_p[centIter]->GetBinContent(rawIter + 1) << ", " << jetOverGenMean_Pt_p[centIter]->GetBinError(rawIter+1) << std::endl;
	  }
	}
      }

      if(temp_jetOverGen_Eta_p[centIter][rawIter]->GetEntries() != 0){
        temp_jetOverGen_Eta_p[centIter][rawIter]->Fit("gaus", "Q WL");
        jetOverGenMean_Eta_p[centIter]->SetBinContent(rawIter+1, temp_jetOverGen_Eta_p[centIter][rawIter]->GetFunction("gaus")->GetParameter(1));
        jetOverGenRes_Eta_p[centIter]->SetBinContent(rawIter+1, temp_jetOverGen_Eta_p[centIter][rawIter]->GetFunction("gaus")->GetParameter(2));

        if(temp_jetOverGen_Eta_p[centIter][rawIter]->GetEntries() != 1){
          jetOverGenMean_Eta_p[centIter]->SetBinError(rawIter+1, temp_jetOverGen_Eta_p[centIter][rawIter]->GetFunction("gaus")->GetParError(1));
          jetOverGenRes_Eta_p[centIter]->SetBinError(rawIter+1, temp_jetOverGen_Eta_p[centIter][rawIter]->GetFunction("gaus")->GetParError(2));


          if(jetOverGenMean_Eta_p[centIter]->GetBinError(rawIter+1) > .05){
	    std::cout << rawIter << ", " << ptBins[rawIter] << ", " << ptBins[rawIter+1] << ", " << jetOverGenMean_Eta_p[centIter]->GetBinContent(rawIter + 1) << ", " << jetOverGenMean_Eta_p[centIter]->GetBinError(rawIter+1) << std::endl;
          }
        }
      }

      if(temp_genEff_Pt_p[centIter][rawIter]->GetEntries() != 0){
	genEff_Pt_p[centIter]->SetBinContent(rawIter+1, temp_genEff_Pt_p[centIter][rawIter]->GetMean());

	if(temp_genEff_Pt_p[centIter][rawIter]->GetEntries() != 1)
	  genEff_Pt_p[centIter]->SetBinError(rawIter+1, temp_genEff_Pt_p[centIter][rawIter]->GetMeanError());
      }

      if(temp_genEff_Eta_p[centIter][rawIter]->GetEntries() != 0){
	genEff_Eta_p[centIter]->SetBinContent(rawIter+1, temp_genEff_Eta_p[centIter][rawIter]->GetMean());

	if(temp_genEff_Eta_p[centIter][rawIter]->GetEntries() != 1)
	  genEff_Eta_p[centIter]->SetBinError(rawIter+1, temp_genEff_Eta_p[centIter][rawIter]->GetMeanError());
      }

    }
  }


  TFile* out = new TFile(outName.c_str(), "UPDATE");
  std::cout << outName << std::endl;
  if(!strcmp(alg.c_str(), "akVs3PF")){
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
    jetOverGenMean_Pt_p[centIter]->Write();
    jetOverGenMean_Eta_p[centIter]->Write();
    jetOverGenRes_Pt_p[centIter]->Write();
    jetOverGenRes_Eta_p[centIter]->Write();
    genEff_Pt_p[centIter]->Write();
    genEff_Eta_p[centIter]->Write();
  }
  out->Close();
  delete out;

  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    delete jetOverGenMean_Pt_p[centIter];
    jetOverGenMean_Pt_p[centIter] = 0;

    delete jetOverGenMean_Eta_p[centIter];
    jetOverGenMean_Eta_p[centIter] = 0;

    delete jetOverGenRes_Pt_p[centIter];
    jetOverGenRes_Pt_p[centIter] = 0;

    delete jetOverGenRes_Eta_p[centIter];
    jetOverGenRes_Eta_p[centIter] = 0;

    delete genEff_Pt_p[centIter];
    genEff_Pt_p[centIter] = 0;

    delete genEff_Eta_p[centIter];
    genEff_Eta_p[centIter] = 0;
  }

  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    for(Int_t rawIter = 0; rawIter < nPtBins; rawIter++){
      delete temp_jetOverGen_Pt_p[centIter][rawIter];
      temp_jetOverGen_Pt_p[centIter][rawIter] = 0;   

      delete temp_jetOverGen_Eta_p[centIter][rawIter];
      temp_jetOverGen_Eta_p[centIter][rawIter] = 0;   

      delete temp_genEff_Pt_p[centIter][rawIter];
      temp_genEff_Pt_p[centIter][rawIter] = 0;

      delete temp_genEff_Eta_p[centIter][rawIter];
      temp_genEff_Eta_p[centIter][rawIter] = 0;
    }
  }

  delete pthatWeight_p;
  pthatWeight_p = 0;

  delete pthatNoWeight_p;
  pthatNoWeight_p = 0;

  CleanChain();

  return;
}


void drawLine(const std::string MeanRes, const std::string PtEta)
{
  if(!strcmp(MeanRes.c_str(), "Mean")){
    TLine* oneLine_p;
    if(!strcmp(PtEta.c_str(), "Pt")) oneLine_p = new TLine(20.0, 1.0, 700.0, 1.0);
    else oneLine_p = new TLine(-2.0, 1.0, 2.0, 1.0);

    oneLine_p->SetLineColor(1);
    oneLine_p->SetLineStyle(2);
    oneLine_p->Draw("SAME");

    TLine* oneUpLine_p;
    if(!strcmp(PtEta.c_str(), "Pt")) oneUpLine_p = new TLine(20.0, 1.01, 700.0, 1.01);
    else oneUpLine_p = new TLine(-2.0, 1.01, 2.0, 1.01);

    oneUpLine_p->SetLineColor(1);
    oneUpLine_p->SetLineStyle(2);
    oneUpLine_p->Draw("SAME");

    TLine* oneDownLine_p;
    if(!strcmp(PtEta.c_str(), "Pt")) oneDownLine_p = new TLine(20.0, 0.99, 700.0, 0.99);
    else oneDownLine_p = new TLine(-2.0, 0.99, 2.0, 0.99);

    oneDownLine_p->SetLineColor(1);
    oneDownLine_p->SetLineStyle(2);
    oneDownLine_p->Draw("SAME");
  }    

  TLine* fiftyLine_p;
  if(!strcmp(MeanRes.c_str(), "Mean")) fiftyLine_p = new TLine(50.0, 0.95, 50.0, 1.05);
  else fiftyLine_p = new TLine(50.0, 0.00, 50.0, 0.50);

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
  //  drawLine();

  return;
}


void plotDiJetIniHistRatPFCalo(const char* histFileName, const char* VsPu = "Vs")
{
  TH1::SetDefaultSumw2();

  TFile *f = new TFile(histFileName, "READ");
  TH1F* getHist_p = (TH1F*)f->Get(Form("ak3%s_PFOverCalo_h", VsPu));
  getHist_p->DrawCopy();
  //  drawLine();

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


Float_t getHistMax(const std::string MeanRes)
{
  if(!strcmp(MeanRes.c_str(), "Mean")) return 1.04999;
  else return 0.34999;
}


Float_t getHistMin(const std::string MeanRes)
{
  if(!strcmp(MeanRes.c_str(), "Mean")) return 0.95001;
  else return 0.00001;
}


void drawGenHist(TH1F* inHist_p, const std::string MeanRes, const std::string PtEta)
{
  if(!strcmp(PtEta.c_str(), "Pt")) gPad->SetLogx();
  inHist_p->Draw();
  drawLine(MeanRes, PtEta);
  return;
}


void plotDiJetIniHist_PYTH(const std::string histFileName, const std::string VsPu, const std::string PFCalo, const std::string MeanRes, const std::string PtEta)
{
  TH1::SetDefaultSumw2();
  TFile *f = new TFile(histFileName.c_str(), "UPDATE");

  std::cout << Form("ak%s3%s_jetOverGen%s_%s_PP_h", VsPu.c_str(), PFCalo.c_str(), MeanRes.c_str(), PtEta.c_str()) << std::endl;

  Float_t max = getHistMax(MeanRes);
  Float_t min = getHistMin(MeanRes);

  TH1F* getHist2_p = (TH1F*)f->Get(Form("ak%s2%s_jetOverGen%s_%s_PP_h", VsPu.c_str(), PFCalo.c_str(), MeanRes.c_str(), PtEta.c_str()));
  getHist2_p->SetMaximum(max);
  getHist2_p->SetMinimum(min);

  TH1F* getHist3_p = (TH1F*)f->Get(Form("ak%s3%s_jetOverGen%s_%s_PP_h", VsPu.c_str(), PFCalo.c_str(), MeanRes.c_str(), PtEta.c_str()));
  getHist3_p->SetMaximum(max);
  getHist3_p->SetMinimum(min);

  TH1F* getHist4_p = (TH1F*)f->Get(Form("ak%s4%s_jetOverGen%s_%s_PP_h", VsPu.c_str(), PFCalo.c_str(), MeanRes.c_str(), PtEta.c_str()));
  getHist4_p->SetMaximum(max);
  getHist4_p->SetMinimum(min);

  TH1F* getHist5_p = (TH1F*)f->Get(Form("ak%s5%s_jetOverGen%s_%s_PP_h", VsPu.c_str(), PFCalo.c_str(), MeanRes.c_str(), PtEta.c_str()));
  getHist5_p->SetMaximum(max);
  getHist5_p->SetMinimum(min);

  TCanvas* plotCanv_p = new TCanvas(Form("ak%s%s_jetOverGen%s_%s_PP_c", VsPu.c_str(), PFCalo.c_str(), MeanRes.c_str(), PtEta.c_str()), Form("ak%s%s_jetOverGen%s_%s_PP_c", VsPu.c_str(), PFCalo.c_str(), MeanRes.c_str(), PtEta.c_str()), 4*300, 1*350);
  plotCanv_p->Divide(4, 1, 0.0, 0.0);

  plotInitHist(getHist2_p);
  plotInitHist(getHist3_p);
  plotInitHist(getHist4_p);
  plotInitHist(getHist5_p);

  plotCanv_p->cd(1);
  drawGenHist(getHist2_p, MeanRes, PtEta);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(23);

  label_p->DrawLatex(.5, .9, Form("ak%s2%s", VsPu.c_str(), PFCalo.c_str()));

  if(!strcmp(PtEta.c_str(), "Pt")) label_p->DrawLatex(.5, .80, Form("p_{T}^{gen} > 20 GeV/c"));
  else label_p->DrawLatex(.5, .80, Form("p_{T}^{gen} > 50 GeV/c"));

  plotCanv_p->cd(2);
  drawGenHist(getHist3_p, MeanRes, PtEta);
  label_p->DrawLatex(.5, .9, Form("ak%s3%s", VsPu.c_str(), PFCalo.c_str()));
  label_p->DrawLatex(.5, .80, Form("|#eta| < 2.0"));

  plotCanv_p->cd(3);
  drawGenHist(getHist4_p, MeanRes, PtEta);
  label_p->DrawLatex(.6, .9, Form("ak%s4%s", VsPu.c_str(), PFCalo.c_str()));
  label_p->DrawLatex(.6, .80, "PYTHIA");

  plotCanv_p->cd(4);
  drawGenHist(getHist5_p, MeanRes, PtEta);
  label_p->DrawLatex(.6, .9, Form("ak%s5%s", VsPu.c_str(), PFCalo.c_str()));
  label_p->DrawLatex(.6, .80, "PYTHIA");

  plotCanv_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(plotCanv_p, Form("pdfDir/ak%s%s_jetOverGen%s_%s_PP", VsPu.c_str(), PFCalo.c_str(), MeanRes.c_str(), PtEta.c_str()), "pdf");

  delete label_p;
  delete plotCanv_p;
  f->Close();
  delete f;
}

void plotDiJetIniHist_PYTHHYD(const std::string histFileName, const std::string alg, const std::string MeanRes, const std::string PtEta)
{
  TH1::SetDefaultSumw2();
  TFile *f = new TFile(histFileName.c_str(), "UPDATE");

  Float_t max = getHistMax(MeanRes);
  Float_t min = getHistMin(MeanRes);

  TH1F* getHist_p[8];

  for(Int_t iter = 0; iter < 8; iter++){
    getHist_p[iter] = (TH1F*)f->Get(Form("%s_jetOverGen%s_%s_%d%d_h", alg.c_str(), MeanRes.c_str(), PtEta.c_str(), centCutArray[iter]/2, centCutArray[iter+1]/2));
    getHist_p[iter]->SetMaximum(max);
    getHist_p[iter]->SetMinimum(min);
  }

  TCanvas* plotCanv_p = new TCanvas(Form("%s_jetOverGen%s_%s_PbPb_c", alg.c_str(), MeanRes.c_str(), PtEta.c_str()), Form("%s_jetOverGen%s_%s_PbPb_c", alg.c_str(), MeanRes.c_str(), PtEta.c_str()), 4*300, 2*350);
  plotCanv_p->Divide(4, 2, 0.0, 0.0);

  for(Int_t iter = 0; iter < 8; iter++){
    plotInitHist(getHist_p[iter]);
  }

  plotCanv_p->cd(1);
  drawGenHist(getHist_p[7], MeanRes, PtEta);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(23);

  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[7]/2, centCutArray[8]/2));
  if(!strcmp(PtEta.c_str(), "Pt")) label_p->DrawLatex(.5, .70, Form("p_{T}^{gen} > 20 GeV/c"));
  else label_p->DrawLatex(.5, .70, Form("p_{T}^{gen} > 50 GeV/c"));

  plotCanv_p->cd(2);
  drawGenHist(getHist_p[6], MeanRes, PtEta);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[6]/2, centCutArray[7]/2));
  label_p->DrawLatex(.5, .70, Form("|#eta| < 2.0"));

  plotCanv_p->cd(3);
  drawGenHist(getHist_p[5], MeanRes, PtEta);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[5]/2, centCutArray[6]/2));
  label_p->DrawLatex(.5, .70, "PYT+HYD");

  plotCanv_p->cd(4);
  drawGenHist(getHist_p[4], MeanRes, PtEta);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[4]/2, centCutArray[5]/2));

  plotCanv_p->cd(5);
  drawGenHist(getHist_p[3], MeanRes, PtEta);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[3]/2, centCutArray[4]/2));

  plotCanv_p->cd(6);
  drawGenHist(getHist_p[2], MeanRes, PtEta);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[2]/2, centCutArray[3]/2));

  plotCanv_p->cd(7);
  drawGenHist(getHist_p[1], MeanRes, PtEta);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[1]/2, centCutArray[2]/2));

  plotCanv_p->cd(8);
  drawGenHist(getHist_p[0], MeanRes, PtEta);
  label_p->DrawLatex(.5, .90, Form("%s", alg.c_str()));
  label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[0]/2, centCutArray[1]/2));


  plotCanv_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(plotCanv_p, Form("pdfDir/%s_jetOverGen%s_%s", alg.c_str(), MeanRes.c_str(), PtEta.c_str()), "pdf");

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
  const std::string PtEta[2] = {"Pt", "Eta"};
  const std::string MeanRes[2] = {"Mean", "Res"};

  for(Int_t ptEtaIter = 0; ptEtaIter < 2; ptEtaIter++){
    for(Int_t meanResIter = 0; meanResIter < 2; meanResIter++){
      plotDiJetIniHist_PYTH(histFileName, "", "PF", MeanRes[meanResIter], PtEta[ptEtaIter]);
      plotDiJetIniHist_PYTH(histFileName, "", "Calo", MeanRes[meanResIter], PtEta[ptEtaIter]);
      
      plotDiJetIniHist_PYTH(histFileName, "Vs", "PF", MeanRes[meanResIter], PtEta[ptEtaIter]);
      plotDiJetIniHist_PYTH(histFileName, "Vs", "Calo", MeanRes[meanResIter], PtEta[ptEtaIter]);
      
      plotDiJetIniHist_PYTH(histFileName, "Pu", "PF", MeanRes[meanResIter], PtEta[ptEtaIter]);
      plotDiJetIniHist_PYTH(histFileName, "Pu", "Calo", MeanRes[meanResIter], PtEta[ptEtaIter]);
    }
  }

  return;
}


void runPlotDiJetIniHist_PYTHHYD(const std::string histFileName)
{
  const std::string PtEta[2] = {"Pt", "Eta"};
  const std::string MeanRes[2] = {"Mean", "Res"};

  for(Int_t ptEtaIter = 0; ptEtaIter < 2; ptEtaIter++){
    for(Int_t meanResIter = 0; meanResIter < 2; meanResIter++){
      plotDiJetIniHist_PYTHHYD(histFileName, "akVs3PF", MeanRes[meanResIter], PtEta[ptEtaIter]);
      plotDiJetIniHist_PYTHHYD(histFileName, "akVs3Calo", MeanRes[meanResIter], PtEta[ptEtaIter]);      
      plotDiJetIniHist_PYTHHYD(histFileName, "akPu3PF", MeanRes[meanResIter], PtEta[ptEtaIter]);
      plotDiJetIniHist_PYTHHYD(histFileName, "akPu3Calo", MeanRes[meanResIter], PtEta[ptEtaIter]);

      plotDiJetIniHist_PYTHHYD(histFileName, "akVs3PF", MeanRes[meanResIter], PtEta[ptEtaIter]);
      plotDiJetIniHist_PYTHHYD(histFileName, "akVs3Calo", MeanRes[meanResIter], PtEta[ptEtaIter]);      
      plotDiJetIniHist_PYTHHYD(histFileName, "akPu3PF", MeanRes[meanResIter], PtEta[ptEtaIter]);
      plotDiJetIniHist_PYTHHYD(histFileName, "akPu3Calo", MeanRes[meanResIter], PtEta[ptEtaIter]);

      plotDiJetIniHist_PYTHHYD(histFileName, "akVs4PF", MeanRes[meanResIter], PtEta[ptEtaIter]);
      plotDiJetIniHist_PYTHHYD(histFileName, "akVs4Calo", MeanRes[meanResIter], PtEta[ptEtaIter]);      
      plotDiJetIniHist_PYTHHYD(histFileName, "akPu4PF", MeanRes[meanResIter], PtEta[ptEtaIter]);
      plotDiJetIniHist_PYTHHYD(histFileName, "akPu4Calo", MeanRes[meanResIter], PtEta[ptEtaIter]);
      
      plotDiJetIniHist_PYTHHYD(histFileName, "akVs5PF", MeanRes[meanResIter], PtEta[ptEtaIter]);
      plotDiJetIniHist_PYTHHYD(histFileName, "akVs5Calo", MeanRes[meanResIter], PtEta[ptEtaIter]);      
      plotDiJetIniHist_PYTHHYD(histFileName, "akPu5PF", MeanRes[meanResIter], PtEta[ptEtaIter]);
      plotDiJetIniHist_PYTHHYD(histFileName, "akPu5Calo", MeanRes[meanResIter], PtEta[ptEtaIter]);
    }    
  }

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
    makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "ak2PF", isPbPb);
    makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "ak2Calo", isPbPb);
  }
  
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs2PF", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu2PF", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akVs2Calo", isPbPb);
  makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), "akPu2Calo", isPbPb);

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
