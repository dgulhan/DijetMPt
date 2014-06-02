//=============================================                                                
// Author: Chris McGinn                                                                        
//                                                                                             
// DiJet Histogram Maker                                                              
//                                                                                             
//=============================================     

#include "commonUtility.h"
#include "TTree.h"
#include "TDatime.h"
#include "TFile.h"
#include "diJetFileTag.h"

TFile* inFile_p = 0;
TFile* outFile_p = 0;

Int_t leadJtCut = 120;
Int_t subLeadJtCut = 50;

void makeImbAsymmHist(TTree* anaTree_p, const char* outName, const char* gorr, Int_t setNum, const char* CNC, const char* FPT, Int_t centLow, Int_t centHi, Int_t histLow, Int_t histHi, const char* Corr = "", Bool_t montecarlo = false)
{
  inFile_p->cd();

  Int_t setCorrNum = setNum;
  if(!strcmp("Corr", Corr))
    setCorrNum = setNum + 3;

  const char* title = Form("%s%sImbAsymmProjA%s%s%s_%d%d_%s_h", gorr, algType[setNum], CNC, FPT, Corr, (Int_t)(centLow*.5), (Int_t)((centHi + 1)*.5), fileTag);

  Float_t xArr[5] = {.0, .11, .22, .33, .50};

  TH1F* imbAsymmHist_p = new TH1F("imbAsymmHist_p", "imbAsymmHist_p", 4, xArr);
  imbAsymmHist_p->GetXaxis()->SetLimits(0.00, 0.50);
  niceTH1(imbAsymmHist_p, histHi, histLow, 505, 406);

  TH1F* getHist_p;

  TString var = Form("%sAlgImbProjA%s%s[%d]", gorr, CNC, FPT, setCorrNum);

  if(!strcmp("F", FPT) && !strcmp("", CNC) && centHi == 19 && !strcmp(Corr, "Corr"))
    std::cout << var << std::endl;

  TCut setCut = makeSetCut(setNum);
  TCut centCut = makeCentCut(centLow, centHi);
  TCut etaCut = makeEtaCut(setNum, 1.6);
  TCut phiCut = makeDelPhiCut(setNum, 5*TMath::Pi()/6);

  TCut jetLCut = Form("AlgLeadJtPt[%d] > %d", setNum, leadJtCut);
  TCut jetSLCut = Form("AlgSubLeadJtPt[%d] > %d", setNum, subLeadJtCut);
  TCut pthat = "pthat > 80";

  const char* name1[4] = {"0_1(10000, -10000, 10000)", "1_2(10000, -10000, 10000)", "2_3(10000, -10000, 10000)", "3_5(10000, -10000, 10000)"};
  const char* name2[4] = {"0_1", "1_2", "2_3", "3_5"};
  Float_t asymmBins[5] = {.00, .11, .22, .33, 1.00};

  if(!strcmp("F", FPT) && !strcmp("", CNC) && centHi == 19 && !strcmp(Corr, "Corr"))
    std::cout << setCut << ", " << centCut << ", " << etaCut << ", " << phiCut << ", " << jetLCut << ", " << jetSLCut << std::endl;

  for(Int_t binIter = 0; binIter < 4; binIter++){
    TCut asymmCut = makeAsymmCut(setNum, asymmBins[binIter], asymmBins[binIter + 1]);

    if(!strcmp("F", FPT) && !strcmp("", CNC) && centHi == 19 && !strcmp(Corr, "Corr"))
      std::cout << asymmCut << std::endl;

    if(montecarlo)
      anaTree_p->Project(name1[binIter], var, setCut && centCut && etaCut && phiCut && jetLCut && jetSLCut && asymmCut, "");
    else
      anaTree_p->Project(name1[binIter], var, setCut && centCut && etaCut && phiCut && jetLCut && jetSLCut && asymmCut, "");

    getHist_p = (TH1F*)inFile_p->Get(name2[binIter]);

    imbAsymmHist_p->SetBinContent(binIter + 1, getHist_p->GetMean());
    imbAsymmHist_p->SetBinError(binIter + 1, getHist_p->GetMeanError());

    if(!strcmp("F", FPT) && !strcmp("", CNC) && centHi == 19 && !strcmp(Corr, "Corr"))
      std::cout << Corr << ", " << centHi << ", " << CNC << ", " << FPT << ", " << "Bin " << binIter << ": " << getHist_p->GetMean() << std::endl;
  }

  outFile_p = new TFile(outName, "UPDATE");
  imbAsymmHist_p->Write(title);
  outFile_p->Close();

  delete outFile_p;
  delete imbAsymmHist_p;

  return;
}


void makeDiJetHists(const char* inName, const char* outName, Bool_t montecarlo = false)
{
  TH1::SetDefaultSumw2();

  setFileTag(inName);

  inFile_p = new TFile(inName, "READ");
  TTree* anaTree_p = (TTree*)inFile_p->Get("jetTreeAna");
  anaTree_p->AddFriend("trackTreeAna");

  Int_t jetAlgMax = 2;

  if(montecarlo){
    anaTree_p->AddFriend("genTreeAna");
    jetAlgMax = 3;
  }

  const char* corr[2] = {"", "Corr"};
  const char* CNC[3] = {"", "C", "NC"};
  const char* FPT[6] = {"F", "0_1", "1_2", "2_4", "4_8", "8_100"};

  Int_t centLow[6] = {0, 20, 60, 100, 0, 60};
  Int_t centHi[6] = {19, 59, 99, 199, 59, 199};

  for(Int_t algIter = 0; algIter < jetAlgMax; algIter++){
    std::cout << "Algorithm: " << algType[algIter] << std::endl;

    for(Int_t corrIter = 0; corrIter < 2; corrIter++){
      for(Int_t CNCIter = 0; CNCIter < 3; CNCIter++){
	for(Int_t centIter = 0; centIter < 6; centIter++){
	  for(Int_t FPTIter = 0; FPTIter < 6; FPTIter++){

	    if((CNCIter == 0 && centIter < 4) || (CNCIter != 0 && centIter >= 4)){
	      makeImbAsymmHist(anaTree_p, outName, "r", algIter, CNC[CNCIter], FPT[FPTIter], centLow[centIter], centHi[centIter], -59.999, 59.999, corr[corrIter], montecarlo);

	      if(montecarlo && corrIter == 0)
		makeImbAsymmHist(anaTree_p, outName, "g", algIter, CNC[CNCIter], FPT[FPTIter], centLow[centIter], centHi[centIter], -59.999, 59.999, corr[corrIter], montecarlo);
	    }

	  }
	}
      }
    }

  }

  inFile_p->Close();
  delete inFile_p;

  return;
}
