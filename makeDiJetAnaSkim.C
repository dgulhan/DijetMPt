//=============================================
// Author: Chris McGinn
// 
// DiJet Analysis Skim Class (MC)
//
//=============================================

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "cfmDiJetAnaSkim.h"
#include "stdlib.h"
#include <fstream>
#include "TComplex.h"

const Int_t pthatCuts_PYTH_HITrk[10] = {15, 30, 50, 80, 120, 170, 220, 280, 370, 10000000};
const Float_t pthatWeights_PYTH_HITrk[9] = {.551019, .034814, .00242254, .000304825, .0000426788, .00000492814, .000000879673, .00000017353, .0000000292439};

const Int_t pthatCuts_PYTH_PPTrk[7] = {15, 30, 50, 80, 120, 170, 1000000};
const Float_t pthatWeights_PYTH_PPTrk[6] = {.161482, .00749461, .000752396, .0000837038, .0000101988, .00000175206};

const Int_t pthatCuts_PYTH_HYD[11] = {15, 30, 50, 80, 100, 120, 170, 220, 280, 370, 10000000};
const Float_t pthatWeights_PYTH_HYD[10] = {.611066, .0374106, .00232016, .00014917, .0000822379, .0000142819, .00000296162, .00000102099, .000000522123, .000000232907};

int makeDiJetAnaSkim(std::string fList = "", sampleType sType = kHIDATA, Int_t num = 0, Bool_t justJt = false, Bool_t isHITrk = false)
{
  //Define MC or Data
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

  std::cout << sType << std::endl;
  std::cout << montecarlo << std::endl;

  std::string buffer;
  std::vector<std::string> listOfFiles;
  int nLines = 0;
  ifstream inFile(fList.data());

  std::cout << fList << std::endl;
  std::cout << inFile.is_open() << std::endl;

  if(!inFile.is_open()){
    std::cout << "Error opening file. Exiting." <<std::endl;
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

  std::cout << "FileList Loaded" << std::endl;

  std::cout << "FileJob: " << listOfFiles[num] << std::endl;

  TFile* iniSkim_p = new TFile(listOfFiles[num].data(), "READ");

  GetDiJetIniSkim(iniSkim_p, sType, justJt);

  std::cout << "IniSkim Loaded" << std::endl;

  //Setup correction tables

  InitFactCorrFiles(sType, isHITrk);
  InitFactCorrHists(sType);

  
  TFile *histWeightFile_p = new TFile("histWeightFile.root", "READ");
  TH1F *hist_DataOverMC_p[11];

  if(sType == kHIMC){
    for(Int_t algIter = 0; algIter < 12; algIter++){
      hist_DataOverMC_p[algIter] = (TH1F*)histWeightFile_p->Get(Form("ak%s_dataHiBin_h", algType[algIter].c_str()));
      std::cout << Form("ak%s_dataHiBin_h", algType[algIter].c_str()) << std::endl;
    }
  }
  

  TFile* etaWeightFile_p = new TFile("etaWeightFile.root", "READ");
  TH1F* leadEtaWeightHist_p = (TH1F*)etaWeightFile_p->Get("leadEtaWeightHist_h");
  TH1F* subleadEtaWeightHist_p = (TH1F*)etaWeightFile_p->Get("subleadEtaWeightHist_h");
  
  std::string outName = listOfFiles[num];
  const std::string cutString = "/";
  const std::string iniString = "Ini";
  std::size_t strIndex = 0;

  std::cout << "Cull string" << std::endl;

  while(true){
    strIndex = outName.find(cutString);

    if(strIndex == std::string::npos) break;

    outName.replace(0, strIndex + 1, "");
  }

  std::cout << "Replace string" << std::endl;

  strIndex = outName.find(iniString);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, iniString.length(), "Ana"); 
  }

  std::cout << "Output name: " << outName.c_str() << std::endl;

  TFile *outFile = new TFile(outName.c_str(), "RECREATE");

  InitDiJetAnaSkim(sType, justJt);

  Long64_t nentries = jetTreeIni_p->GetEntries();

  std::cout << nentries << std::endl;

  for(Long64_t jentry = 0; jentry < nentries; jentry++){
    jetTreeIni_p->GetEntry(jentry);

    if(!justJt){
      trackTreeIni_p->GetEntry(jentry);
      
      if(montecarlo)
	genTreeIni_p->GetEntry(jentry);
    }      

    if(jentry%1000 == 0) std::cout << jentry << std::endl;

    InitJetVar(sType);

    getJtVar(nPu3Calo_, Pu3CaloPt_, Pu3CaloPhi_, Pu3CaloEta_, Pu3CaloTrkMax_, Pu3CaloRawPt_, Pu3CaloRefPt_, Pu3CaloRefPhi_, Pu3CaloRefEta_, 0, montecarlo, false);
    getJtVar(nPu4Calo_, Pu4CaloPt_, Pu4CaloPhi_, Pu4CaloEta_, Pu4CaloTrkMax_, Pu4CaloRawPt_, Pu4CaloRefPt_, Pu4CaloRefPhi_, Pu4CaloRefEta_, 1, montecarlo, false);
    getJtVar(nPu5Calo_, Pu5CaloPt_, Pu5CaloPhi_, Pu5CaloEta_, Pu5CaloTrkMax_, Pu5CaloRawPt_, Pu5CaloRefPt_, Pu5CaloRefPhi_, Pu5CaloRefEta_, 2, montecarlo, false);
    getJtVar(nVs2Calo_, Vs2CaloPt_, Vs2CaloPhi_, Vs2CaloEta_, Vs2CaloTrkMax_, Vs2CaloRawPt_, Vs2CaloRefPt_, Vs2CaloRefPhi_, Vs2CaloRefEta_, 3, montecarlo, false);
    getJtVar(nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_, Vs3CaloTrkMax_, Vs3CaloRawPt_, Vs3CaloRefPt_, Vs3CaloRefPhi_, Vs3CaloRefEta_, 4, montecarlo, false);
    getJtVar(nVs4Calo_, Vs4CaloPt_, Vs4CaloPhi_, Vs4CaloEta_, Vs4CaloTrkMax_, Vs4CaloRawPt_, Vs4CaloRefPt_, Vs4CaloRefPhi_, Vs4CaloRefEta_, 5, montecarlo, false);
    getJtVar(nVs5Calo_, Vs5CaloPt_, Vs5CaloPhi_, Vs5CaloEta_, Vs5CaloTrkMax_, Vs5CaloRawPt_, Vs5CaloRefPt_, Vs5CaloRefPhi_, Vs5CaloRefEta_, 6, montecarlo, false);

    getJtVar(nVs2CaloFrag_, Vs2CaloFragPt_, Vs2CaloFragPhi_, Vs2CaloFragEta_, Vs2CaloFragTrkMax_, Vs2CaloFragRawPt_, Vs2CaloFragRefPt_, Vs2CaloFragRefPhi_, Vs2CaloFragRefEta_, 7, montecarlo, false);
    getJtVar(nVs3CaloFrag_, Vs3CaloFragPt_, Vs3CaloFragPhi_, Vs3CaloFragEta_, Vs3CaloFragTrkMax_, Vs3CaloFragRawPt_, Vs3CaloFragRefPt_, Vs3CaloFragRefPhi_, Vs3CaloFragRefEta_, 8, montecarlo, false);
    getJtVar(nVs4CaloFrag_, Vs4CaloFragPt_, Vs4CaloFragPhi_, Vs4CaloFragEta_, Vs4CaloFragTrkMax_, Vs4CaloFragRawPt_, Vs4CaloFragRefPt_, Vs4CaloFragRefPhi_, Vs4CaloFragRefEta_, 9, montecarlo, false);
    getJtVar(nVs5CaloFrag_, Vs5CaloFragPt_, Vs5CaloFragPhi_, Vs5CaloFragEta_, Vs5CaloFragTrkMax_, Vs5CaloFragRawPt_, Vs5CaloFragRefPt_, Vs5CaloFragRefPhi_, Vs5CaloFragRefEta_, 10, montecarlo, false);
    getJtVar(nVs3CaloRes_, Vs3CaloResPt_, Vs3CaloResPhi_, Vs3CaloResEta_, Vs3CaloResTrkMax_, Vs3CaloResRawPt_, Vs3CaloResRefPt_, Vs3CaloResRefPhi_, Vs3CaloResRefEta_, 11, montecarlo, false);

    Float_t dummyArray2[nT2_];
    Float_t dummyArray3[nT3_];
    Float_t dummyArray4[nT4_];
    Float_t dummyArray5[nT5_];
    getJtVar(nT2_, T2Pt_, T2Phi_, T2Eta_, dummyArray2, dummyArray2, dummyArray2, dummyArray2, dummyArray2, 12, montecarlo, true);
    getJtVar(nT3_, T3Pt_, T3Phi_, T3Eta_, dummyArray3, dummyArray3, dummyArray3, dummyArray3, dummyArray3, 13, montecarlo, true);
    getJtVar(nT4_, T4Pt_, T4Phi_, T4Eta_, dummyArray4, dummyArray4, dummyArray4, dummyArray4, dummyArray4, 14, montecarlo, true);
    getJtVar(nT5_, T5Pt_, T5Phi_, T5Eta_, dummyArray5, dummyArray5, dummyArray5, dummyArray5, dummyArray5, 15, montecarlo, true);
    getJtVar(nPu3PF_, Pu3PFPt_, Pu3PFPhi_, Pu3PFEta_, Pu3PFTrkMax_, Pu3PFRawPt_, Pu3PFRefPt_, Pu3PFRefPhi_, Pu3PFRefEta_, 16, montecarlo, false);
    getJtVar(nVs3PF_, Vs3PFPt_, Vs3PFPhi_, Vs3PFEta_, Vs3PFTrkMax_, Vs3PFRawPt_, Vs3PFRefPt_, Vs3PFRefPhi_, Vs3PFRefEta_, 17, montecarlo, false);

    if(eventSet_[Pu3Calo] == false && eventSet_[Pu4Calo] == false && eventSet_[Pu5Calo] == false && eventSet_[Vs2Calo] == false && eventSet_[Vs3Calo] == false && eventSet_[Vs4Calo] == false && eventSet_[Vs5Calo] == false && eventSet_[Vs2CaloFrag] == false && eventSet_[Vs3CaloFrag] == false && eventSet_[Vs4CaloFrag] == false && eventSet_[Vs5CaloFrag] == false && eventSet_[Vs3CaloRes] == false && eventSet_[T2] == false && eventSet_[T3] == false && eventSet_[T4] == false && eventSet_[T5] == false && eventSet_[PuPF] == false && eventSet_[VsPF] == false){
      std::cout << "No event pass after IniSkim; Potential bug" << std::endl;

      for(Int_t iter = 0; iter < 10; iter++){
	std::cout << "eventSet " << iter << ": " << eventSet_[iter] << std::endl;
      }

      std::cout << "Pu3Calo Lead Jt: " << Pu3CaloPt_[0] << ", " << Pu3CaloPhi_[0] << ", " << Pu3CaloEta_[0] << std::endl;
      std::cout << "Pu3Calo subLead Jt: " << Pu3CaloPt_[1] << ", " << Pu3CaloPhi_[1] << ", " << Pu3CaloEta_[1] << std::endl;

      std::cout << "Pu4Calo Lead Jt: " << Pu4CaloPt_[0] << ", " << Pu4CaloPhi_[0] << ", " << Pu4CaloEta_[0] << std::endl;
      std::cout << "Pu4Calo subLead Jt: " << Pu4CaloPt_[1] << ", " << Pu4CaloPhi_[1] << ", " << Pu4CaloEta_[1] << std::endl;

      std::cout << "Pu5Calo Lead Jt: " << Pu5CaloPt_[0] << ", " << Pu5CaloPhi_[0] << ", " << Pu5CaloEta_[0] << std::endl;
      std::cout << "Pu5Calo subLead Jt: " << Pu5CaloPt_[1] << ", " << Pu5CaloPhi_[1] << ", " << Pu5CaloEta_[1] << std::endl;

      std::cout << "Vs2Calo Lead Jt: " << Vs2CaloPt_[0] << ", " << Vs2CaloPhi_[0] << ", " << Vs2CaloEta_[0] << std::endl;
      std::cout << "Vs2Calo subLead Jt: " << Vs2CaloPt_[1] << ", " << Vs2CaloPhi_[1] << ", " << Vs2CaloEta_[1] << std::endl;

      std::cout << "Vs3Calo Lead Jt: " << Vs3CaloPt_[0] << ", " << Vs3CaloPhi_[0] << ", " << Vs3CaloEta_[0] << std::endl;
      std::cout << "Vs3Calo subLead Jt: " << Vs3CaloPt_[1] << ", " << Vs3CaloPhi_[1] << ", " << Vs3CaloEta_[1] << std::endl;

      std::cout << "Vs4Calo Lead Jt: " << Vs4CaloPt_[0] << ", " << Vs4CaloPhi_[0] << ", " << Vs4CaloEta_[0] << std::endl;
      std::cout << "Vs4Calo subLead Jt: " << Vs4CaloPt_[1] << ", " << Vs4CaloPhi_[1] << ", " << Vs4CaloEta_[1] << std::endl;

      std::cout << "Vs5Calo Lead Jt: " << Vs5CaloPt_[0] << ", " << Vs5CaloPhi_[0] << ", " << Vs5CaloEta_[0] << std::endl;
      std::cout << "Vs5Calo subLead Jt: " << Vs5CaloPt_[1] << ", " << Vs5CaloPhi_[1] << ", " << Vs5CaloEta_[1] << std::endl;

      std::cout << "T Lead Jt: " << T3Pt_[0] << ", " << T3Phi_[0] << ", " << T3Eta_[0] << std::endl;
      std::cout << "T subLead Jt: " << T3Pt_[1] << ", " << T3Phi_[1] << ", " << T3Eta_[1] << std::endl;

      std::cout << "PuPF Lead Jt: " << Pu3PFPt_[0] << ", " << Pu3PFPhi_[0] << ", " << Pu3PFEta_[0] << std::endl;
      std::cout << "PuPF subLead Jt: " << Pu3PFPt_[1] << ", " << Pu3PFPhi_[1] << ", " << Pu3PFEta_[1] << std::endl;

      std::cout << "VsPF Lead Jt: " << Vs3PFPt_[0] << ", " << Vs3PFPhi_[0] << ", " << Vs3PFEta_[0] << std::endl;
      std::cout << "VsPF subLead Jt: " << Vs3PFPt_[1] << ", " << Vs3PFPhi_[1] << ", " << Vs3PFEta_[1] << std::endl;

      continue;
    }

    run_ = runIni_;
    evt_ = evtIni_;
    lumi_ = lumiIni_;

    if(montecarlo)
      pthat_ = pthatIni_;
    
    if(hi){
      hiBin_ = hiBinIni_;
      hiEvtPlane_ = hiEvtPlaneIni_;
      psin_ = psinIni_;
    }

    if(montecarlo){
      if(isHITrk){
	for(Int_t hatIter = 0; hatIter < 9; hatIter++){
	  if(pthat_ >= pthatCuts_PYTH_HITrk[hatIter] && pthat_ < pthatCuts_PYTH_HITrk[hatIter + 1]){
	    pthatWeight_ = pthatWeights_PYTH_HITrk[hatIter];
	    break;
	  }
	}
      }
      else if(hi){
	for(Int_t hatIter = 0; hatIter < 10; hatIter++){
          if(pthat_ >= pthatCuts_PYTH_HYD[hatIter] && pthat_ < pthatCuts_PYTH_HYD[hatIter + 1]){
            pthatWeight_ = pthatWeights_PYTH_HYD[hatIter];
            break;
          }
        }
      }
      else{
        for(Int_t hatIter = 0; hatIter < 6; hatIter++){
          if(pthat_ >= pthatCuts_PYTH_PPTrk[hatIter] && pthat_ < pthatCuts_PYTH_PPTrk[hatIter + 1]){
            pthatWeight_ = pthatWeights_PYTH_PPTrk[hatIter];
            break;
          }
        }
      }
    }

    if(TMath::Abs(AlgJtEta_[4][0]) && TMath::Abs(AlgJtEta_[4][1]) < 0.6){
	leadEtaWeight_ = leadEtaWeightHist_p->GetBinContent(leadEtaWeightHist_p->FindBin(AlgJtEta_[4][0]));
	subleadEtaWeight_ = subleadEtaWeightHist_p->GetBinContent(subleadEtaWeightHist_p->FindBin(AlgJtEta_[4][1]));
    }

    if(sType == kHIMC){
      for(Int_t algIter = 0; algIter < 12; algIter++){
	centWeight_[algIter] = hist_DataOverMC_p[algIter]->GetBinContent(hist_DataOverMC_p[algIter]->FindBin(hiBin_));
      }
      centWeight_[12] = hist_DataOverMC_p[3]->GetBinContent(hist_DataOverMC_p[3]->FindBin(hiBin_));
      centWeight_[13] = hist_DataOverMC_p[4]->GetBinContent(hist_DataOverMC_p[4]->FindBin(hiBin_));
      centWeight_[14] = hist_DataOverMC_p[5]->GetBinContent(hist_DataOverMC_p[5]->FindBin(hiBin_));
      centWeight_[15] = hist_DataOverMC_p[6]->GetBinContent(hist_DataOverMC_p[6]->FindBin(hiBin_));
    }
    
   //Iterate over tracks

    InitProjPerp(sType);

    //Switch below to iterated OR EDIT HERE

    if((eventSet_[Pu3Calo] || eventSet_[Pu4Calo] || eventSet_[Pu5Calo] || eventSet_[Vs2Calo] || eventSet_[Vs3Calo] || eventSet_[Vs4Calo] || eventSet_[Vs5Calo] || eventSet_[Vs2CaloFrag] || eventSet_[Vs3CaloFrag] || eventSet_[Vs4CaloFrag] || eventSet_[Vs5CaloFrag] || eventSet_[Vs3CaloRes] || eventSet_[T2] || eventSet_[T3] || eventSet_[T4] || eventSet_[T5]) && !justJt){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
        
	//Grab proj. Pt Spectra For Tracks in each Event Subset
	
	for(Int_t jtIter = 0; jtIter < nSumAlg; jtIter++){
	  if(eventSet_[jtIter]) GetTrkProjPerp(jtIter, jtIter, trkPt_[trkEntry], trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry]);
	}
      }

      if(hi) InitPosArrPbPb(hiBin_);

      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
	Int_t ptPos = getPtBin(trkPt_[trkEntry], sType);

	Float_t tempRMin[nSumAlg];
	Float_t tempFact[nSumAlg];
	Float_t tempCorr[nSumAlg];

	for(Int_t tempIter = 0; tempIter < nSumAlg; tempIter++){
	  tempRMin[tempIter] = 199.;
	  tempFact[tempIter] = 0.;
	  tempCorr[tempIter] = 0.;
	}

	if(eventSet_[Pu3Calo]) tempRMin[Pu3Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nPu3Calo_, Pu3CaloPt_, Pu3CaloPhi_, Pu3CaloEta_);
	if(eventSet_[Pu4Calo]) tempRMin[Pu4Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nPu3Calo_, Pu3CaloPt_, Pu3CaloPhi_, Pu3CaloEta_);
	if(eventSet_[Pu5Calo]) tempRMin[Pu5Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nPu3Calo_, Pu3CaloPt_, Pu3CaloPhi_, Pu3CaloEta_);
	if(eventSet_[Vs2Calo]) tempRMin[Vs2Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_);
	if(eventSet_[Vs3Calo]) tempRMin[Vs3Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_);
	if(eventSet_[Vs4Calo]) tempRMin[Vs4Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_);
	if(eventSet_[Vs5Calo]) tempRMin[Vs5Calo] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_);
	if(eventSet_[Vs2CaloFrag]) tempRMin[Vs2CaloFrag] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3CaloFrag_, Vs3CaloFragPt_, Vs3CaloFragPhi_, Vs3CaloFragEta_);
	if(eventSet_[Vs3CaloFrag]) tempRMin[Vs3CaloFrag] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3CaloFrag_, Vs3CaloFragPt_, Vs3CaloFragPhi_, Vs3CaloFragEta_);
	if(eventSet_[Vs4CaloFrag]) tempRMin[Vs4CaloFrag] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3CaloFrag_, Vs3CaloFragPt_, Vs3CaloFragPhi_, Vs3CaloFragEta_);
	if(eventSet_[Vs5CaloFrag]) tempRMin[Vs5CaloFrag] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3CaloFrag_, Vs3CaloFragPt_, Vs3CaloFragPhi_, Vs3CaloFragEta_);
	if(eventSet_[Vs3CaloRes]) tempRMin[Vs3CaloRes] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3CaloRes_, Vs3CaloResPt_, Vs3CaloResPhi_, Vs3CaloResEta_);
	if(eventSet_[T2]) tempRMin[T2] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_);
	if(eventSet_[T3]) tempRMin[T3] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_);
	if(eventSet_[T4]) tempRMin[T4] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_);
	if(eventSet_[T5]) tempRMin[T5] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_);

	if(eventSet_[Pu3Calo]){
	  tempFact[Pu3Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Pu3Calo], sType);
	  tempCorr[Pu3Calo] = trkPt_[trkEntry]*tempFact[Pu3Calo];
	}
	if(eventSet_[Pu4Calo]){
	  tempFact[Pu4Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Pu4Calo], sType);
	  tempCorr[Pu4Calo] = trkPt_[trkEntry]*tempFact[Pu4Calo];
	}
	if(eventSet_[Pu5Calo]){
	  tempFact[Pu5Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Pu5Calo], sType);
	  tempCorr[Pu5Calo] = trkPt_[trkEntry]*tempFact[Pu5Calo];
	}
	if(eventSet_[Vs2Calo]){
	  tempFact[Vs2Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs2Calo], sType);
	  tempCorr[Vs2Calo] = trkPt_[trkEntry]*tempFact[Vs2Calo];
	}
	if(eventSet_[Vs3Calo]){
	  tempFact[Vs3Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs3Calo], sType);
	  tempCorr[Vs3Calo] = trkPt_[trkEntry]*tempFact[Vs3Calo];
	}
	if(eventSet_[Vs4Calo]){
	  tempFact[Vs4Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs4Calo], sType);
	  tempCorr[Vs4Calo] = trkPt_[trkEntry]*tempFact[Vs4Calo];
	}
	if(eventSet_[Vs5Calo]){
	  tempFact[Vs5Calo] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs5Calo], sType);
	  tempCorr[Vs5Calo] = trkPt_[trkEntry]*tempFact[Vs5Calo];
	}
	if(eventSet_[Vs2CaloFrag]){
	  tempFact[Vs2CaloFrag] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs2CaloFrag], sType);
	  tempCorr[Vs2CaloFrag] = trkPt_[trkEntry]*tempFact[Vs2CaloFrag];
	}
	if(eventSet_[Vs3CaloFrag]){
	  tempFact[Vs3CaloFrag] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs3CaloFrag], sType);
	  tempCorr[Vs3CaloFrag] = trkPt_[trkEntry]*tempFact[Vs3CaloFrag];
	}
	if(eventSet_[Vs4CaloFrag]){
	  tempFact[Vs4CaloFrag] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs4CaloFrag], sType);
	  tempCorr[Vs4CaloFrag] = trkPt_[trkEntry]*tempFact[Vs4CaloFrag];
	}
	if(eventSet_[Vs5CaloFrag]){
	  tempFact[Vs5CaloFrag] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs5CaloFrag], sType);
	  tempCorr[Vs5CaloFrag] = trkPt_[trkEntry]*tempFact[Vs5CaloFrag];
	}
	if(eventSet_[Vs3CaloRes]){
	  tempFact[Vs3CaloRes] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[Vs3CaloRes], sType);
	  tempCorr[Vs3CaloRes] = trkPt_[trkEntry]*tempFact[Vs3CaloRes];
	}
	if(montecarlo){
	  if(eventSet_[T2]){
	    tempFact[T2] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[T2], sType);
	    tempCorr[T2] = trkPt_[trkEntry]*tempFact[T2];
	  }

	  if(eventSet_[T3]){
	    tempFact[T3] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[T3], sType);
	    tempCorr[T3] = trkPt_[trkEntry]*tempFact[T3];
	  }

	  if(eventSet_[T4]){
	    tempFact[T4] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[T4], sType);
	    tempCorr[T4] = trkPt_[trkEntry]*tempFact[T4];
	  }

	  if(eventSet_[T5]){
	    tempFact[T5] = factorizedPtCorr(ptPos, hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], tempRMin[T5], sType);
	    tempCorr[T5] = trkPt_[trkEntry]*tempFact[T5];
	  }
	}
	
	for(Int_t jtIter = 0; jtIter < nSumAlg; jtIter++){
	  if(eventSet_[jtIter]) GetTrkProjPerp(jtIter, jtIter+nSumAlg, trkPt_[trkEntry], tempCorr[jtIter], trkPhi_[trkEntry], trkEta_[trkEntry]);

	}	
      }

     if(montecarlo){
	//Iterate over Truth
	for(Int_t genEntry = 0; genEntry < nGen_; genEntry++){
	  for(Int_t jtIter = 0; jtIter < nSumAlg; jtIter++){
	    if(eventSet_[jtIter]) GetGenProjPerp(jtIter, genPt_[genEntry], genPhi_[genEntry], genEta_[genEntry]);
	  }
	}
      }

    }
    jetTreeAna_p->Fill();

    if(!justJt){
      trackTreeAna_p->Fill();
    
      if(montecarlo) genTreeAna_p->Fill();
    }
  }
  
  outFile->cd();
  jetTreeAna_p->Write("", TObject::kOverwrite);

  if(!justJt){
    trackTreeAna_p->Write("", TObject::kOverwrite);

    if(montecarlo) genTreeAna_p->Write("", TObject::kOverwrite);
  }

  /*  
  histWeightFile_p->Close();
  delete histWeightFile_p;
  */

  CleanupDiJetAnaSkim();
  outFile->Close();
  delete outFile;

  iniSkim_p->Close();
  delete iniSkim_p;

  printf("Done.\n");
  return(0);
}


int main(int argc, char *argv[])
{
  if(argc != 6)
    {
      std::cout << "Usage: sortForest <inputFile> <sampleType> <#> <justJtBool> <isHITrk>" << std::endl;
      return 1;
    }

  int rStatus = -1;

  rStatus = makeDiJetAnaSkim(argv[1], sampleType(atoi(argv[2])), atoi(argv[3]), Bool_t(atoi(argv[4])), Bool_t(atoi(argv[5])));

  return rStatus;
}
