//=============================================
// Author McGinn
// 
// DiJet Initial Skim Class (MC)
//
// !!NOTE: Written for jets sorted by pt, tracks unsorted!!
//
//=============================================

#include "/net/hisrv0001/home/cfmcginn/emDiJet/CMSSW_5_3_12_patch3/tempHIFA/HiForestAnalysis/hiForest.h"
#include "cfmDiJetIniSkim.h"
#include "stdlib.h"
#include <fstream>

#include "TComplex.h"

#include <vector>
#include "TLorentzVector.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"

const Float_t leadJtPtCut = 120.;
const Float_t subLeadJtPtCut = 50.;
const Float_t jtDelPhiCut = 0;
const Float_t jtEtaCut = 2.0; // Default Max at 2.4 to avoid transition junk, otherwise vary as needed

collisionType getCType(sampleType sType);

Double_t R = 0.3;
fastjet::JetAlgorithm algorithm = fastjet::antikt_algorithm;
fastjet::JetDefinition jetDef(algorithm, R, fastjet::E_scheme, fastjet::Best);

Bool_t passesDijet(Jets jtCollection, Int_t &lPtCut, Int_t &sLPtCut)
{
  Int_t leadJtIndex = -1;
  Int_t subLeadJtIndex = -1;

  if(jtCollection.nref == 0){
    lPtCut++;
    return false;
  }
  else if(jtCollection.nref == 1){
    sLPtCut++;
    return false;
  }

  for(Int_t jtEntry = 0; jtEntry < jtCollection.nref; jtEntry++){
    if(leadJtIndex < 0){
      if(jtCollection.jtpt[jtEntry] > leadJtPtCut){
	if(TMath::Abs(jtCollection.jteta[jtEntry]) < jtEtaCut)
	  leadJtIndex = jtEntry;
      }
      else{
	lPtCut++;
	return false;
      }
    }
    else if(subLeadJtIndex < 0){
      if(jtCollection.jtpt[jtEntry] > subLeadJtPtCut){
	if(TMath::Abs(jtCollection.jteta[jtEntry]) < jtEtaCut)
	  subLeadJtIndex = jtEntry;
      }
      else{
	sLPtCut++;
	return false;
      }
    }
    else
      return true;
  }

  if(leadJtIndex >= 0){
    if(subLeadJtIndex >= 0)
      return true;
    else{
      sLPtCut++;
      return false;
    }
  }
  else{
    lPtCut++;
    return false;
  }
}


void GetTrkJts(Int_t jtNum, sampleType sType, fastjet::PseudoJet inJtVect, Int_t &nConst, Float_t constPt[], Float_t constPhi[], Float_t constEta[], Float_t constCorr[])
{
  if(nConst != 0) return;

  if(getDR(inJtVect.eta(), inJtVect.phi_std(), Vs3CaloEta_[jtNum], Vs3CaloPhi_[jtNum]) > 0.3) return;

  std::vector<fastjet::PseudoJet> jtConst = inJtVect.constituents();

  nConst = (Int_t)(jtConst.size());
  TrkJtPt_[jtNum] = inJtVect.perp();
  TrkJtPhi_[jtNum] = inJtVect.phi_std();
  TrkJtEta_[jtNum] = inJtVect.eta();

  for(Int_t constIter = 0; constIter < nConst; constIter++){
    constPt[constIter] = jtConst[constIter].perp();
    constPhi[constIter] = jtConst[constIter].phi_std();
    constEta[constIter] = jtConst[constIter].eta();

    Int_t ptPos = getPtBin(jtConst[constIter].perp(), sType);
    Float_t tempRMin = getTrkRMin(jtConst[constIter].phi_std(), jtConst[constIter].eta(), nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_);

    constCorr[constIter] = factorizedPtCorr(ptPos, hiBinIni_, jtConst[constIter].perp(), jtConst[constIter].phi_std(), jtConst[constIter].eta(), tempRMin, sType);
  }
  return;
}



int makeDiJetIniSkim(string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_DIJETINISKIM.root", Int_t num = 0, Bool_t justJt = false)
{
  //Define MC or Data
  Bool_t montecarlo = isMonteCarlo(sType);
  Bool_t hi = isHI(sType);

  std::cout << sType << std::endl;
  std::cout << montecarlo << std::endl;

  collisionType cType = getCType(sType);

  string buffer;
  std::vector<string> listOfFiles;
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

  for(Int_t iter = 0; iter < (Int_t)(listOfFiles.size()); iter++){
    std::cout << listOfFiles[iter] << std::endl;
  }

  //Setup correction tables

  TFile *outFile = new TFile(Form("%s_%d.root", outName, num), "RECREATE");

  InitDiJetIniSkim(sType, justJt);

  if(justJt){
    InitCorrFiles(sType);
    InitCorrHists(sType);
  }

  std::cout << "JobFile: " << listOfFiles[num] << std::endl;

  HiForest *c = new HiForest(listOfFiles[num].data(), "Forest", cType, montecarlo);

  c->InitTree();

  c->LoadNoTrees();

  c->hasSkimTree = true;
  c->hasEvtTree = true;

  if(hi){
    c->hasAkPu3JetTree = true;
    c->hasAkVs3PFJetTree = true;
    c->hasAkPu3CaloJetTree = true;
    c->hasAkPu4CaloJetTree = true;
    c->hasAkPu5CaloJetTree = true;
    c->hasAkVs2CaloJetTree = true;
    c->hasAkVs3CaloJetTree = true;
    c->hasAkVs4CaloJetTree = true;
    c->hasAkVs5CaloJetTree = true;
  }
  else{
    c->hasAk2CaloJetTree = true;
    c->hasAk3CaloJetTree = true;
    c->hasAk4CaloJetTree = true;
    c->hasAk5CaloJetTree = true;
    c->hasAkPu3CaloJetTree = true;
    c->hasAkPu4CaloJetTree = true;
    c->hasAkPu5CaloJetTree = true;
  }

  std::cout << "TreeTruth: " << c->hasAk3CaloJetTree << std::endl;

  Float_t meanVz = 0;

  if(sType == kHIMC){
    //mean mc .16458, mean data -.337
    meanVz = .16458 + .337;
  }
  else if(sType == kPPMC){
    //MC vz = .4205,  Data vz = .6953
    meanVz = .4205 - .6953;
  }

  Long64_t nentries;

  if(!hi)
    nentries = c->ak3CaloJetTree->GetEntries();
  else
    nentries = c->GetEntries();

  std::cout << nentries << std::endl;

  Int_t totEv = 0;
  Int_t selectCut = 0;
  Int_t vzCut = 0;

  Int_t AlgLeadJtPtCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Int_t AlgSubLeadJtPtCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  std::cout << "Cuts, Lead/Sublead Pt, delphi, eta: " << leadJtPtCut << ", " << subLeadJtPtCut << ", " << jtDelPhiCut << ", " << jtEtaCut << std::endl; 

  for(Long64_t jentry = 0; jentry < nentries; jentry++){
    c->GetEntry(jentry);

    totEv++;

    if(jentry%1000 == 0)
      std::cout << jentry << std::endl;

    if(!c->selectEvent()){
      selectCut++;
      continue;
    }

    if(TMath::Abs(c->evt.vz - meanVz) > 15){
      vzCut++;
      continue;
    }

    //particle flow

    Jets AlgJtCollection[10];
    Int_t algMax = 10;

    if(hi){
      AlgJtCollection[0] = c->akPu3Calo;
      AlgJtCollection[1] = c->akPu4Calo;
      AlgJtCollection[2] = c->akPu5Calo;
      AlgJtCollection[3] = c->akVs2Calo;
      AlgJtCollection[4] = c->akVs3Calo;
      AlgJtCollection[5] = c->akVs4Calo;
      AlgJtCollection[6] = c->akVs5Calo;
      AlgJtCollection[8] = c->akPu3PF;
      AlgJtCollection[9] = c->akVs3PF;
    }
    else{
      AlgJtCollection[0] = c->akPu3Calo;
      AlgJtCollection[1] = c->akPu4Calo;
      AlgJtCollection[2] = c->akPu5Calo;
      AlgJtCollection[3] = c->ak2Calo;
      AlgJtCollection[4] = c->ak3Calo;
      AlgJtCollection[5] = c->ak4Calo;
      AlgJtCollection[6] = c->ak5Calo;
      algMax = 7;
    }

    Bool_t algPasses[10] = {false, false, false, false, false, false, false, false, false, false};

    for(Int_t algIter = 0; algIter < algMax; algIter++){
      if(algIter == 7)
	continue;

      algPasses[algIter] = passesDijet(AlgJtCollection[algIter], AlgLeadJtPtCut[algIter], AlgSubLeadJtPtCut[algIter]);
    }

    //truth, doesn't work w/ getLeadJt because truth doesnt get its own tree

    if(montecarlo){
      Int_t leadJtIndex = -1;
      Int_t subLeadJtIndex = -1;

      if(c->akPu3PF.ngen == 0){
	AlgLeadJtPtCut[7]++;
	algPasses[7] = false;
      }
      else if(c->akPu3PF.ngen == 1){
	AlgSubLeadJtPtCut[7]++;
	algPasses[7] = false;
      }
      else{
	for(Int_t jtEntry = 0; jtEntry < c->akPu3PF.ngen; jtEntry++){
	  if(leadJtIndex < 0){
	    if(c->akPu3PF.genpt[jtEntry] > leadJtPtCut){
	      if(TMath::Abs(c->akPu3PF.geneta[jtEntry]) < jtEtaCut)
		leadJtIndex = jtEntry;
	    }
	    else{
	      algPasses[7] = false;
	      break;
	    }
	  }
	  else if(subLeadJtIndex < 0){
	    if(c->akPu3PF.genpt[jtEntry] > subLeadJtPtCut){
	      if(TMath::Abs(c->akPu3PF.geneta[jtEntry]) < jtEtaCut)
		subLeadJtIndex = jtEntry;
	    }
	    else{
	      algPasses[7] = false;
	      break;
	    }
	  }
	  else{
	    algPasses[7] = true;
	    break;
	  }
	}
      }

      if(leadJtIndex >= 0){
	if(subLeadJtIndex >= 0)
	  algPasses[7] = true;
	else
	  algPasses[7] = false;
      }
      else{
	AlgLeadJtPtCut[7]++;
	algPasses[7] = false;
      }
    }
    
    if(algPasses[0] == false && algPasses[1] == false && algPasses[2] == false && algPasses[3] == false && algPasses[4] == false && algPasses[5] == false && algPasses[6] == false && algPasses[7] == false && algPasses[8] == false && algPasses[9] == false)
      continue;

    c->hasTrackTree = true;
    c->hasPFTree = true;
    if(montecarlo) c->hasGenParticleTree = true;

    /*
    Bool_t novaBool = false;
    for(Int_t contIter = 0; contIter < algMax; contIter++){
      if(AlgJtCollection[contIter].nref > 60) {
	novaBool = true;
	break;
      }
    }

    if(novaBool && !montecarlo) continue;
*/

    if(kHIMC == sType) pthatIni_ = c->akPu3PF.pthat;
    else if(kPPMC == sType) pthatIni_ = c->ak3Calo.pthat;

    if(hi){
      hiEvtPlaneIni_ = c->evt.hiEvtPlanes[21];                                                  
      TComplex cn1((c->pf.sumpt[0])*(c->pf.vn[2][0]), c->pf.psin[2][0], true);                    
      TComplex cn2((c->pf.sumpt[14])*(c->pf.vn[2][14]), c->pf.psin[2][14], true);                
      TComplex cn = cn1+cn2;                                                                    
      psinIni_ = cn.Theta();      
    }      

    runIni_ = c->evt.run;
    evtIni_ = c->evt.evt;
    lumiIni_ = c->evt.lumi;

    if(hi)
      hiBinIni_ = c->evt.hiBin;

    //Iterate over jets

    nPu3Calo_ = 0;
    nPu4Calo_ = 0;
    nPu5Calo_ = 0;
    nVs2Calo_ = 0;
    nVs3Calo_ = 0;
    nVs4Calo_ = 0;
    nVs5Calo_ = 0;
    nT3_ = 0;
    nPu3PF_ = 0;
    nVs3PF_ = 0;

    for(Int_t Pu3CaloIter = 0; Pu3CaloIter < AlgJtCollection[0].nref; Pu3CaloIter++){
      if(AlgJtCollection[0].jtpt[Pu3CaloIter] < 30.0)
	break;
      else if(TMath::Abs(AlgJtCollection[0].jteta[Pu3CaloIter]) > jtEtaCut)
	continue;

      Pu3CaloPt_[nPu3Calo_] = AlgJtCollection[0].jtpt[Pu3CaloIter];
      Pu3CaloPhi_[nPu3Calo_] = AlgJtCollection[0].jtphi[Pu3CaloIter];
      Pu3CaloEta_[nPu3Calo_] = AlgJtCollection[0].jteta[Pu3CaloIter];

      Pu3CaloTrkMax_[nPu3Calo_] = AlgJtCollection[0].trackMax[Pu3CaloIter];
      Pu3CaloRawPt_[nPu3Calo_] = AlgJtCollection[0].rawpt[Pu3CaloIter];

      if(montecarlo){
	Pu3CaloRefPt_[nPu3Calo_] = AlgJtCollection[0].refpt[Pu3CaloIter];
	Pu3CaloRefPhi_[nPu3Calo_] = AlgJtCollection[0].refphi[Pu3CaloIter];
	Pu3CaloRefEta_[nPu3Calo_] = AlgJtCollection[0].refeta[Pu3CaloIter];
	Pu3CaloRefPart_[nPu3Calo_] = AlgJtCollection[0].refparton_flavor[Pu3CaloIter];
      }

      nPu3Calo_++;
    }

    for(Int_t Pu4CaloIter = 0; Pu4CaloIter < AlgJtCollection[1].nref; Pu4CaloIter++){
      if(AlgJtCollection[1].jtpt[Pu4CaloIter] < 30.0)
	break;
      else if(TMath::Abs(AlgJtCollection[1].jteta[Pu4CaloIter]) > jtEtaCut)
	continue;

      Pu4CaloPt_[nPu4Calo_] = AlgJtCollection[1].jtpt[Pu4CaloIter];
      Pu4CaloPhi_[nPu4Calo_] = AlgJtCollection[1].jtphi[Pu4CaloIter];
      Pu4CaloEta_[nPu4Calo_] = AlgJtCollection[1].jteta[Pu4CaloIter];

      Pu4CaloTrkMax_[nPu4Calo_] = AlgJtCollection[1].trackMax[Pu4CaloIter];
      Pu4CaloRawPt_[nPu4Calo_] = AlgJtCollection[1].rawpt[Pu4CaloIter];

      if(montecarlo){
	Pu4CaloRefPt_[nPu4Calo_] = AlgJtCollection[1].refpt[Pu4CaloIter];
	Pu4CaloRefPhi_[nPu4Calo_] = AlgJtCollection[1].refphi[Pu4CaloIter];
	Pu4CaloRefEta_[nPu4Calo_] = AlgJtCollection[1].refeta[Pu4CaloIter];
	Pu4CaloRefPart_[nPu4Calo_] = AlgJtCollection[1].refparton_flavor[Pu4CaloIter];
      }

      nPu4Calo_++;
    }

    for(Int_t Pu5CaloIter = 0; Pu5CaloIter < AlgJtCollection[2].nref; Pu5CaloIter++){
      if(AlgJtCollection[2].jtpt[Pu5CaloIter] < 30.0)
	break;
      else if(TMath::Abs(AlgJtCollection[2].jteta[Pu5CaloIter]) > jtEtaCut)
	continue;

      Pu5CaloPt_[nPu5Calo_] = AlgJtCollection[2].jtpt[Pu5CaloIter];
      Pu5CaloPhi_[nPu5Calo_] = AlgJtCollection[2].jtphi[Pu5CaloIter];
      Pu5CaloEta_[nPu5Calo_] = AlgJtCollection[2].jteta[Pu5CaloIter];

      Pu5CaloTrkMax_[nPu5Calo_] = AlgJtCollection[2].trackMax[Pu5CaloIter];
      Pu5CaloRawPt_[nPu5Calo_] = AlgJtCollection[2].rawpt[Pu5CaloIter];

      if(montecarlo){
	Pu5CaloRefPt_[nPu5Calo_] = AlgJtCollection[2].refpt[Pu5CaloIter];
	Pu5CaloRefPhi_[nPu5Calo_] = AlgJtCollection[2].refphi[Pu5CaloIter];
	Pu5CaloRefEta_[nPu5Calo_] = AlgJtCollection[2].refeta[Pu5CaloIter];
	Pu5CaloRefPart_[nPu5Calo_] = AlgJtCollection[2].refparton_flavor[Pu5CaloIter];
      }

      nPu5Calo_++;
    }
    
    for(Int_t Vs2CaloIter = 0; Vs2CaloIter < AlgJtCollection[3].nref; Vs2CaloIter++){
      if(AlgJtCollection[3].jtpt[Vs2CaloIter] < 30.0)
	break;
      else if(TMath::Abs(AlgJtCollection[3].jteta[Vs2CaloIter]) > jtEtaCut)
	continue;
      
      Vs2CaloPt_[nVs2Calo_] = AlgJtCollection[3].jtpt[Vs2CaloIter];
      Vs2CaloPhi_[nVs2Calo_] = AlgJtCollection[3].jtphi[Vs2CaloIter];
      Vs2CaloEta_[nVs2Calo_] = AlgJtCollection[3].jteta[Vs2CaloIter];
      
      Vs2CaloTrkMax_[nVs2Calo_] = AlgJtCollection[3].trackMax[Vs2CaloIter];
      Vs2CaloRawPt_[nVs2Calo_] = AlgJtCollection[3].rawpt[Vs2CaloIter];
      
      if(montecarlo){
	Vs2CaloRefPt_[nVs2Calo_] = AlgJtCollection[3].refpt[Vs2CaloIter];
	Vs2CaloRefPhi_[nVs2Calo_] = AlgJtCollection[3].refphi[Vs2CaloIter];
	Vs2CaloRefEta_[nVs2Calo_] = AlgJtCollection[3].refeta[Vs2CaloIter];
	Vs2CaloRefPart_[nVs2Calo_] = AlgJtCollection[3].refparton_flavor[Vs2CaloIter];
      }
      
      nVs2Calo_++;
    }
       
    for(Int_t Vs3CaloIter = 0; Vs3CaloIter < AlgJtCollection[4].nref; Vs3CaloIter++){
      if(AlgJtCollection[4].jtpt[Vs3CaloIter] < 30.0)
	break;
      else if(TMath::Abs(AlgJtCollection[4].jteta[Vs3CaloIter]) > jtEtaCut)
	continue;
      
      Vs3CaloPt_[nVs3Calo_] = AlgJtCollection[4].jtpt[Vs3CaloIter];
      Vs3CaloPhi_[nVs3Calo_] = AlgJtCollection[4].jtphi[Vs3CaloIter];
      Vs3CaloEta_[nVs3Calo_] = AlgJtCollection[4].jteta[Vs3CaloIter];
      
      Vs3CaloTrkMax_[nVs3Calo_] = AlgJtCollection[4].trackMax[Vs3CaloIter];
      Vs3CaloRawPt_[nVs3Calo_] = AlgJtCollection[4].rawpt[Vs3CaloIter];

      if(montecarlo){
	Vs3CaloRefPt_[nVs3Calo_] = AlgJtCollection[4].refpt[Vs3CaloIter];
	Vs3CaloRefPhi_[nVs3Calo_] = AlgJtCollection[4].refphi[Vs3CaloIter];
	Vs3CaloRefEta_[nVs3Calo_] = AlgJtCollection[4].refeta[Vs3CaloIter];
	Vs3CaloRefPart_[nVs3Calo_] = AlgJtCollection[4].refparton_flavor[Vs3CaloIter];
      }
      
      nVs3Calo_++;
    }

    for(Int_t Vs4CaloIter = 0; Vs4CaloIter < AlgJtCollection[5].nref; Vs4CaloIter++){
      if(AlgJtCollection[5].jtpt[Vs4CaloIter] < 30.0)
	break;
      else if(TMath::Abs(AlgJtCollection[5].jteta[Vs4CaloIter]) > jtEtaCut)
	continue;

      Vs4CaloPt_[nVs4Calo_] = AlgJtCollection[5].jtpt[Vs4CaloIter];
      Vs4CaloPhi_[nVs4Calo_] = AlgJtCollection[5].jtphi[Vs4CaloIter];
      Vs4CaloEta_[nVs4Calo_] = AlgJtCollection[5].jteta[Vs4CaloIter];

      Vs4CaloTrkMax_[nVs4Calo_] = AlgJtCollection[5].trackMax[Vs4CaloIter];
      Vs4CaloRawPt_[nVs4Calo_] = AlgJtCollection[5].rawpt[Vs4CaloIter];

      if(montecarlo){
	Vs4CaloRefPt_[nVs4Calo_] = AlgJtCollection[5].refpt[Vs4CaloIter];
	Vs4CaloRefPhi_[nVs4Calo_] = AlgJtCollection[5].refphi[Vs4CaloIter];
	Vs4CaloRefEta_[nVs4Calo_] = AlgJtCollection[5].refeta[Vs4CaloIter];
	Vs4CaloRefPart_[nVs4Calo_] = AlgJtCollection[5].refparton_flavor[Vs4CaloIter];
      }

      nVs4Calo_++;
    }

    for(Int_t Vs5CaloIter = 0; Vs5CaloIter < AlgJtCollection[6].nref; Vs5CaloIter++){
      if(AlgJtCollection[6].jtpt[Vs5CaloIter] < 30.0)
	break;
      else if(TMath::Abs(AlgJtCollection[6].jteta[Vs5CaloIter]) > jtEtaCut)
	continue;

      Vs5CaloPt_[nVs5Calo_] = AlgJtCollection[6].jtpt[Vs5CaloIter];
      Vs5CaloPhi_[nVs5Calo_] = AlgJtCollection[6].jtphi[Vs5CaloIter];
      Vs5CaloEta_[nVs5Calo_] = AlgJtCollection[6].jteta[Vs5CaloIter];

      Vs5CaloTrkMax_[nVs5Calo_] = AlgJtCollection[6].trackMax[Vs5CaloIter];
      Vs5CaloRawPt_[nVs5Calo_] = AlgJtCollection[6].rawpt[Vs5CaloIter];

      if(montecarlo){
	Vs5CaloRefPt_[nVs5Calo_] = AlgJtCollection[6].refpt[Vs5CaloIter];
	Vs5CaloRefPhi_[nVs5Calo_] = AlgJtCollection[6].refphi[Vs5CaloIter];
	Vs5CaloRefEta_[nVs5Calo_] = AlgJtCollection[6].refeta[Vs5CaloIter];
	Vs5CaloRefPart_[nVs5Calo_] = AlgJtCollection[6].refparton_flavor[Vs5CaloIter];
      }

      nVs5Calo_++;
    }

    if(montecarlo){
      for(Int_t T3Iter = 0; T3Iter < c->akPu3PF.ngen; T3Iter++){
	if(c->akPu3PF.genpt[T3Iter] < 30.0)
	  break;
	else if(TMath::Abs(c->akPu3PF.geneta[T3Iter]) > jtEtaCut)
	  continue;
	
	T3Pt_[nT3_] = c->akPu3PF.genpt[T3Iter];
	T3Phi_[nT3_] = c->akPu3PF.genphi[T3Iter];
	T3Eta_[nT3_] = c->akPu3PF.geneta[T3Iter];
	T3Part_[nT3_] = c->akPu3PF.refparton_flavor[c->akPu3PF.genmatchindex[T3Iter]];
	
	nT3_++;
      }
    }
    
    if(hi){
      for(Int_t Pu3PFIter = 0; Pu3PFIter < AlgJtCollection[8].nref; Pu3PFIter++){
	if(AlgJtCollection[8].jtpt[Pu3PFIter] < 30.0)
	  break;
	else if(TMath::Abs(AlgJtCollection[8].jteta[Pu3PFIter]) > jtEtaCut)
	  continue;
	
	Pu3PFPt_[nPu3PF_] = AlgJtCollection[8].jtpt[Pu3PFIter];
	Pu3PFPhi_[nPu3PF_] = AlgJtCollection[8].jtphi[Pu3PFIter];
	Pu3PFEta_[nPu3PF_] = AlgJtCollection[8].jteta[Pu3PFIter];
	
	Pu3PFTrkMax_[nPu3PF_] = AlgJtCollection[8].trackMax[Pu3PFIter];
	Pu3PFRawPt_[nPu3PF_] = AlgJtCollection[8].rawpt[Pu3PFIter];
	
	if(montecarlo){
	  Pu3PFRefPt_[nPu3PF_] = AlgJtCollection[8].refpt[Pu3PFIter];
	  Pu3PFRefPhi_[nPu3PF_] = AlgJtCollection[8].refphi[Pu3PFIter];
	  Pu3PFRefEta_[nPu3PF_] = AlgJtCollection[8].refeta[Pu3PFIter];
	  Pu3PFRefPart_[nPu3PF_] = AlgJtCollection[8].refparton_flavor[Pu3PFIter];
	}
	
	nPu3PF_++;
      }
      
      for(Int_t Vs3PFIter = 0; Vs3PFIter < AlgJtCollection[9].nref; Vs3PFIter++){
	if(AlgJtCollection[9].jtpt[Vs3PFIter] < 30.0)
	  break;
	else if(TMath::Abs(AlgJtCollection[9].jteta[Vs3PFIter]) > jtEtaCut)
	  continue;
	
	Vs3PFPt_[nVs3PF_] = AlgJtCollection[9].jtpt[Vs3PFIter];
	Vs3PFPhi_[nVs3PF_] = AlgJtCollection[9].jtphi[Vs3PFIter];
	Vs3PFEta_[nVs3PF_] = AlgJtCollection[9].jteta[Vs3PFIter];
	
	Vs3PFTrkMax_[nVs3PF_] = AlgJtCollection[9].trackMax[Vs3PFIter];
	Vs3PFRawPt_[nVs3PF_] = AlgJtCollection[9].rawpt[Vs3PFIter];

	if(montecarlo){
	  Vs3PFRefPt_[nVs3PF_] = AlgJtCollection[9].refpt[Vs3PFIter];
	  Vs3PFRefPhi_[nVs3PF_] = AlgJtCollection[9].refphi[Vs3PFIter];
	  Vs3PFRefEta_[nVs3PF_] = AlgJtCollection[9].refeta[Vs3PFIter];
	  Vs3PFRefPart_[nVs3PF_] = AlgJtCollection[9].refparton_flavor[Vs3PFIter];
	}
	
	nVs3PF_++;
      }
    }

    //Iterate over tracks

    nTrk_ = 0;

    Tracks trkCollection;
    trkCollection = c->track;

    std::vector<fastjet::PseudoJet>* jtVect_p = new std::vector<fastjet::PseudoJet>;
    TLorentzVector tempTL;
    
    for(Int_t trkEntry = 0; trkEntry < trkCollection.nTrk; trkEntry++){
      
      if(TMath::Abs(trkCollection.trkEta[trkEntry]) > 2.4)
	continue;
      
      if(trkCollection.trkPt[trkEntry] <= 0.5)
	continue;
	
      if(!trkCollection.highPurity[trkEntry])
	continue;
	
      if(TMath::Abs(trkCollection.trkDz1[trkEntry]/trkCollection.trkDzError1[trkEntry]) > 3)
	continue;
	
      if(TMath::Abs(trkCollection.trkDxy1[trkEntry]/trkCollection.trkDxyError1[trkEntry]) > 3)
	continue;
	
      if(trkCollection.trkPtError[trkEntry]/trkCollection.trkPt[trkEntry] > 0.1)
	continue;
	
      trkPt_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPhi_[nTrk_] = trkCollection.trkPhi[trkEntry];
      trkEta_[nTrk_] = trkCollection.trkEta[trkEntry];
      tempTL.SetPtEtaPhiM(trkCollection.trkPt[trkEntry], trkCollection.trkEta[trkEntry], trkCollection.trkPhi[trkEntry], 0);
      jtVect_p->push_back(tempTL);
	
      //Grab proj. Pt Spectra For Tracks in each Event Subset    
	
      nTrk_++;
      if(nTrk_ > maxTracks - 1){
	printf("ERROR: Trk arrays not large enough.\n");
	return(1);
      }
    }

    nPF_ = 0;

    PFs pf = c->pf;

    for(Int_t pfIter = 0; pfIter < pf.nPFpart; pfIter++){
      if(TMath::Abs(pf.pfEta[pfIter]) > 2.3) continue;

      if(pf.pfPt[pfIter] < 0.5) continue;

      pfPt_[nPF_] = pf.pfPt[pfIter];
      pfVsIniPt_[nPF_] = pf.pfVsPtInitial[pfIter];
      pfPhi_[nPF_] = pf.pfPhi[pfIter];
      pfEta_[nPF_] = pf.pfEta[pfIter];

      nPF_++;
      if(nPF_ > maxPF - 1){
	printf("ERROR: PF arrays not large enough.\n");
	return(1);
      }
    }

    if(justJt){
      nLeadJtConst_ = 0;
      nSubLeadJtConst_ = 0;
      nThirdJtConst_ = 0;
      nFourthJtConst_ = 0;
      nFifthJtConst_ = 0;
    }

    InitTrkJts(justJt);
    if(hi)
      InitPosArrPbPb(hiBinIni_);

    if(justJt && algPasses[1] == true){
      fastjet::ClusterSequence cs(*jtVect_p, jetDef);
      std::vector<fastjet::PseudoJet> jtVectSort = sorted_by_pt(cs.inclusive_jets());

      for(Int_t jtIter = 0; jtIter < (Int_t)(jtVectSort.size()); jtIter++){
	GetTrkJts(0, sType, jtVectSort[jtIter], nLeadJtConst_, TrkLeadJtConstPt_, TrkLeadJtConstPhi_, TrkLeadJtConstEta_, TrkLeadJtConstCorr_);
	GetTrkJts(1, sType, jtVectSort[jtIter], nSubLeadJtConst_, TrkSubLeadJtConstPt_, TrkSubLeadJtConstPhi_, TrkSubLeadJtConstEta_, TrkSubLeadJtConstCorr_);
	GetTrkJts(2, sType, jtVectSort[jtIter], nThirdJtConst_, TrkThirdJtConstPt_, TrkThirdJtConstPhi_, TrkThirdJtConstEta_, TrkThirdJtConstCorr_);
	GetTrkJts(3, sType, jtVectSort[jtIter], nFourthJtConst_, TrkFourthJtConstPt_, TrkFourthJtConstPhi_, TrkFourthJtConstEta_, TrkFourthJtConstCorr_);
	GetTrkJts(4, sType, jtVectSort[jtIter], nFifthJtConst_, TrkFifthJtConstPt_, TrkFifthJtConstPhi_, TrkFifthJtConstEta_, TrkFifthJtConstCorr_);
      }
      jtVectSort.clear();
    }
    
    if(!justJt){
      if(montecarlo){
	//Iterate over truth
	nGen_ = 0;
	
	GenParticles genCollection;
	genCollection = c->genparticle;
	
	for(Int_t genEntry = 0; genEntry < genCollection.mult; genEntry++){
	  if(genCollection.chg[genEntry] == 0)
	    continue;
	
	  if(TMath::Abs(genCollection.eta[genEntry]) > 2.4)
	    continue;
	  
	  if(genCollection.pt[genEntry] < 0.5)
	    continue;
	
	  genPt_[nGen_] = genCollection.pt[genEntry];
	  genPhi_[nGen_] = genCollection.phi[genEntry];
	  genEta_[nGen_] = genCollection.eta[genEntry];
	  
	  nGen_++;
	  if(nGen_ > maxEntrySim - 1){
	    printf("ERROR: Gen arrays not large enough.\n");
	    return(1);
	  }
	}
      }
    }
    
    jetTreeIni_p->Fill();

    if(!justJt){
      trackTreeIni_p->Fill();
      pfCandTreeIni_p->Fill();
    
      if(montecarlo) genTreeIni_p->Fill();
    }

    jtVect_p->clear();
    delete jtVect_p;
  }

  std::cout << "totEv: " << totEv << std::endl;
  Int_t tempTot = totEv - selectCut;
  std::cout << "selectCut: " << tempTot << std::endl;
  tempTot = tempTot - vzCut;
  std::cout << "vzCut: " << tempTot << std::endl;

  for(Int_t cutIter = 0; cutIter < 10; cutIter++){
    std::cout << std::endl;
    tempTot = totEv - selectCut - vzCut - AlgLeadJtPtCut[cutIter];
    std::cout << "AlgLeadJtPtCut[" << cutIter << "]: " << tempTot << std::endl;
    tempTot = tempTot - AlgSubLeadJtPtCut[cutIter];
    std::cout << "AlgSubLeadJtPtCut[" << cutIter << "]: " << tempTot << std::endl;
  }

  outFile->cd();

  jetTreeIni_p->Write("", TObject::kOverwrite);

  if(!justJt){
    trackTreeIni_p->Write("", TObject::kOverwrite);
    pfCandTreeIni_p->Write("", TObject::kOverwrite);
    
    if(montecarlo) genTreeIni_p->Write("", TObject::kOverwrite);
  }

  delete c;
  CleanupDiJetIniSkim();
  outFile->Close();
  delete outFile;

  printf("Done.\n");
  return(0);
}


collisionType getCType(sampleType sType)
{
  switch (sType)
    {
    case kPPDATA:
    case kPPMC:
      return cPP;
    case kPADATA:
    case kPAMC:
      return cPPb;
    case kHIDATA:
    case kHIMC:
      return cPbPb;
    }
  return cPbPb; //probably a bad guess
}


int main(int argc, char *argv[])
{
  if(argc != 6)
    {
      std::cout << "Usage: makeDiJetIniSkim <inputFile> <sType> <outputFile> <#> <justJtBool>" << std::endl;
      std::cout << argc << std::endl;
      std::cout << argv[0] << ", "  << argv[1] << ", " << argv[2] << ", " << argv[3]<< ", " << argv[4] << ", " << argv[5] << std::endl;
      return 1;
    }

  int rStatus = -1;

  rStatus = makeDiJetIniSkim(argv[1], sampleType(atoi(argv[2])), argv[3], atoi(argv[4]), Bool_t(atoi(argv[5])));

  return rStatus;
}
