//=============================================
// Author: Chris McGinn
// 
// DiJet Initial Skim Class (MC)
//
// !!NOTE: Written for jets sorted by pt, tracks unsorted!!
//
//=============================================

#include "TFile.h"
#include "TTree.h"
#include "/net/hisrv0001/home/cfmcginn/emDiJet/CMSSW_5_3_12_patch3/tempHIFA/HiForestAnalysis/hiForest.h"
#include "cfmDiJetIniSkim.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>

#include "TComplex.h"

const Float_t leadJtPtCut = 120.;
const Float_t subLeadJtPtCut = 50.;
const Float_t jtDelPhiCut = 0;
const Float_t jtEtaCut = 2.0; // Default Max at 2.4 to avoid transition junk, otherwise vary as needed

collisionType getCType(sampleType sType);

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


int makeDiJetIniSkim(string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_DIJETINISKIM.root", Int_t num = 0)
{
  //Define MC or Data
  Bool_t montecarlo = false;
  if(sType == kPPMC || sType == kPAMC || sType == kHIMC)
    montecarlo = true;

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
    while(!inFile.eof()){
      inFile >> buffer;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }

  std::cout << "FileList Loaded" << std::endl;

  //Setup correction tables

  TFile *outFile = new TFile(Form("%s_%d.root", outName, num), "RECREATE");

  InitDiJetIniSkim(montecarlo, sType);

  HiForest *c = new HiForest(listOfFiles[0].data(), "Forest", cType, montecarlo);

  c->InitTree();

  c->LoadNoTrees();

  c->hasSkimTree = true;
  c->hasTrackTree = true;
  c->hasEvtTree = true;

  if(sType == kHIDATA || sType == kHIMC){
    c->hasAkPu3JetTree = true;
    c->hasAkVs3PFJetTree = true;
    c->hasAkPu3CaloJetTree = true;
    c->hasAkVs3CaloJetTree = true;
  }
  else if(sType == kPPDATA || sType == kPPMC){
    c->hasAk3CaloJetTree = true;
  }

  std::cout << "TreeTruth: " << c->hasAk3CaloJetTree << std::endl;

  Float_t meanVz = 0;

  if(sType == kHIMC){
    c->hasGenParticleTree = true;
    //mean mc .16458, mean data -.337
    meanVz = .16458 + .337;
  }
  else if(sType == kPPMC){
    c->hasGenParticleTree = true;
    //MC vz = .4205,  Data vz = .6953
    meanVz = .4205 - .6953;
  }

  Long64_t nentries;

  if(sType == kPPDATA || sType == kPPMC)
    nentries = c->ak3CaloJetTree->GetEntries();
  else
    nentries = c->GetEntries();

  std::cout << nentries << std::endl;

  Int_t totEv = 0;
  Int_t selectCut = 0;
  Int_t vzCut = 0;

  Int_t AlgLeadJtPtCut[5] = {0, 0, 0, 0, 0};
  Int_t AlgSubLeadJtPtCut[5] = {0, 0, 0, 0, 0};

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

    Jets AlgJtCollection[5];
    Int_t algMax = 5;

    if(sType == kHIDATA || sType == kHIMC){
      AlgJtCollection[0] = c->akPu3Calo;
      AlgJtCollection[1] = c->akVs3Calo;
      AlgJtCollection[3] = c->akPu3PF;
      AlgJtCollection[4] = c->akVs3PF;
    }
    else if(sType == kPPDATA || sType == kPPMC){
      AlgJtCollection[0] = c->akPu3Calo;
      AlgJtCollection[1] = c->ak3Calo;
      algMax = 2;
    }

    Bool_t algPasses[5] = {false, false, false, false, false};

    for(Int_t algIter = 0; algIter < algMax; algIter++){
      if(algIter == 2)
	continue;

      algPasses[algIter] = passesDijet(AlgJtCollection[algIter], AlgLeadJtPtCut[algIter], AlgSubLeadJtPtCut[algIter]);
    }

    //truth, doesn't work w/ getLeadJt because truth doesnt get its own tree

    if(montecarlo){
      Int_t leadJtIndex = -1;
      Int_t subLeadJtIndex = -1;

      if(c->akPu3PF.ngen == 0){
	AlgLeadJtPtCut[2]++;
	algPasses[2] = false;
      }
      else if(c->akPu3PF.ngen == 1){
	AlgSubLeadJtPtCut[2]++;
	algPasses[2] = false;
      }
      else{
	for(Int_t jtEntry = 0; jtEntry < c->akPu3PF.ngen; jtEntry++){
	  if(leadJtIndex < 0){
	    if(c->akPu3PF.genpt[jtEntry] > leadJtPtCut){
	      if(TMath::Abs(c->akPu3PF.geneta[jtEntry]) < jtEtaCut)
		leadJtIndex = jtEntry;
	    }
	    else{
	      algPasses[2] = false;
	      break;
	    }
	  }
	  else if(subLeadJtIndex < 0){
	    if(c->akPu3PF.genpt[jtEntry] > subLeadJtPtCut){
	      if(TMath::Abs(c->akPu3PF.geneta[jtEntry]) < jtEtaCut)
		subLeadJtIndex = jtEntry;
	    }
	    else{
	      algPasses[2] = false;
	      break;
	    }
	  }
	  else{
	    algPasses[2] = true;
	    break;
	  }
	}
      }

      if(leadJtIndex >= 0){
	if(subLeadJtIndex >= 0)
	  algPasses[2] = true;
	else
	  algPasses[2] = false;
      }
      else{
	AlgLeadJtPtCut[2]++;
	algPasses[2] = false;
      }
    }
    
    if(algPasses[0] == false && algPasses[1] == false && algPasses[2] == false && algPasses[3] == false && algPasses[4] == false)
      continue;

    if(kHIMC) pthatIni_ = c->akPu3PF.pthat;
    else if(kPPMC) pthatIni_ = c->ak3Calo.pthat;

    if(sType == kHIDATA || sType == kHIMC){
      hiEvtPlaneIni_ = c->evt.hiEvtPlanes[21];                                                  
      TComplex cn1((c->pf.sumpt[0])*(c->pf.vn[2][0]), c->pf.psin[2][0], true);                    
      TComplex cn2((c->pf.sumpt[14])*(c->pf.vn[2][14]), c->pf.psin[2][14], true);                
      TComplex cn = cn1+cn2;                                                                    
      psinIni_ = cn.Theta();      
    }      

    runIni_ = c->evt.run;
    evtIni_ = c->evt.evt;
    lumiIni_ = c->evt.lumi;

    if(sType == kHIDATA || sType == kHIMC)
      hiBinIni_ = c->evt.hiBin;

    //Iterate over jets

    nPu3Calo_ = 0;
    nVs3Calo_ = 0;
    nT3_ = 0;
    nPu3PF_ = 0;
    nVs3PF_ = 0;

    for(Int_t Pu3CaloIter = 0; Pu3CaloIter < AlgJtCollection[0].nref; Pu3CaloIter++){
      if(AlgJtCollection[0].jtpt[Pu3CaloIter] < subLeadJtPtCut)
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
      }

      nPu3Calo_++;
    }

    for(Int_t Vs3CaloIter = 0; Vs3CaloIter < AlgJtCollection[1].nref; Vs3CaloIter++){
      if(AlgJtCollection[1].jtpt[Vs3CaloIter] < subLeadJtPtCut)
	break;
      else if(TMath::Abs(AlgJtCollection[1].jteta[Vs3CaloIter]) > jtEtaCut)
	continue;

      Vs3CaloPt_[nVs3Calo_] = AlgJtCollection[1].jtpt[Vs3CaloIter];
      Vs3CaloPhi_[nVs3Calo_] = AlgJtCollection[1].jtphi[Vs3CaloIter];
      Vs3CaloEta_[nVs3Calo_] = AlgJtCollection[1].jteta[Vs3CaloIter];

      Vs3CaloTrkMax_[nVs3Calo_] = AlgJtCollection[1].trackMax[Vs3CaloIter];
      Vs3CaloRawPt_[nVs3Calo_] = AlgJtCollection[1].rawpt[Vs3CaloIter];

      if(montecarlo){
	Vs3CaloRefPt_[nVs3Calo_] = AlgJtCollection[1].refpt[Vs3CaloIter];
	Vs3CaloRefPhi_[nVs3Calo_] = AlgJtCollection[1].refphi[Vs3CaloIter];
	Vs3CaloRefEta_[nVs3Calo_] = AlgJtCollection[1].refeta[Vs3CaloIter];
      }

      nVs3Calo_++;
    }

    if(montecarlo){
      for(Int_t T3Iter = 0; T3Iter < c->akPu3PF.ngen; T3Iter++){
	if(c->akPu3PF.genpt[T3Iter] < subLeadJtPtCut)
	  break;
	else if(TMath::Abs(c->akPu3PF.geneta[T3Iter]) > jtEtaCut)
	  continue;
	
	T3Pt_[nT3_] = c->akPu3PF.genpt[T3Iter];
	T3Phi_[nT3_] = c->akPu3PF.genphi[T3Iter];
	T3Eta_[nT3_] = c->akPu3PF.geneta[T3Iter];
	
	nT3_++;
      }
    }

    for(Int_t Pu3PFIter = 0; Pu3PFIter < AlgJtCollection[3].nref; Pu3PFIter++){
      if(AlgJtCollection[3].jtpt[Pu3PFIter] < subLeadJtPtCut)
	break;
      else if(TMath::Abs(AlgJtCollection[3].jteta[Pu3PFIter]) > jtEtaCut)
	continue;

      Pu3PFPt_[nPu3PF_] = AlgJtCollection[3].jtpt[Pu3PFIter];
      Pu3PFPhi_[nPu3PF_] = AlgJtCollection[3].jtphi[Pu3PFIter];
      Pu3PFEta_[nPu3PF_] = AlgJtCollection[3].jteta[Pu3PFIter];

      Pu3PFTrkMax_[nPu3PF_] = AlgJtCollection[3].trackMax[Pu3PFIter];
      Pu3PFRawPt_[nPu3PF_] = AlgJtCollection[3].rawpt[Pu3PFIter];

      if(montecarlo){
	Pu3PFRefPt_[nPu3PF_] = AlgJtCollection[3].refpt[Pu3PFIter];
	Pu3PFRefPhi_[nPu3PF_] = AlgJtCollection[3].refphi[Pu3PFIter];
	Pu3PFRefEta_[nPu3PF_] = AlgJtCollection[3].refeta[Pu3PFIter];
      }

      nPu3PF_++;
    }

    for(Int_t Vs3PFIter = 0; Vs3PFIter < AlgJtCollection[4].nref; Vs3PFIter++){
      if(AlgJtCollection[4].jtpt[Vs3PFIter] < subLeadJtPtCut)
	break;
      else if(TMath::Abs(AlgJtCollection[4].jteta[Vs3PFIter]) > jtEtaCut)
	continue;

      Vs3PFPt_[nVs3PF_] = AlgJtCollection[4].jtpt[Vs3PFIter];
      Vs3PFPhi_[nVs3PF_] = AlgJtCollection[4].jtphi[Vs3PFIter];
      Vs3PFEta_[nVs3PF_] = AlgJtCollection[4].jteta[Vs3PFIter];

      Vs3PFTrkMax_[nVs3PF_] = AlgJtCollection[4].trackMax[Vs3PFIter];
      Vs3PFRawPt_[nVs3PF_] = AlgJtCollection[4].rawpt[Vs3PFIter];

      if(montecarlo){
	Vs3PFRefPt_[nVs3PF_] = AlgJtCollection[4].refpt[Vs3PFIter];
	Vs3PFRefPhi_[nVs3PF_] = AlgJtCollection[4].refphi[Vs3PFIter];
	Vs3PFRefEta_[nVs3PF_] = AlgJtCollection[4].refeta[Vs3PFIter];
      }

      nVs3PF_++;
    }

    //Iterate over tracks

    nTrk_ = 0;

    Tracks trkCollection;
    trkCollection = c->track;
    
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
      
      //Grab proj. Pt Spectra For Tracks in each Event Subset    
      
      nTrk_++;
      if(nTrk_ > MAXTRKS - 1){
	printf("ERROR: Trk arrays not large enough.\n");
	return(1);
      }
    }
    
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
	if(nGen_ > MAXGEN - 1){
	  printf("ERROR: Gen arrays not large enough.\n");
	  return(1);
	}
      }
    }
    
    jetTreeIni_p->Fill();
    trackTreeIni_p->Fill();
    
    if(montecarlo)
      genTreeIni_p->Fill();
  }

  std::cout << "totEv: " << totEv << std::endl;
  Int_t tempTot = totEv - selectCut;
  std::cout << "selectCut: " << tempTot << std::endl;
  tempTot = tempTot - vzCut;
  std::cout << "vzCut: " << tempTot << std::endl;

  for(Int_t cutIter = 0; cutIter < 3; cutIter++){
    std::cout << std::endl;
    tempTot = totEv - selectCut - vzCut - AlgLeadJtPtCut[cutIter];
    std::cout << "AlgLeadJtPtCut[" << cutIter << "]: " << tempTot << std::endl;
    tempTot = tempTot - AlgSubLeadJtPtCut[cutIter];
    std::cout << "AlgSubLeadJtPtCut[" << cutIter << "]: " << tempTot << std::endl;
  }

  outFile->cd();

  jetTreeIni_p->Write("", TObject::kOverwrite);
  trackTreeIni_p->Write("", TObject::kOverwrite);

  if(montecarlo)
    genTreeIni_p->Write("", TObject::kOverwrite);

  outFile->Close();

  delete c;
  CleanupDiJetIniSkim(montecarlo);
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
  if(argc != 5)
    {
      std::cout << "Usage: makeDiJetIniSkim <inputFile> <MCBool> <outputFile> <#>" << std::endl;
      return 1;
    }

  int rStatus = -1;

  rStatus = makeDiJetIniSkim(argv[1], sampleType(atoi(argv[2])), argv[3], atoi(argv[4]));

  return rStatus;
}
