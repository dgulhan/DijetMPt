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

const Int_t nEvtPerFile = 10000;

const Float_t jtAlgR[17] = {0.3, 0.4, 0.5, 0.2, 0.3, 0.4, 0.5, 0.2, 0.3, 0.4, 0.5, 0.2, 0.3, 0.4, 0.5, 0.3, 0.3};
const Float_t jtAlgRBin[17] = {1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 1, 1};

collisionType getCType(sampleType sType);

Double_t R = 0.3;
fastjet::JetAlgorithm algorithm = fastjet::antikt_algorithm;
fastjet::JetDefinition jetDef(algorithm, R, fastjet::E_scheme, fastjet::Best);

void setJtBranches(TTree* inJtTree, Bool_t montecarlo = false, Bool_t isGen = false)
{
  inJtTree->SetBranchStatus("*", 0);
  inJtTree->SetBranchStatus("nref", 1);
  inJtTree->SetBranchStatus("jtpt", 1);
  inJtTree->SetBranchStatus("rawpt", 1);
  inJtTree->SetBranchStatus("jteta", 1);
  inJtTree->SetBranchStatus("jtphi", 1);
  inJtTree->SetBranchStatus("trackMax", 1);
  if(montecarlo){
    inJtTree->SetBranchStatus("pthat", 1);
    inJtTree->SetBranchStatus("refparton_flavor", 1);
    inJtTree->SetBranchStatus("refpt", 1);
    inJtTree->SetBranchStatus("refeta", 1);
    inJtTree->SetBranchStatus("refphi", 1);
    if(isGen){
      inJtTree->SetBranchStatus("ngen", 1);
      inJtTree->SetBranchStatus("genpt", 1);
      inJtTree->SetBranchStatus("geneta", 1);
      inJtTree->SetBranchStatus("genphi", 1);
      inJtTree->SetBranchStatus("genmatchindex", 1);
    }
  }
  return;
}


void getCorrJtCollection(Int_t algNum, Jets* inJt, PFs *pf, sampleType sType, Int_t hiBin)
{
  for(Int_t jtIter = 0; jtIter < inJt->nref; jtIter++){
    if(inJt->jtpt[jtIter] < 20) break;

    Int_t nPF = 0;
    if(sType == kHIDATA || sType == kHIMC) nPF = Get2PFCand(jtAlgR[algNum], inJt->jtphi[jtIter], inJt->jteta[jtIter], pf->nPFpart, pf->pfVsPt, pf->pfId, pf->pfPhi, pf->pfEta);
    else if(sType == kPPDATA || sType == kPPMC) nPF = Get2PFCand(jtAlgR[algNum], inJt->jtphi[jtIter], inJt->jteta[jtIter], pf->nPFpart, pf->pfPt, pf->pfId, pf->pfPhi, pf->pfEta);

    inJt->jtpt[jtIter] = inJt->jtpt[jtIter]*GetJtCorrPt(sType, jtAlgRBin[algNum], hiBin, inJt->jtpt[jtIter], nPF);
  }

  Int_t sortIter = 0;

  Float_t tempPt = 0;
  Float_t tempPhi = 0;
  Float_t tempEta = 0;
  Float_t tempTrkMax = 0;
  Float_t tempRawPt = 0;
  Float_t tempRefPt = 0;
  Float_t tempRefPhi = 0;
  Float_t tempRefEta = 0;
  Int_t tempRefPart = 0;

  while(sortIter < inJt->nref){
    Int_t maxIndex = 0;
    Float_t maxPt = 0;

    for(Int_t jtIter = sortIter; jtIter < inJt->nref; jtIter++){
      if(inJt->jtpt[jtIter] > maxPt){
	maxPt = inJt->jtpt[jtIter];
	maxIndex = jtIter;
      }
    }

    if(maxPt < 30) break;

    tempPt = inJt->jtpt[sortIter];
    tempPhi = inJt->jtphi[sortIter];
    tempEta = inJt->jteta[sortIter];
    tempTrkMax = inJt->trackMax[sortIter];
    tempRawPt = inJt->rawpt[sortIter];
    tempRefPt = inJt->refpt[sortIter];
    tempRefPhi = inJt->refphi[sortIter];
    tempRefEta = inJt->refeta[sortIter];
    tempRefPart = inJt->refparton_flavor[sortIter];

    inJt->jtpt[sortIter] = inJt->jtpt[maxIndex];
    inJt->jtphi[sortIter] = inJt->jtphi[maxIndex];
    inJt->jteta[sortIter] = inJt->jteta[maxIndex];
    inJt->trackMax[sortIter] = inJt->trackMax[maxIndex];
    inJt->rawpt[sortIter] = inJt->rawpt[maxIndex];
    inJt->refpt[sortIter] = inJt->refpt[maxIndex];
    inJt->refphi[sortIter] = inJt->refphi[maxIndex];
    inJt->refeta[sortIter] = inJt->refeta[maxIndex];
    inJt->refparton_flavor[sortIter] = inJt->refparton_flavor[maxIndex];

    inJt->jtpt[maxIndex] = tempPt;
    inJt->jtphi[maxIndex] = tempPhi;
    inJt->jteta[maxIndex] = tempEta;
    inJt->trackMax[maxIndex] = tempTrkMax;
    inJt->rawpt[maxIndex] = tempRawPt;
    inJt->refpt[maxIndex] = tempRefPt;
    inJt->refphi[maxIndex] = tempRefPhi;
    inJt->refeta[maxIndex] = tempRefEta;
    inJt->refparton_flavor[maxIndex] = tempRefPart;

    sortIter++;
  }

  return;
}


Bool_t passesDijet(Jets jtCollection, Int_t algNum, Int_t &lPtCut, Int_t &sLPtCut)
{
  Int_t leadJtIndex = -1;
  Int_t subLeadJtIndex = -1;

  if(algNum < 11 || algNum > 14){
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
	  if(TMath::Abs(jtCollection.jteta[jtEntry]) < jtEtaCut) leadJtIndex = jtEntry;
	}
	else{
	  lPtCut++;
	  return false;
	}
      }
      else if(subLeadJtIndex < 0){
	if(jtCollection.jtpt[jtEntry] > subLeadJtPtCut){
	  if(TMath::Abs(jtCollection.jteta[jtEntry]) < jtEtaCut) subLeadJtIndex = jtEntry;
	}
	else{
	  sLPtCut++;
	  return false;
	}
      }
      else return true;
    }
    
    if(leadJtIndex >= 0){
      if(subLeadJtIndex >= 0) return true;
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
  else{
    if(jtCollection.ngen == 0){
      lPtCut++;
      return false;
    }
    else if(jtCollection.ngen == 1){
      sLPtCut++;
      return false;
    }

    for(Int_t jtEntry = 0; jtEntry < jtCollection.ngen; jtEntry++){
      if(leadJtIndex < 0){
	if(jtCollection.genpt[jtEntry] > leadJtPtCut){
	  if(TMath::Abs(jtCollection.geneta[jtEntry]) < jtEtaCut) leadJtIndex = jtEntry;
	}
	else{
	  lPtCut++;
	  return false;
	}
      }
      else if(subLeadJtIndex < 0){
	if(jtCollection.genpt[jtEntry] > subLeadJtPtCut){
	  if(TMath::Abs(jtCollection.geneta[jtEntry]) < jtEtaCut) subLeadJtIndex = jtEntry;
	}
	else{
	  sLPtCut++;
	  return false;
	}
      }
      else return true;
    }
    
    if(leadJtIndex >= 0){
      if(subLeadJtIndex >= 0) return true;
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



int makeDiJetIniSkim(string fList = "", sampleType sType = kHIDATA, const char *outName = "defaultName_DIJETINISKIM.root", Int_t num = 0, Bool_t justJt = false, Bool_t oldSample = false)
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

  HiForest *c = new HiForest(listOfFiles[num].data(), "Forest", cType, montecarlo);
  std::cout << "JobFile: " << listOfFiles[num] << std::endl;

  c->skimTree->SetBranchStatus("*", 0);
  c->skimTree->SetBranchStatus("pHBHENoiseFilter", 1);
  
  c->evtTree->SetBranchStatus("*", 0);
  c->evtTree->SetBranchStatus("vz", 1);
  if(hi) c->evtTree->SetBranchStatus("hiEvtPlanes", 1);
  c->evtTree->SetBranchStatus("run", 1);
  c->evtTree->SetBranchStatus("evt", 1);
  c->evtTree->SetBranchStatus("lumi", 1);  

  c->trackTree->SetBranchStatus("*", 0);
  c->trackTree->SetBranchStatus("nTrk", 1);
  c->trackTree->SetBranchStatus("trkPt", 1);
  c->trackTree->SetBranchStatus("trkPhi", 1);
  c->trackTree->SetBranchStatus("trkEta", 1);
  c->trackTree->SetBranchStatus("highPurity", 1);
  c->trackTree->SetBranchStatus("trkDz1", 1);
  c->trackTree->SetBranchStatus("trkDzError1", 1);
  c->trackTree->SetBranchStatus("trkDxy1", 1);
  c->trackTree->SetBranchStatus("trkDxyError1", 1);
  c->trackTree->SetBranchStatus("trkPtError", 1);

  if(montecarlo){
    c->genParticleTree->SetBranchStatus("*", 0);
    c->genParticleTree->SetBranchStatus("b", 1);
    c->genParticleTree->SetBranchStatus("mult", 1);
    c->genParticleTree->SetBranchStatus("chg", 1);
    c->genParticleTree->SetBranchStatus("pt", 1);
    c->genParticleTree->SetBranchStatus("phi", 1);
    c->genParticleTree->SetBranchStatus("eta", 1);
    c->genParticleTree->SetBranchStatus("sube", 1);
  }

  if(hi){
    c->skimTree->SetBranchStatus("pcollisionEventSelection", 1);

    c->evtTree->SetBranchStatus("hiBin", 1);

    setJtBranches(c->akPu3CaloJetTree, montecarlo);
    setJtBranches(c->akPu4CaloJetTree, montecarlo);
    setJtBranches(c->akPu5CaloJetTree, montecarlo);
    if(!oldSample) setJtBranches(c->akVs2CaloJetTree, montecarlo, true);
    setJtBranches(c->akVs3CaloJetTree, montecarlo, true);
    setJtBranches(c->akVs4CaloJetTree, montecarlo, true);
    setJtBranches(c->akVs5CaloJetTree, montecarlo, true);
    setJtBranches(c->akPu3PFJetTree, montecarlo);
    setJtBranches(c->akVs3PFJetTree, montecarlo);
  }
  else{
    c->skimTree->SetBranchStatus("pPAcollisionEventSelectionPA", 1);

    if(!oldSample) setJtBranches(c->akPu3CaloJetTree, montecarlo);
    if(!oldSample) setJtBranches(c->akPu4CaloJetTree, montecarlo);
    if(!oldSample) setJtBranches(c->akPu5CaloJetTree, montecarlo);
    if(!oldSample) setJtBranches(c->ak2CaloJetTree, montecarlo, true);
    setJtBranches(c->ak3CaloJetTree, montecarlo, true);
    if(!oldSample) setJtBranches(c->ak4CaloJetTree, montecarlo, true);
    if(!oldSample) setJtBranches(c->ak5CaloJetTree, montecarlo, true);
    if(!oldSample) setJtBranches(c->akPu3PFJetTree, montecarlo);
  }

  c->LoadNoTrees();

  c->hasSkimTree = true;
  c->hasEvtTree = true;

  if(hi){
    c->hasAkPu3JetTree = true;
    c->hasAkVs3PFJetTree = true;
    c->hasAkPu3CaloJetTree = true;
    c->hasAkPu4CaloJetTree = true;
    c->hasAkPu5CaloJetTree = true;
    if(!oldSample) c->hasAkVs2CaloJetTree = true;
    c->hasAkVs3CaloJetTree = true;
    c->hasAkVs4CaloJetTree = true;
    c->hasAkVs5CaloJetTree = true;
  }
  else{
    if(!oldSample) c->hasAk2CaloJetTree = true;
    c->hasAk3CaloJetTree = true;
    if(!oldSample) c->hasAk4CaloJetTree = true;
    if(!oldSample) c->hasAk5CaloJetTree = true;
    if(!oldSample) c->hasAkPu3CaloJetTree = true;
    if(!oldSample) c->hasAkPu4CaloJetTree = true;
    if(!oldSample) c->hasAkPu5CaloJetTree = true;
    if(!oldSample) c->hasAkPu3JetTree = true;
  }

  c->hasTrackTree = true;
  c->hasPFTree = true;
  if(montecarlo) c->hasGenParticleTree = true;

  c->InitTree();

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

  Int_t nEvtsOutTag = 0;

  if(!hi)
    nentries = c->ak3CaloJetTree->GetEntries();
  else
    nentries = c->GetEntries();

  std::cout << nentries << std::endl;

  TFile *outFile = new TFile(Form("%s_%d_%d.root", outName, num, 0), "RECREATE");
  InitDiJetIniSkim(sType, justJt);

  InitFactCorrFiles(sType);
  InitFactCorrHists(sType);

  InitJECCorrFiles(sType);
  InitJECCorrHists(sType);

  Int_t totEv = 0;
  Int_t selectCut = 0;
  Int_t vzCut = 0;

  Int_t AlgLeadJtPtCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Int_t AlgSubLeadJtPtCut[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  for(Long64_t jentry = 0; jentry < nentries; jentry++){
    c->hasTrackTree = false;
    if(montecarlo) c->hasGenParticleTree = false;

    c->GetEntry(jentry);

    totEv++;

    if(jentry%1000 == 0)
      std::cout << jentry << ", " << nEvtsOutTag << std::endl;

    if(!c->selectEvent()){
      selectCut++;
      continue;
    }

    if(TMath::Abs(c->evt.vz - meanVz) > 15){
      vzCut++;
      continue;
    }

    //particle flow

    Jets AlgJtCollection[17];
    Int_t algMax = 10;

    if(hi){
      AlgJtCollection[0] = c->akPu3Calo;
      AlgJtCollection[1] = c->akPu4Calo;
      AlgJtCollection[2] = c->akPu5Calo;
      if(!oldSample) AlgJtCollection[3] = c->akVs2Calo;
      AlgJtCollection[4] = c->akVs3Calo;
      AlgJtCollection[5] = c->akVs4Calo;
      AlgJtCollection[6] = c->akVs5Calo;
      if(!oldSample){
	AlgJtCollection[7] = c->akVs2Calo;
        getCorrJtCollection(7, &(AlgJtCollection[7]), &(c->pf), sType, c->evt.hiBin);
      }
      AlgJtCollection[8] = c->akVs3Calo;
      getCorrJtCollection(8, &(AlgJtCollection[8]), &(c->pf), sType, c->evt.hiBin);
      AlgJtCollection[9] = c->akVs4Calo;
      getCorrJtCollection(9, &(AlgJtCollection[9]), &(c->pf), sType, c->evt.hiBin);
      AlgJtCollection[10] = c->akVs5Calo;
      getCorrJtCollection(10, &(AlgJtCollection[10]), &(c->pf), sType, c->evt.hiBin);
      if(!oldSample) AlgJtCollection[11] = c->akVs2Calo;
      AlgJtCollection[12] = c->akVs3Calo;
      AlgJtCollection[13] = c->akVs4Calo;
      AlgJtCollection[14] = c->akVs5Calo;
      AlgJtCollection[15] = c->akPu3PF;
      AlgJtCollection[16] = c->akVs3PF;
    }
    else{
      if(!oldSample) AlgJtCollection[0] = c->akPu3Calo;
      if(!oldSample) AlgJtCollection[1] = c->akPu4Calo;
      if(!oldSample) AlgJtCollection[2] = c->akPu5Calo;
      if(!oldSample) AlgJtCollection[3] = c->ak2Calo;
      AlgJtCollection[4] = c->ak3Calo;
      if(!oldSample) AlgJtCollection[5] = c->ak4Calo;
      if(!oldSample) AlgJtCollection[6] = c->ak5Calo;
      if(!oldSample){
	AlgJtCollection[7] = c->ak2Calo;
	getCorrJtCollection(7, &(AlgJtCollection[7]), &(c->pf), sType, c->evt.hiBin);
      }
      AlgJtCollection[8] = c->ak3Calo;
      getCorrJtCollection(8, &(AlgJtCollection[8]), &(c->pf), sType, c->evt.hiBin);
      if(!oldSample){
	AlgJtCollection[9] = c->ak4Calo;
	getCorrJtCollection(9, &(AlgJtCollection[9]), &(c->pf), sType, c->evt.hiBin);
      }
      if(!oldSample){
	AlgJtCollection[10] = c->ak5Calo;
	getCorrJtCollection(10, &(AlgJtCollection[10]), &(c->pf), sType, c->evt.hiBin);
      }
      if(!oldSample) AlgJtCollection[11] = c->ak2Calo;
      AlgJtCollection[12] = c->ak3Calo;
      if(!oldSample) AlgJtCollection[13] = c->ak4Calo;
      if(!oldSample) AlgJtCollection[14] = c->ak5Calo;
      algMax = 15;
    }

    Bool_t algPasses[17] = {false, false, false, false, false, false, false, false, false, false, false, false, false, false};

    for(Int_t algIter = 0; algIter < algMax; algIter++){
      if(oldSample && algIter == 3) continue;
      else if(oldSample && !hi && algIter != 4) continue;

      algPasses[algIter] = passesDijet(AlgJtCollection[algIter], algIter, AlgLeadJtPtCut[algIter], AlgSubLeadJtPtCut[algIter]);
    }

    //truth, doesn't work w/ getLeadJt because truth doesnt get its own tree

    Bool_t passOneAlg = false;
    for(Int_t algIter = 0; algIter < 17; algIter++){
      if(algPasses[algIter]){
	passOneAlg = true;
	break;
      }
    }
    if(!passOneAlg) continue;

    c->hasTrackTree = true;
    if(montecarlo) c->hasGenParticleTree = true;

    c->GetEntry(jentry);

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

    if(hi) hiBinIni_ = c->evt.hiBin;

    InitPosArrPbPb(hiBinIni_);

    //Iterate over jets

    nPu3Calo_ = 0;
    nPu4Calo_ = 0;
    nPu5Calo_ = 0;
    nVs2Calo_ = 0;
    nVs3Calo_ = 0;
    nVs4Calo_ = 0;
    nVs5Calo_ = 0;
    nVs2CaloCorr_ = 0;
    nVs3CaloCorr_ = 0;
    nVs4CaloCorr_ = 0;
    nVs5CaloCorr_ = 0;
    nT2_ = 0;
    nT3_ = 0;
    nT4_ = 0;
    nT5_ = 0;
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



    for(Int_t Vs2CaloCorrIter = 0; Vs2CaloCorrIter < AlgJtCollection[7].nref; Vs2CaloCorrIter++){
      if(AlgJtCollection[7].jtpt[Vs2CaloCorrIter] < 30.0)
	break;
      else if(TMath::Abs(AlgJtCollection[7].jteta[Vs2CaloCorrIter]) > jtEtaCut)
	continue;
      
      Vs2CaloCorrPt_[nVs2CaloCorr_] = AlgJtCollection[7].jtpt[Vs2CaloCorrIter];
      Vs2CaloCorrPhi_[nVs2CaloCorr_] = AlgJtCollection[7].jtphi[Vs2CaloCorrIter];
      Vs2CaloCorrEta_[nVs2CaloCorr_] = AlgJtCollection[7].jteta[Vs2CaloCorrIter];
      
      Vs2CaloCorrTrkMax_[nVs2CaloCorr_] = AlgJtCollection[7].trackMax[Vs2CaloCorrIter];
      Vs2CaloCorrRawPt_[nVs2CaloCorr_] = AlgJtCollection[7].rawpt[Vs2CaloCorrIter];
      
      if(montecarlo){
	Vs2CaloCorrRefPt_[nVs2CaloCorr_] = AlgJtCollection[7].refpt[Vs2CaloCorrIter];
	Vs2CaloCorrRefPhi_[nVs2CaloCorr_] = AlgJtCollection[7].refphi[Vs2CaloCorrIter];
	Vs2CaloCorrRefEta_[nVs2CaloCorr_] = AlgJtCollection[7].refeta[Vs2CaloCorrIter];
	Vs2CaloCorrRefPart_[nVs2CaloCorr_] = AlgJtCollection[7].refparton_flavor[Vs2CaloCorrIter];
      }
      
      nVs2CaloCorr_++;
    }
       
    for(Int_t Vs3CaloCorrIter = 0; Vs3CaloCorrIter < AlgJtCollection[8].nref; Vs3CaloCorrIter++){
      if(AlgJtCollection[8].jtpt[Vs3CaloCorrIter] < 30.0)
	break;
      else if(TMath::Abs(AlgJtCollection[8].jteta[Vs3CaloCorrIter]) > jtEtaCut)
	continue;
      
      Vs3CaloCorrPt_[nVs3CaloCorr_] = AlgJtCollection[8].jtpt[Vs3CaloCorrIter];
      Vs3CaloCorrPhi_[nVs3CaloCorr_] = AlgJtCollection[8].jtphi[Vs3CaloCorrIter];
      Vs3CaloCorrEta_[nVs3CaloCorr_] = AlgJtCollection[8].jteta[Vs3CaloCorrIter];
      
      Vs3CaloCorrTrkMax_[nVs3CaloCorr_] = AlgJtCollection[8].trackMax[Vs3CaloCorrIter];
      Vs3CaloCorrRawPt_[nVs3CaloCorr_] = AlgJtCollection[8].rawpt[Vs3CaloCorrIter];

      if(montecarlo){
	Vs3CaloCorrRefPt_[nVs3CaloCorr_] = AlgJtCollection[8].refpt[Vs3CaloCorrIter];
	Vs3CaloCorrRefPhi_[nVs3CaloCorr_] = AlgJtCollection[8].refphi[Vs3CaloCorrIter];
	Vs3CaloCorrRefEta_[nVs3CaloCorr_] = AlgJtCollection[8].refeta[Vs3CaloCorrIter];
	Vs3CaloCorrRefPart_[nVs3CaloCorr_] = AlgJtCollection[8].refparton_flavor[Vs3CaloCorrIter];
      }
      
      nVs3CaloCorr_++;
    }

    for(Int_t Vs4CaloCorrIter = 0; Vs4CaloCorrIter < AlgJtCollection[9].nref; Vs4CaloCorrIter++){
      if(AlgJtCollection[9].jtpt[Vs4CaloCorrIter] < 30.0)
	break;
      else if(TMath::Abs(AlgJtCollection[9].jteta[Vs4CaloCorrIter]) > jtEtaCut)
	continue;

      Vs4CaloCorrPt_[nVs4CaloCorr_] = AlgJtCollection[9].jtpt[Vs4CaloCorrIter];
      Vs4CaloCorrPhi_[nVs4CaloCorr_] = AlgJtCollection[9].jtphi[Vs4CaloCorrIter];
      Vs4CaloCorrEta_[nVs4CaloCorr_] = AlgJtCollection[9].jteta[Vs4CaloCorrIter];

      Vs4CaloCorrTrkMax_[nVs4CaloCorr_] = AlgJtCollection[9].trackMax[Vs4CaloCorrIter];
      Vs4CaloCorrRawPt_[nVs4CaloCorr_] = AlgJtCollection[9].rawpt[Vs4CaloCorrIter];

      if(montecarlo){
	Vs4CaloCorrRefPt_[nVs4CaloCorr_] = AlgJtCollection[9].refpt[Vs4CaloCorrIter];
	Vs4CaloCorrRefPhi_[nVs4CaloCorr_] = AlgJtCollection[9].refphi[Vs4CaloCorrIter];
	Vs4CaloCorrRefEta_[nVs4CaloCorr_] = AlgJtCollection[9].refeta[Vs4CaloCorrIter];
	Vs4CaloCorrRefPart_[nVs4CaloCorr_] = AlgJtCollection[9].refparton_flavor[Vs4CaloCorrIter];
      }

      nVs4CaloCorr_++;
    }

    for(Int_t Vs5CaloCorrIter = 0; Vs5CaloCorrIter < AlgJtCollection[10].nref; Vs5CaloCorrIter++){
      if(AlgJtCollection[10].jtpt[Vs5CaloCorrIter] < 30.0)
	break;
      else if(TMath::Abs(AlgJtCollection[10].jteta[Vs5CaloCorrIter]) > jtEtaCut)
	continue;

      Vs5CaloCorrPt_[nVs5CaloCorr_] = AlgJtCollection[10].jtpt[Vs5CaloCorrIter];
      Vs5CaloCorrPhi_[nVs5CaloCorr_] = AlgJtCollection[10].jtphi[Vs5CaloCorrIter];
      Vs5CaloCorrEta_[nVs5CaloCorr_] = AlgJtCollection[10].jteta[Vs5CaloCorrIter];

      Vs5CaloCorrTrkMax_[nVs5CaloCorr_] = AlgJtCollection[10].trackMax[Vs5CaloCorrIter];
      Vs5CaloCorrRawPt_[nVs5CaloCorr_] = AlgJtCollection[10].rawpt[Vs5CaloCorrIter];

      if(montecarlo){
	Vs5CaloCorrRefPt_[nVs5CaloCorr_] = AlgJtCollection[10].refpt[Vs5CaloCorrIter];
	Vs5CaloCorrRefPhi_[nVs5CaloCorr_] = AlgJtCollection[10].refphi[Vs5CaloCorrIter];
	Vs5CaloCorrRefEta_[nVs5CaloCorr_] = AlgJtCollection[10].refeta[Vs5CaloCorrIter];
	Vs5CaloCorrRefPart_[nVs5CaloCorr_] = AlgJtCollection[10].refparton_flavor[Vs5CaloCorrIter];
      }

      nVs5CaloCorr_++;
    }


    if(montecarlo){
      for(Int_t T2Iter = 0; T2Iter < AlgJtCollection[11].ngen; T2Iter++){
	if(AlgJtCollection[11].genpt[T2Iter] < 30.0)
	  break;
	else if(TMath::Abs(AlgJtCollection[11].geneta[T2Iter]) > jtEtaCut)
	  continue;
	
	T2Pt_[nT2_] = AlgJtCollection[11].genpt[T2Iter];
	T2Phi_[nT2_] = AlgJtCollection[11].genphi[T2Iter];
	T2Eta_[nT2_] = AlgJtCollection[11].geneta[T2Iter];
	T2Part_[nT2_] = AlgJtCollection[11].refparton_flavor[AlgJtCollection[11].genmatchindex[T2Iter]];
	
	nT2_++;
      }

      for(Int_t T3Iter = 0; T3Iter < AlgJtCollection[12].ngen; T3Iter++){
	if(AlgJtCollection[12].genpt[T3Iter] < 30.0)
	  break;
	else if(TMath::Abs(AlgJtCollection[12].geneta[T3Iter]) > jtEtaCut)
	  continue;
	
	T3Pt_[nT3_] = AlgJtCollection[12].genpt[T3Iter];
	T3Phi_[nT3_] = AlgJtCollection[12].genphi[T3Iter];
	T3Eta_[nT3_] = AlgJtCollection[12].geneta[T3Iter];
	T3Part_[nT3_] = AlgJtCollection[12].refparton_flavor[AlgJtCollection[12].genmatchindex[T3Iter]];
	
	nT3_++;
      }

      for(Int_t T4Iter = 0; T4Iter < AlgJtCollection[13].ngen; T4Iter++){
	if(AlgJtCollection[13].genpt[T4Iter] < 30.0)
	  break;
	else if(TMath::Abs(AlgJtCollection[13].geneta[T4Iter]) > jtEtaCut)
	  continue;
	
	T4Pt_[nT4_] = AlgJtCollection[13].genpt[T4Iter];
	T4Phi_[nT4_] = AlgJtCollection[13].genphi[T4Iter];
	T4Eta_[nT4_] = AlgJtCollection[13].geneta[T4Iter];
	T4Part_[nT4_] = AlgJtCollection[13].refparton_flavor[AlgJtCollection[13].genmatchindex[T4Iter]];
	
	nT4_++;
      }

      for(Int_t T5Iter = 0; T5Iter < AlgJtCollection[14].ngen; T5Iter++){
	if(AlgJtCollection[14].genpt[T5Iter] < 30.0)
	  break;
	else if(TMath::Abs(AlgJtCollection[14].geneta[T5Iter]) > jtEtaCut)
	  continue;
	
	T5Pt_[nT5_] = AlgJtCollection[14].genpt[T5Iter];
	T5Phi_[nT5_] = AlgJtCollection[14].genphi[T5Iter];
	T5Eta_[nT5_] = AlgJtCollection[14].geneta[T5Iter];
	T5Part_[nT5_] = AlgJtCollection[14].refparton_flavor[AlgJtCollection[14].genmatchindex[T5Iter]];
	
	nT5_++;
      }
    }
    
    if(hi){
      for(Int_t Pu3PFIter = 0; Pu3PFIter < AlgJtCollection[15].nref; Pu3PFIter++){
	if(AlgJtCollection[15].jtpt[Pu3PFIter] < 30.0)
	  break;
	else if(TMath::Abs(AlgJtCollection[15].jteta[Pu3PFIter]) > jtEtaCut)
	  continue;
	
	Pu3PFPt_[nPu3PF_] = AlgJtCollection[15].jtpt[Pu3PFIter];
	Pu3PFPhi_[nPu3PF_] = AlgJtCollection[15].jtphi[Pu3PFIter];
	Pu3PFEta_[nPu3PF_] = AlgJtCollection[15].jteta[Pu3PFIter];
	
	Pu3PFTrkMax_[nPu3PF_] = AlgJtCollection[15].trackMax[Pu3PFIter];
	Pu3PFRawPt_[nPu3PF_] = AlgJtCollection[15].rawpt[Pu3PFIter];
	
	if(montecarlo){
	  Pu3PFRefPt_[nPu3PF_] = AlgJtCollection[15].refpt[Pu3PFIter];
	  Pu3PFRefPhi_[nPu3PF_] = AlgJtCollection[15].refphi[Pu3PFIter];
	  Pu3PFRefEta_[nPu3PF_] = AlgJtCollection[15].refeta[Pu3PFIter];
	  Pu3PFRefPart_[nPu3PF_] = AlgJtCollection[15].refparton_flavor[Pu3PFIter];
	}
	
	nPu3PF_++;
      }
      
      for(Int_t Vs3PFIter = 0; Vs3PFIter < AlgJtCollection[16].nref; Vs3PFIter++){
	if(AlgJtCollection[16].jtpt[Vs3PFIter] < 30.0)
	  break;
	else if(TMath::Abs(AlgJtCollection[16].jteta[Vs3PFIter]) > jtEtaCut)
	  continue;
	
	Vs3PFPt_[nVs3PF_] = AlgJtCollection[16].jtpt[Vs3PFIter];
	Vs3PFPhi_[nVs3PF_] = AlgJtCollection[16].jtphi[Vs3PFIter];
	Vs3PFEta_[nVs3PF_] = AlgJtCollection[16].jteta[Vs3PFIter];
	
	Vs3PFTrkMax_[nVs3PF_] = AlgJtCollection[16].trackMax[Vs3PFIter];
	Vs3PFRawPt_[nVs3PF_] = AlgJtCollection[16].rawpt[Vs3PFIter];

	if(montecarlo){
	  Vs3PFRefPt_[nVs3PF_] = AlgJtCollection[16].refpt[Vs3PFIter];
	  Vs3PFRefPhi_[nVs3PF_] = AlgJtCollection[16].refphi[Vs3PFIter];
	  Vs3PFRefEta_[nVs3PF_] = AlgJtCollection[16].refeta[Vs3PFIter];
	  Vs3PFRefPart_[nVs3PF_] = AlgJtCollection[16].refparton_flavor[Vs3PFIter];
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
      
      if(TMath::Abs(trkCollection.trkEta[trkEntry]) > 2.4) continue;
      
      if(trkCollection.trkPt[trkEntry] <= 0.5) continue;
	
      if(!trkCollection.highPurity[trkEntry]) continue;
	
      if(TMath::Abs(trkCollection.trkDz1[trkEntry]/trkCollection.trkDzError1[trkEntry]) > 3) continue;
	
      if(TMath::Abs(trkCollection.trkDxy1[trkEntry]/trkCollection.trkDxyError1[trkEntry]) > 3) continue;
	
      if(trkCollection.trkPtError[trkEntry]/trkCollection.trkPt[trkEntry] > 0.1) continue;
	
      trkPt_[nTrk_] = trkCollection.trkPt[trkEntry];
      trkPhi_[nTrk_] = trkCollection.trkPhi[trkEntry];
      trkEta_[nTrk_] = trkCollection.trkEta[trkEntry];
      tempTL.SetPtEtaPhiM(trkCollection.trkPt[trkEntry], trkCollection.trkEta[trkEntry], trkCollection.trkPhi[trkEntry], 0);
      jtVect_p->push_back(tempTL);
	
      Float_t tempRMin = getTrkRMin(trkPhi_[nTrk_], trkEta_[nTrk_], nVs3Calo_, Vs3CaloPt_, Vs3CaloPhi_, Vs3CaloEta_);
      trkCorr_[nTrk_] = factorizedPtCorr(getPtBin(trkPt_[nTrk_], sType), hiBinIni_, trkPt_[nTrk_], trkPhi_[nTrk_], trkEta_[nTrk_], tempRMin, sType);

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
    if(hi) InitPosArrPbPb(hiBinIni_);

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

    Int_t nEvtsPassed = jetTreeIni_p->GetEntries();

    if(nEvtsPassed%nEvtPerFile == 0){
      outFile->cd();

      nEvtsOutTag++;

      jetTreeIni_p->Write("", TObject::kOverwrite);
      if(!justJt){
	trackTreeIni_p->Write("", TObject::kOverwrite);
	pfCandTreeIni_p->Write("", TObject::kOverwrite);

	if(montecarlo) genTreeIni_p->Write("", TObject::kOverwrite);
      }

      CleanupDiJetIniSkim();
      outFile->Close();
      delete outFile;
      outFile = new TFile(Form("%s_%d_%d.root", outName, num, nEvtsOutTag), "RECREATE");
      InitDiJetIniSkim(sType, justJt);
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

  CleanupDiJetIniSkim();
  outFile->Close();
  delete outFile;

  delete c;

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
  if(argc != 7)
    {
      std::cout << "Usage: makeDiJetIniSkim <inputFile> <sType> <outputFile> <#> <justJtBool> <oldSample>" << std::endl;
      std::cout << argc << std::endl;
      std::cout << argv[0] << ", "  << argv[1] << ", " << argv[2] << ", " << argv[3]<< ", " << argv[4] << ", " << argv[5] << ", " << argv[6] << std::endl;
      return 1;
    }

  int rStatus = -1;

  rStatus = makeDiJetIniSkim(argv[1], sampleType(atoi(argv[2])), argv[3], atoi(argv[4]), Bool_t(atoi(argv[5])), Bool_t(atoi(argv[6])));

  return rStatus;
}
