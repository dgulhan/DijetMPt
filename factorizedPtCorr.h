//=============================================                                                         // Author: Chris McGinn                                                                                 //                                                                                                      // DiJet Analysis Class (MC)                                                                            //                                                                                                     
//=============================================                                                                         

#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"

#include <iostream>

enum sample{
  kHIDATA, //0
  kHIMC,   //1
  kPPDATA, //2
  kPPMC,   //3
  kPADATA, //4
  kPAMC    //5
}

//Current # of correction histograms

const Int_t nHistPbPb = 29;
const Int_t nHistPP = 4;

//Vs Calo File and Hist array

TFile* VsCaloFile_p[nHistPbPb];
TProfile* VsCalocent_p[nHistPbPb];
TProfile2D* VsCalophiEta_p[nHistPbPb];
TProfile* VsCalopt_p[nHistPbPb];
TProfile* VsCalodelR_p[nHistPbPb];


//Fake Vs Calo File and Hist array

TFile* FakeVsCaloFile_p[nHistPbPb];
TProfile* FakeVsCalocent_p[nHistPbPb];
TProfile2D* FakeVsCalophiEta_p[nHistPbPb];
TProfile* FakeVsCalopt_p[nHistPbPb];
TProfile* FakeVsCalodelR_p[nHistPbPb];

void InitCorrFiles(sample sType = kHIDATA)
{
  //File names w/ various binnings, ordered by pt and then centrality. Each Jet Algorithm gets a file array

  if(sType == kHIDATA || sType == kHIMC){
    VsCaloFile_p[0] = new TFile("eff_pt50_55_cent0_10.root", "READ");
    VsCaloFile_p[1] = new TFile("eff_pt50_55_cent10_20.root", "READ");
    VsCaloFile_p[2] = new TFile("eff_pt50_55_cent20_30.root", "READ");
    VsCaloFile_p[3] = new TFile("eff_pt50_55_cent30_50.root", "READ");
    VsCaloFile_p[4] = new TFile("eff_pt50_55_cent50_100.root", "READ");
    
    VsCaloFile_p[5] = new TFile("eff_pt55_65_cent0_10.root", "READ");
    VsCaloFile_p[6] = new TFile("eff_pt55_65_cent10_20.root", "READ");
    VsCaloFile_p[7] = new TFile("eff_pt55_65_cent20_30.root", "READ");
    VsCaloFile_p[8] = new TFile("eff_pt55_65_cent30_50.root", "READ");
    VsCaloFile_p[9] = new TFile("eff_pt55_65_cent50_100.root", "READ");
    
    VsCaloFile_p[10] = new TFile("eff_pt65_80_cent0_10.root", "READ");
    VsCaloFile_p[11] = new TFile("eff_pt65_80_cent10_20.root", "READ");
    VsCaloFile_p[12] = new TFile("eff_pt65_80_cent20_30.root", "READ");
    VsCaloFile_p[13] = new TFile("eff_pt65_80_cent30_50.root", "READ");
    VsCaloFile_p[14] = new TFile("eff_pt65_80_cent50_100.root", "READ");
    
    VsCaloFile_p[15] = new TFile("eff_pt80_100_cent0_10.root", "READ");
    VsCaloFile_p[16] = new TFile("eff_pt80_100_cent10_20.root", "READ");
    VsCaloFile_p[17] = new TFile("eff_pt80_100_cent20_30.root", "READ");
    VsCaloFile_p[18] = new TFile("eff_pt80_100_cent30_50.root", "READ");
    VsCaloFile_p[19] = new TFile("eff_pt80_100_cent50_100.root", "READ");
    
    VsCaloFile_p[20] = new TFile("eff_pt100_300_cent0_10.root", "READ");
    VsCaloFile_p[21] = new TFile("eff_pt100_300_cent10_20.root", "READ");
    VsCaloFile_p[22] = new TFile("eff_pt100_300_cent20_30.root", "READ");
    VsCaloFile_p[23] = new TFile("eff_pt100_300_cent30_50.root", "READ");
    VsCaloFile_p[24] = new TFile("eff_pt100_300_cent50_100.root", "READ");

    VsCaloFile_p[25] = new TFile("eff_pt300_800_cent0_10.root", "READ");
    VsCaloFile_p[26] = new TFile("eff_pt300_800_cent10_20.root", "READ");
    VsCaloFile_p[27] = new TFile("eff_pt300_800_cent20_100.root", "READ");
    VsCaloFile_p[28] = new TFile("eff_pt800_30000_cent0_100.root", "READ");
    
    //Fakes VsCalo
    
    FakeVsCaloFile_p[0] = new TFile("fake_pt50_55_cent0_10.root", "READ");
    FakeVsCaloFile_p[1] = new TFile("fake_pt50_55_cent10_20.root", "READ");
    FakeVsCaloFile_p[2] = new TFile("fake_pt50_55_cent20_30.root", "READ");
    FakeVsCaloFile_p[3] = new TFile("fake_pt50_55_cent30_50.root", "READ");
    FakeVsCaloFile_p[4] = new TFile("fake_pt50_55_cent50_100.root", "READ");
    
    FakeVsCaloFile_p[5] = new TFile("fake_pt55_65_cent0_10.root", "READ");
    FakeVsCaloFile_p[6] = new TFile("fake_pt55_65_cent10_20.root", "READ");
    FakeVsCaloFile_p[7] = new TFile("fake_pt55_65_cent20_30.root", "READ");
    FakeVsCaloFile_p[8] = new TFile("fake_pt55_65_cent30_50.root", "READ");
    FakeVsCaloFile_p[9] = new TFile("fake_pt55_65_cent50_100.root", "READ");
    
    FakeVsCaloFile_p[10] = new TFile("fake_pt65_80_cent0_10.root", "READ");
    FakeVsCaloFile_p[11] = new TFile("fake_pt65_80_cent10_20.root", "READ");
    FakeVsCaloFile_p[12] = new TFile("fake_pt65_80_cent20_30.root", "READ");
    FakeVsCaloFile_p[13] = new TFile("fake_pt65_80_cent30_50.root", "READ");
    FakeVsCaloFile_p[14] = new TFile("fake_pt65_80_cent50_100.root", "READ");
    
    FakeVsCaloFile_p[15] = new TFile("fake_pt80_100_cent0_10.root", "READ");
    FakeVsCaloFile_p[16] = new TFile("fake_pt80_100_cent10_20.root", "READ");
    FakeVsCaloFile_p[17] = new TFile("fake_pt80_100_cent20_30.root", "READ");
    FakeVsCaloFile_p[18] = new TFile("fake_pt80_100_cent30_50.root", "READ");
    FakeVsCaloFile_p[19] = new TFile("fake_pt80_100_cent50_100.root", "READ");
    
    FakeVsCaloFile_p[20] = new TFile("fake_pt100_300_cent0_10.root", "READ");
    FakeVsCaloFile_p[21] = new TFile("fake_pt100_300_cent10_20.root", "READ");
    FakeVsCaloFile_p[22] = new TFile("fake_pt100_300_cent20_30.root", "READ");
    FakeVsCaloFile_p[23] = new TFile("fake_pt100_300_cent30_50.root", "READ");
    FakeVsCaloFile_p[24] = new TFile("fake_pt100_300_cent50_100.root", "READ");
    
    FakeVsCaloFile_p[25] = new TFile("fake_pt300_800_cent0_10.root", "READ");
    FakeVsCaloFile_p[26] = new TFile("fake_pt300_800_cent10_20.root", "READ");
    FakeVsCaloFile_p[27] = new TFile("fake_pt300_800_cent20_100.root", "READ");
    FakeVsCaloFile_p[28] = new TFile("fake_pt800_30000_cent0_100.root", "READ");
  }
  else if(sType == kPPDATA || sType == kHIMC){
    return;
  }

  return;
}


void InitCorrHists(sample sType = kHIDATA)
{
  if(sType == kHIDATA || sType == kHIMC){
    for(Int_t hIter = 0; hIter < nHistPbPb; hIter++){
      VsCalocent_p[hIter] = (TProfile*)VsCaloFile_p[hIter]->Get("p_eff_cent");
      VsCalophiEta_p[hIter] = (TProfile2D*)VsCaloFile_p[hIter]->Get("p_eff_acceptance");
      VsCalopt_p[hIter] = (TProfile*)VsCaloFile_p[hIter]->Get("p_eff_pt");
      VsCalodelR_p[hIter] = (TProfile*)VsCaloFile_p[hIter]->Get("p_eff_rmin");
      
      FakeVsCalocent_p[hIter] = (TProfile*)FakeVsCaloFile_p[hIter]->Get("p_fake_cent");
      FakeVsCalophiEta_p[hIter] = (TProfile2D*)FakeVsCaloFile_p[hIter]->Get("p_fake_acceptance");
      FakeVsCalopt_p[hIter] = (TProfile*)FakeVsCaloFile_p[hIter]->Get("p_fake_pt");
      FakeVsCalodelR_p[hIter] = (TProfile*)FakeVsCaloFile_p[hIter]->Get("p_fake_rmin");
    }    
  }
  else if(sType == kPPDATA || sType == kPPMC){
    return;
  }

  return;
}


Int_t getPtBin(Float_t pt, Int_t hiSet1, Int_t hiSet2, Int_t hiSet3, Int_t hiSet4, Int_t hiSet5, Int_t hiSet6, Int_t hiSet7)
{
  Int_t ptPos = -1;

  if(.5 <= pt && pt < .55)
    ptPos = hiSet1;
  else if(.55 <= pt && pt < .65)
    ptPos = hiSet2;
  else if(.65 <= pt && pt < .80)
    ptPos = hiSet3;
  else if(.80 <= pt && pt < 1.00)
    ptPos = hiSet4;
  else if(1.00 <= pt && pt < 3.00)
    ptPos = hiSet5;
  else if(3.00 <= pt && pt < 8.00)
    ptPos = hiSet6;
  else if(8 <= pt)
    ptPos = hiSet7;

  return ptPos;
}


//Feed variables and the histograms along w/ the appropriate rmincut (currently rmin only defined to 3 for PuPF and 5 for calo)

Float_t factorizedPtCorr(Int_t hiBin, Float_t pt, Float_t phi, Float_t eta, Float_t rmin, TProfile* centProf_p, TProfile2D* etaPhiProf_p, TProfile* ptProf_p, TProfile* rminProf_p, Bool_t eff = true)
{
  Float_t corrFactor = 1;

  if(hiBin < 0 || hiBin > 200){
    if(eff)
      return 1;
    else
      return 0;
  }

  if(pt < .5){
    if(eff)
      return 1;
    else
      return 0;
  }

  if(eff){
    corrFactor = corrFactor*(centProf_p->GetBinContent(centProf_p->FindBin(hiBin)));
    corrFactor = corrFactor*(etaPhiProf_p->GetBinContent(etaPhiProf_p->FindBin(phi, eta)));
    corrFactor = corrFactor*(ptProf_p->GetBinContent(ptProf_p->FindBin(pt)));
    corrFactor = corrFactor*(rminProf_p->GetBinContent(rminProf_p->FindBin(rmin)));

    if(corrFactor == 0){
      if(pt > 100)
	corrFactor = .8;
      else
	corrFactor = 1;
    }
  }
  else{
    corrFactor = 0;
    corrFactor = corrFactor + centProf_p->GetBinContent(centProf_p->FindBin(hiBin));
    corrFactor = corrFactor + etaPhiProf_p->GetBinContent(etaPhiProf_p->FindBin(phi, eta));
    corrFactor = corrFactor + ptProf_p->GetBinContent(ptProf_p->FindBin(pt)) ;
    corrFactor = corrFactor +  rminProf_p->GetBinContent(rminProf_p->FindBin(rmin));
  }

  return corrFactor;
}



/*
Above implemented in code as follows:

 for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
      trkRMinPuPF_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[0]);
      trkRMinPuCalo_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[1]);
      trkRMinVsPF_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[2]);
      trkRMinVsCalo_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[3]);
    }

    if(montecarlo){
      for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
        trkRMinT_[trkEntry] = getTrkRMin(trkPhi_[trkEntry], trkEta_[trkEntry], AlgJtCollection[3], true);
      }
    }

    Int_t hiBinDiv[5] = {20, 40, 60, 100, 200};
    Int_t hiSetEff[15] = {0, 5, 10, 1, 6, 11, 2, 7, 12, 3, 8, 12, 4, 9, 12};

    for(Int_t hiBinIter = 0; hiBinIter < 5; hiBinIter++){
      if(hiBin_ < hiBinDiv[hiBinIter]){
        for(Int_t trkEntry = 0; trkEntry < nTrk_; trkEntry++){
          Int_t ptPos = getPtBin(trkPt_[trkEntry], hiSetEff[hiBinIter*3], hiSetEff[hiBinIter*3 + 1], hiSetEff[hiBinIter*3 + 2], 13);

          trkPtFactPuPF_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPuPF_[trkEntry], PuPFcent_p[ptPos], PuPFphiEta_p[ptPos], PuPFpt_p[ptPos], PuPFdelR_p[ptPos], 3);
          trkPtCorrPuPF_[trkEntry] = trkPt_[trkEntry]/trkPtFactPuPF_[trkEntry];

          trkPtFactPuCalo_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPuCalo_[trkEntry], PuCalocent_p[ptPos], PuCalophiEta_p[ptPos], PuCalopt_p[ptPos], PuCalodelR_p[ptPos], 5);
          trkPtCorrPuCalo_[trkEntry] = trkPt_[trkEntry]*(1 - factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinPuCalo_[trkEntry], FakePuCalocent_p[ptPos], FakePuCalophiEta_p[ptPos], FakePuCalopt_p[ptPos], FakePuCalodelR_p[ptPos], 5, false))/trkPtFactPuCalo_[trkEntry];

          trkPtFactVsCalo_[trkEntry] = factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinVsCalo_[trkEntry], VsCalocent_p[ptPos], VsCalophiEta_p[ptPos], VsCalopt_p[ptPos], VsCalodelR_p[ptPos], 5);
          trkPtCorrVsCalo_[trkEntry] = trkPt_[trkEntry]*(1 - factorizedPtCorr(hiBin_, trkPt_[trkEntry], trkPhi_[trkEntry], trkEta_[trkEntry], trkRMinVsCalo_[trkEntry], FakeVsCalocent_p[ptPos], FakeVsCalophiEta_p[ptPos], FakeVsCalopt_p[ptPos], FakeVsCalodelR_p[ptPos], 5, false))/trkPtFactVsCalo_[trkEntry];


	}
	break;
      }
    }
*/
