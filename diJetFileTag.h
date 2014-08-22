#ifndef diJetFileTag_h
#define diJetFileTag_h

#include "TMath.h"
#include "iostream"
#include "TH1F.h"
#include "TCut.h"

const char* fileTag;

const char* Di30a = "HydjetDrum_Pyquen_Dijet30_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMSKIM_20140423_0.root";

const char* Di50a = "HydjetDrum_Pyquen_Dijet50_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMSKIM_20140423_0.root";

const char* Di80a = "Pythia80_HydjetDrum_mix01_HiForest2_v20_CFMSKIM.root";
const char* Di80b = "Dijet80_HydjetDrum_v27_mergedV1_CFMSKIM.root";
const char* Di80c = "HydjetDrum_Pyquen_Dijet80_Embedded_d20140122_Track7_v2_CFMSKIM.root";
const char* Di80d = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0_CFMSKIM.root";
const char* Di80e = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet14_STARTHI53_LV1_merged_forest_0_CFMSKIM.root";
const char* Di80f = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_50k_CFMSKIM.root";
const char* Di80g = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_300k_CFMSKIM.root";
const char* Di80h = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0_300k_CFMSKIM.root";
const char* Di80i = "HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt_0_CFMSKIM.root";
const char* Di80j = "HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt_0_CFMSKIM_20140414_0.root";

const char* Di100a = "Dijet100_HydjetDrum_v27_mergedV1_CFMSKIM.root";
const char* Di100b = "HydjetDrum_Pyquen_Dijet100_FOREST_Track8_Jet24_FixedPtHat_v0_0_CFMSKIM_20140323_4_0.root";
const char* Di100c = "HydjetDrum_Pyquen_Dijet100_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMSKIM_20140423_0.root";

const char* Di120a = "HydjetDrum_Pyquen_Dijet120_Embedded_RECO_STARTHI53_LV1_Track8_Jet21_300k_v0_merged_0_CFMSKIM.root";
const char* Di120b = "HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHat_v0_0_CFMSKIM.root";
const char* Di120c = "HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMSKIM_20140423_0.root";

const char* DiAlla = "HydjetDrum_Pyquen_DijetAll_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMSKIM_20140425_0.root";

const char* DiAllc = "HydjetDrum_Pyquen_DijetAll_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMSKIM_20140508_0.root";

const char* EmDi80a = "PbPb_pythiaHYDJET_forest_EmEnrichedDijet80_CFMSKIM.root";

const char* DataA = "Track8_Jet17_GR_R_53_LV6_SUB_0_CFMSKIM.root";
const char* DataB = "hiForest_Jet80or95_GR_R_53_LV6_02Mar2014_1300CET_Track8_Jet15_0_1200k_CFMSKIM.root";
const char* DataC = "hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21_0_700k_CFMSKIM.root";
const char* DataD = "hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21_0_1200k_CFMSKIM.root";
const char* DataE = "hiForest_Jet80or95_GR_R_53_LV6_03Mar2014_1600CET_CMSSW_5_3_16_merged_0_CFMSKIM.root";
const char* DataF = "HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4-merged_0_CFMSKIM_20140430_0.root";
const char* DataG = "HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4ANDv9-merged_0_CFMSKIM_20140505_0.root";

const char* DataG_5 = "HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4ANDv9-merged_0_CFMSKIM_20140511_RAD5_0.root";

const char* DataF_JtCutDown = "HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4-merged_0_CFMSKIM_20140430_JtCutDown_0.root";
const char* DataG_JtCutDown = "HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4ANDv9-merged_0_CFMSKIM_20140508_jtCutDown_0.root";


const char* TESTMixing = "testJtANA_30CHECK_3rdJet_0.root";
const char* TESTPbPb = "HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4ANDv9-merged_0_DIJETANASKIM_20140619_CHECK_0.root";
const char* TESTPP = "HiForest_pp_Jet80_v8_PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_merged_forest_0_CFMANASKIM_20140612_0.root";

const char* TESTPbPbMC = "HydjetDrum_Pyquen_DijetMerge_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMANASKIM_20140618_0.root";
const char* TESTPPMC = "HiForest_pt80_PYTHIA_ppReco_JECv85_merged_forest_0_DIJETANASKIM_20140618_0.root";
const char* TESTPPMCGEN = "pt80_pp2013_P01_prod22_v81_merged_forest_0_DIJETANASKIM_20140703_0.root";

const char* name1[10] = {"0(10000, -10000, 10000)", "1(10000, -10000, 10000)", "2(10000, -10000, 10000)", "3(10000, -10000, 10000)", "4(10000, -10000, 10000)", "5(10000, -10000, 10000)", "6(10000, -10000, 10000)", "7(10000, -1000, 10000)", "8(10000, -1000, 10000)", "9(10000, -1000, 10000)"};
const char* name2[10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};

void setFileTag(const char* inName)
{
  if(!strcmp(inName, Di30a)){
    std::cout << Di30a << std::endl;
    fileTag = "Di30a";
  }
  else if(!strcmp(inName, Di50a)){
    std::cout << Di50a << std::endl;
    fileTag = "Di50a";
  }
  else if(!strcmp(inName, Di80a)){
    std::cout << Di80a << std::endl;
    fileTag = "Di80a";
  }
  else if(!strcmp(inName, Di80b)){
    std::cout << Di80b << std::endl;
    fileTag = "Di80b";
  }
  else if(!strcmp(inName, Di100a)){
    std::cout << Di100a << std::endl;
    fileTag = "Di100a";
  }
  else if(!strcmp(inName, Di100b)){
    std::cout << Di100b << std::endl;
    fileTag = "Di100b";
  }
  else if(!strcmp(inName, Di100c)){
    std::cout << Di100c << std::endl;
    fileTag = "Di100c";
  }
  else if(!strcmp(inName, Di120a)){
    std::cout << Di120a << std::endl;
    fileTag = "Di120a";
  }
  else if(!strcmp(inName, Di120b)){
    std::cout << Di120b << std::endl;
    fileTag = "Di120b";
  }
  else if(!strcmp(inName, Di120c)){
    std::cout << Di120c << std::endl;
    fileTag = "Di120c";
  }
  else if(!strcmp(inName, EmDi80a)){
    std::cout << EmDi80a << std::endl;
    fileTag = "EmDi80a";
  }
  else if(!strcmp(inName, Di80c)){
    std::cout << Di80c << std::endl;
    fileTag = "Di80c";
  }
  else if(!strcmp(inName, Di80d)){
    std::cout << Di80d << std::endl;
    fileTag = "Di80d";
  }
  else if(!strcmp(inName, Di80e)){
    std::cout << Di80e << std::endl;
    fileTag = "Di80e";
  }
  else if(!strcmp(inName, Di80f)){
    std::cout << Di80f << std::endl;
    fileTag = "Di80f";
  }
  else if(!strcmp(inName, Di80g)){
    std::cout << Di80g << std::endl;
    fileTag = "Di80g";
  }
  else if(!strcmp(inName, Di80h)){
    std::cout << Di80h << std::endl;
    fileTag = "Di80h";
  }
  else if(!strcmp(inName, Di80i)){
    std::cout << Di80i << std::endl;
    fileTag = "Di80i";
  }
  else if(!strcmp(inName, Di80j)){
    std::cout << Di80j << std::endl;
    fileTag = "Di80j";
  }
  else if(!strcmp(inName, DiAlla)){
    std::cout << DiAlla << std::endl;
    fileTag = "DiAlla";
  }
  else if(!strcmp(inName, DiAllc)){
    std::cout << DiAllc << std::endl;
    fileTag = "DiAllc";
  }
  else if(!strcmp(inName, DataA)){
    std::cout << DataA << std::endl;
    fileTag = "DataA";
  }
  else if(!strcmp(inName, DataB)){
    std::cout << DataB << std::endl;
    fileTag = "DataB";
  }
  else if(!strcmp(inName, DataC)){
    std::cout << DataC << std::endl;
    fileTag = "DataC";
  }
  else if(!strcmp(inName, DataD)){
    std::cout << DataD << std::endl;
    fileTag = "DataD";
  }
  else if(!strcmp(inName, DataE)){
    std::cout << DataE << std::endl;
    fileTag = "DataE";
  }
  else if(!strcmp(inName, DataF)){
    std::cout << DataF << std::endl;
    fileTag = "DataF";
  }
  else if(!strcmp(inName, DataF_JtCutDown)){
    std::cout << DataF_JtCutDown << std::endl;
    fileTag = "DataF_JtCutDown";
  }
  else if(!strcmp(inName, DataG)){
    std::cout << DataG << std::endl;
    fileTag = "DataG";
  }
  else if(!strcmp(inName, DataG_5)){
    std::cout << DataG_5 << std::endl;
    fileTag = "DataG_5";
  }
  else if(!strcmp(inName, DataG_JtCutDown)){
    std::cout << DataG_JtCutDown << std::endl;
    fileTag = "DataG_JtCutDown";
  }
  else if(!strcmp(inName, TESTMixing)){
    std::cout << TESTMixing << std::endl;
    fileTag = "TESTMixing";
  }
  else if(!strcmp(inName, TESTPbPb)){
    std::cout << TESTPbPb << std::endl;
    fileTag = "TESTPbPb";
  }
  else if(!strcmp(inName, TESTPP)){
    std::cout << TESTPP << std::endl;
    fileTag = "TESTPP";
  }
  else if(!strcmp(inName, TESTPbPbMC)){
    std::cout << TESTPbPbMC << std::endl;
    fileTag = "TESTPbPbMC";
  }
  else if(!strcmp(inName, TESTPPMC)){
    std::cout << TESTPPMC << std::endl;
    fileTag = "TESTPPMC";
  }
  else if(!strcmp(inName, TESTPPMCGEN)){
    std::cout << TESTPPMCGEN << std::endl;
    fileTag = "TESTPPMCGEN";
  }

  std::cout << "fileTag is: " << fileTag << std::endl;

  return;
}


Bool_t sameSign(Double_t num1, Double_t num2){
  if((num1 > 0 && num2 > 0) || (num1 < 0 && num2 < 0)) return true;

  return false;
}

void handsomeTH1(TH1 *a = 0, Int_t col = 1, Float_t size = 1, Int_t markerstyle = 20)
{
  a->SetMarkerColor(col);
  a->SetMarkerSize(size);
  a->SetMarkerStyle(markerstyle);
  a->SetLineColor(col);
  a->GetYaxis()->SetTitleOffset(1.25);
  a->GetXaxis()->SetTitleOffset(.75);
  a->GetXaxis()->CenterTitle();
  a->GetYaxis()->CenterTitle();
}


void niceTH1(TH1F* uglyTH1, float max , float min, float ndivX, float ndivY, Int_t col = 1, Float_t size = 1, Int_t style = 20)
{
  handsomeTH1(uglyTH1, col, size, style);
  uglyTH1->SetMaximum(max);
  uglyTH1->SetMinimum(min);
  uglyTH1->SetNdivisions(ndivX);
  uglyTH1->SetNdivisions(ndivY, "Y");
}


void niceTH1N(TH1F* uglyTH1, float max, float min, float ndivX, float ndivY, Int_t col = 1, Float_t size = 1, Int_t style = 20)
{
  uglyTH1->Scale(1./uglyTH1->Integral());
  niceTH1(uglyTH1, max, min, ndivX, ndivY, col, size, style);
}


Bool_t checkSetRange(Int_t setNum)
{
  if(setNum > 4 || setNum < 0){
    std::cout << "checkSetRange: setNum must be between 0-2, empty cut returned" << std::endl;
    return false;
  }

  return true;
}


Bool_t checkCentRange(Int_t centLow, Int_t centHi){
  if(centLow >= 0 && centHi >= centLow && centHi <= 199) return true;
  else{
    std::cout << "checkCentRange: centLow/centHi incorrectly specified, empty cut returned" << std::endl;
    return false;
  }
}


#endif
