#ifndef diJetFileTag_h
#define diJetFileTag_h

#include "TMath.h"
#include "iostream"
#include "TH1F.h"
#include "TCut.h"

#include <string>

std::string fileTag;

const std::string Di30a = "HydjetDrum_Pyquen_Dijet30_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMSKIM_20140423_0.root";

const std::string Di50a = "HydjetDrum_Pyquen_Dijet50_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMSKIM_20140423_0.root";

const std::string Di80a = "Pythia80_HydjetDrum_mix01_HiForest2_v20_CFMSKIM.root";
const std::string Di80b = "Dijet80_HydjetDrum_v27_mergedV1_CFMSKIM.root";
const std::string Di80c = "HydjetDrum_Pyquen_Dijet80_Embedded_d20140122_Track7_v2_CFMSKIM.root";
const std::string Di80d = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet6_STARTHI53_LV1_merged_forest_0_CFMSKIM.root";
const std::string Di80e = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet14_STARTHI53_LV1_merged_forest_0_CFMSKIM.root";
const std::string Di80f = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_50k_CFMSKIM.root";
const std::string Di80g = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_300k_CFMSKIM.root";
const std::string Di80h = "HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0_300k_CFMSKIM.root";
const std::string Di80i = "HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt_0_CFMSKIM.root";
const std::string Di80j = "HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt_0_CFMSKIM_20140414_0.root";

const std::string Di100a = "Dijet100_HydjetDrum_v27_mergedV1_CFMSKIM.root";
const std::string Di100b = "HydjetDrum_Pyquen_Dijet100_FOREST_Track8_Jet24_FixedPtHat_v0_0_CFMSKIM_20140323_4_0.root";
const std::string Di100c = "HydjetDrum_Pyquen_Dijet100_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMSKIM_20140423_0.root";

const std::string Di120a = "HydjetDrum_Pyquen_Dijet120_Embedded_RECO_STARTHI53_LV1_Track8_Jet21_300k_v0_merged_0_CFMSKIM.root";
const std::string Di120b = "HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHat_v0_0_CFMSKIM.root";
const std::string Di120c = "HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMSKIM_20140423_0.root";

const std::string DiAlla = "HydjetDrum_Pyquen_DijetAll_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMSKIM_20140425_0.root";

const std::string DiAllc = "HydjetDrum_Pyquen_DijetAll_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMSKIM_20140508_0.root";

const std::string EmDi80a = "PbPb_pythiaHYDJET_forest_EmEnrichedDijet80_CFMSKIM.root";

const std::string DataA = "Track8_Jet17_GR_R_53_LV6_SUB_0_CFMSKIM.root";
const std::string DataB = "hiForest_Jet80or95_GR_R_53_LV6_02Mar2014_1300CET_Track8_Jet15_0_1200k_CFMSKIM.root";
const std::string DataC = "hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21_0_700k_CFMSKIM.root";
const std::string DataD = "hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21_0_1200k_CFMSKIM.root";
const std::string DataE = "hiForest_Jet80or95_GR_R_53_LV6_03Mar2014_1600CET_CMSSW_5_3_16_merged_0_CFMSKIM.root";
const std::string DataF = "HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4-merged_0_CFMSKIM_20140430_0.root";
const std::string DataG = "HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4ANDv9-merged_0_CFMSKIM_20140505_0.root";

const std::string DataG_5 = "HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4ANDv9-merged_0_CFMSKIM_20140511_RAD5_0.root";

const std::string DataF_JtCutDown = "HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4-merged_0_CFMSKIM_20140430_JtCutDown_0.root";
const std::string DataG_JtCutDown = "HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4ANDv9-merged_0_CFMSKIM_20140508_jtCutDown_0.root";


const std::string TESTMixing = "testJtANA_30CHECK_3rdJet_0.root";
const std::string TESTPbPb = "HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4ANDv9-merged_0_DIJETANASKIM_20140619_CHECK_0.root";
const std::string TESTPP = "HiForest_pp_Jet80_v8_PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_merged_forest_0_CFMANASKIM_20140612_0.root";

const std::string TESTPbPbMC = "HydjetDrum_Pyquen_DijetMerge_FOREST_Track8_Jet24_FixedPtHatJES_v0_0_CFMANASKIM_20140618_0.root";
const std::string TESTPPMC = "HiForest_pt80_PYTHIA_ppReco_JECv85_merged_forest_0_DIJETANASKIM_20140618_0.root";
const std::string TESTPPMCGEN = "pt80_pp2013_P01_prod22_v81_merged_forest_0_DIJETANASKIM_20140703_0.root";


const std::string PYTH_HITrk = "allPtHatPythia_HITracking_DijetAnaSkim_20140903_AllHatMerge.root";
// /mnt/hadoop/cms/store/user/cfmcginn/PPMC/DijetAnaSkim/20140903/allPtHatPythia_HITracking_DijetAnaSkim_20140903_AllHatMerge.root

const std::string PYTH_PPTrk = "allPtHatPythia_PPTracking_DijetAnaSkim_20140918_AllHatMerge.root";
// /mnt/hadoop/cms/store/user/cfmcginn/PPMC/DijetAnaSkim/20140918/allPtHatPythia_PPTracking_DijetAnaSkim_20140918_AllHatMerge.root

const std::string PYTH_HYD = "allPtHatPythia_Hydjet_DijetAnaSkim_20140909_AllHatMerge.root";
// /mnt/hadoop/cms/store/user/cfmcginn/HIMC/DijetAnaSkim/20140905/allPtHatPythia_Hydjet_DijetAnaSkim_20140905_AllHatMerge.root 


const std::string DATA_PbPb_OLD = "HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4ANDv9-merged_0_DijetAnaSkim_20140903_0.root";
// /mnt/hadoop/cms/store/user/cfmcginn/HIData/DijetAnaSkim/20140903/HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4ANDv9-merged_0_DijetAnaSkim_20140903_0.root

const std::string DATA_PbPb_NEW = "PbPbForest_MatchEqR_Calo_HIHighPt_HIRun2011-14Mar2014-v4_DijetAnaSkim_20140911_0.root";
// /mnt/hadoop/cms/store/user/cfmcginn/HIData/DijetAnaSkim/20140828/PbPbForest_MatchEqR_Calo_HIHighPt_HIRun2011-14Mar2014-v4_DijetAnaSkim_20140828_0.root

const std::string DATA_PP_OLD = "PPHighPtData_ForestTag_PYTHIA_localdb_ppJEC_merged_forest_0_DijetAnaSkim_20140902_0.root";
// /mnt/hadoop/cms/store/user/cfmcginn/PPData/DijetAnaSkim/20140902/PPHighPtData_ForestTag_PYTHIA_localdb_ppJEC_merged_forest_0_DijetAnaSkim_20140902_0.root

const std::string DATA_PP_NEW = "PPHighPtData_ForestTag_PYTHIA_localdb_ppJEC_merged_forest_0_DijetAnaSkim_20140917_0.root";
// /mnt/hadoop/cms/store/user/cfmcginn/PPData/DijetAnaSkim/20140917/PPHighPtData_ForestTag_PYTHIA_localdb_ppJEC_merged_forest_0_DijetAnaSkim_20140917_0.root


const std::string name1[10] = {"0(10000, -10000, 10000)", "1(10000, -10000, 10000)", "2(10000, -10000, 10000)", "3(10000, -10000, 10000)", "4(10000, -10000, 10000)", "5(10000, -10000, 10000)", "6(10000, -10000, 10000)", "7(10000, -1000, 10000)", "8(10000, -1000, 10000)", "9(10000, -1000, 10000)"};
const std::string name2[10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};

void setFileTag(const std::string inName)
{
  if(!strcmp(inName.c_str(), Di30a.c_str())){
    std::cout << Di30a << std::endl;
    fileTag = "Di30a";
  }
  else if(!strcmp(inName.c_str(), Di50a.c_str())){
    std::cout << Di50a << std::endl;
    fileTag = "Di50a";
  }
  else if(!strcmp(inName.c_str(), Di80a.c_str())){
    std::cout << Di80a << std::endl;
    fileTag = "Di80a";
  }
  else if(!strcmp(inName.c_str(), Di80b.c_str())){
    std::cout << Di80b << std::endl;
    fileTag = "Di80b";
  }
  else if(!strcmp(inName.c_str(), Di100a.c_str())){
    std::cout << Di100a << std::endl;
    fileTag = "Di100a";
  }
  else if(!strcmp(inName.c_str(), Di100b.c_str())){
    std::cout << Di100b << std::endl;
    fileTag = "Di100b";
  }
  else if(!strcmp(inName.c_str(), Di100c.c_str())){
    std::cout << Di100c << std::endl;
    fileTag = "Di100c";
  }
  else if(!strcmp(inName.c_str(), Di120a.c_str())){
    std::cout << Di120a << std::endl;
    fileTag = "Di120a";
  }
  else if(!strcmp(inName.c_str(), Di120b.c_str())){
    std::cout << Di120b << std::endl;
    fileTag = "Di120b";
  }
  else if(!strcmp(inName.c_str(), Di120c.c_str())){
    std::cout << Di120c << std::endl;
    fileTag = "Di120c";
  }
  else if(!strcmp(inName.c_str(), EmDi80a.c_str())){
    std::cout << EmDi80a << std::endl;
    fileTag = "EmDi80a";
  }
  else if(!strcmp(inName.c_str(), Di80c.c_str())){
    std::cout << Di80c << std::endl;
    fileTag = "Di80c";
  }
  else if(!strcmp(inName.c_str(), Di80d.c_str())){
    std::cout << Di80d << std::endl;
    fileTag = "Di80d";
  }
  else if(!strcmp(inName.c_str(), Di80e.c_str())){
    std::cout << Di80e << std::endl;
    fileTag = "Di80e";
  }
  else if(!strcmp(inName.c_str(), Di80f.c_str())){
    std::cout << Di80f << std::endl;
    fileTag = "Di80f";
  }
  else if(!strcmp(inName.c_str(), Di80g.c_str())){
    std::cout << Di80g << std::endl;
    fileTag = "Di80g";
  }
  else if(!strcmp(inName.c_str(), Di80h.c_str())){
    std::cout << Di80h << std::endl;
    fileTag = "Di80h";
  }
  else if(!strcmp(inName.c_str(), Di80i.c_str())){
    std::cout << Di80i << std::endl;
    fileTag = "Di80i";
  }
  else if(!strcmp(inName.c_str(), Di80j.c_str())){
    std::cout << Di80j << std::endl;
    fileTag = "Di80j";
  }
  else if(!strcmp(inName.c_str(), DiAlla.c_str())){
    std::cout << DiAlla << std::endl;
    fileTag = "DiAlla";
  }
  else if(!strcmp(inName.c_str(), DiAllc.c_str())){
    std::cout << DiAllc << std::endl;
    fileTag = "DiAllc";
  }
  else if(!strcmp(inName.c_str(), DataA.c_str())){
    std::cout << DataA << std::endl;
    fileTag = "DataA";
  }
  else if(!strcmp(inName.c_str(), DataB.c_str())){
    std::cout << DataB << std::endl;
    fileTag = "DataB";
  }
  else if(!strcmp(inName.c_str(), DataC.c_str())){
    std::cout << DataC << std::endl;
    fileTag = "DataC";
  }
  else if(!strcmp(inName.c_str(), DataD.c_str())){
    std::cout << DataD << std::endl;
    fileTag = "DataD";
  }
  else if(!strcmp(inName.c_str(), DataE.c_str())){
    std::cout << DataE << std::endl;
    fileTag = "DataE";
  }
  else if(!strcmp(inName.c_str(), DataF.c_str())){
    std::cout << DataF << std::endl;
    fileTag = "DataF";
  }
  else if(!strcmp(inName.c_str(), DataF_JtCutDown.c_str())){
    std::cout << DataF_JtCutDown << std::endl;
    fileTag = "DataF_JtCutDown";
  }
  else if(!strcmp(inName.c_str(), DataG.c_str())){
    std::cout << DataG << std::endl;
    fileTag = "DataG";
  }
  else if(!strcmp(inName.c_str(), DataG_5.c_str())){
    std::cout << DataG_5 << std::endl;
    fileTag = "DataG_5";
  }
  else if(!strcmp(inName.c_str(), DataG_JtCutDown.c_str())){
    std::cout << DataG_JtCutDown << std::endl;
    fileTag = "DataG_JtCutDown";
  }
  else if(!strcmp(inName.c_str(), TESTMixing.c_str())){
    std::cout << TESTMixing << std::endl;
    fileTag = "TESTMixing";
  }
  else if(!strcmp(inName.c_str(), TESTPbPb.c_str())){
    std::cout << TESTPbPb << std::endl;
    fileTag = "TESTPbPb";
  }
  else if(!strcmp(inName.c_str(), TESTPP.c_str())){
    std::cout << TESTPP << std::endl;
    fileTag = "TESTPP";
  }
  else if(!strcmp(inName.c_str(), TESTPbPbMC.c_str())){
    std::cout << TESTPbPbMC << std::endl;
    fileTag = "TESTPbPbMC";
  }
  else if(!strcmp(inName.c_str(), TESTPPMC.c_str())){
    std::cout << TESTPPMC << std::endl;
    fileTag = "TESTPPMC";
  }
  else if(!strcmp(inName.c_str(), TESTPPMCGEN.c_str())){
    std::cout << TESTPPMCGEN << std::endl;
    fileTag = "TESTPPMCGEN";
  }
  else if(!strcmp(inName.c_str(), PYTH_HITrk.c_str())){
    std::cout << PYTH_HITrk << std::endl;
    fileTag = "PYTH_HITrk";
  }
  else if(!strcmp(inName.c_str(), PYTH_PPTrk.c_str())){
    std::cout << PYTH_PPTrk << std::endl;
    fileTag = "PYTH_PPTrk";
  }
  else if(!strcmp(inName.c_str(), PYTH_HYD.c_str())){
    std::cout << PYTH_HYD << std::endl;
    fileTag = "PYTH_HYD";
  }
  else if(!strcmp(inName.c_str(), DATA_PbPb_OLD.c_str())){
    std::cout << DATA_PbPb_OLD << std::endl;
    fileTag = "DATA_PbPb_OLD";
  }
  else if(!strcmp(inName.c_str(), DATA_PbPb_NEW.c_str())){
    std::cout << DATA_PbPb_NEW << std::endl;
    fileTag = "DATA_PbPb_NEW";
  }
  else if(!strcmp(inName.c_str(), DATA_PP_OLD.c_str())){
    std::cout << DATA_PP_OLD << std::endl;
    fileTag = "DATA_PP_OLD";
  }
  else if(!strcmp(inName.c_str(), DATA_PP_NEW.c_str())){
    std::cout << DATA_PP_NEW << std::endl;
    fileTag = "DATA_PP_NEW";
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
