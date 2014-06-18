//=============================================                                                
// Author: Chris McGinn                                                                         
//                                                                                            
// DiJet Plotter                                                              
//                                                                                            
//=============================================                                               

#include "TDatime.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLine.h"
#include "TLatex.h"
#include <iostream>

TFile* histFile_p = 0;

TFile* plotFile_p = 0;

const char* algType[3] = {"PuCalo", "VsCalo", "T"};

enum sampleType{
  kHIDATA, //0                                                                                                                 
  kHIMC,   //1                                                                                                                  
  kPPDATA, //2                                                                                                                 
  kPPMC,   //3                                                                                                                 
  kPADATA, //4                                                                                                                 
  kPAMC    //5                                                                                                             
};

Bool_t sameSign(Double_t num1, Double_t num2)
{
  if((num1 > 0 && num2 > 0) || (num1 < 0 && num2 < 0)) return true;

  return false;
}


Double_t quadSum(Double_t one, Double_t two)
{
  return TMath::Sqrt(one*one + two*two);
}


void drawPatch(float x1, float y1, float x2, float y2){
  TLegend *t1=new TLegend(x1,y1,x2,y2);
  t1->SetFillColor(kWhite);
  t1->SetBorderSize(0);
  t1->SetFillStyle(1001);
  t1->Draw("");

  return;
}


void handsomeTH1( TH1 *a=0, Int_t col =1, Float_t size=1, Int_t markerstyle=20)
{
  a->SetMarkerColor(col);
  a->SetMarkerSize(size);
  a->SetMarkerStyle(markerstyle);
  a->SetLineColor(col);
  a->GetYaxis()->SetTitleOffset(1.25);
  a->GetXaxis()->CenterTitle();
  a->GetYaxis()->CenterTitle();

  return;
}

void handsomeTH1N( TH1 *a=0, Int_t col =1)
{
  handsomeTH1(a,col);
  a->Scale(1./a->GetEntries());

  return;
}


void niceTH1(TH1F* uglyTH1, float max , float min, float ndivX, float ndivY)
{
  handsomeTH1(uglyTH1);
  uglyTH1->SetMaximum(max);
  uglyTH1->SetMinimum(min);
  uglyTH1->SetNdivisions(ndivX);
  uglyTH1->SetNdivisions(ndivY, "Y");

  return;
}


void claverCanvasSaving(TCanvas* c, TString s,TString format="gif")
{
  TDatime* date = new TDatime();
  c->SaveAs(Form("%s_%d.%s",s.Data(),date->GetDate(), format.Data()));
}


void addJtHistToPanel(TFile* file_p, const char* fileTag, TCanvas* canv_p, Int_t setNum, const char* jtVarIn, TLegend* leg, Int_t pos = 1, const char* cent = "50100", const char* opt = "E1", sampleType sType = kHIDATA)
{
  file_p->cd();
  TH1F* hist_p = (TH1F*)file_p->Get(Form("%s%s_%s_%s_h", algType[setNum], jtVarIn, cent, fileTag));

  canv_p->cd(pos);
  hist_p->DrawCopy(opt);

  if(pos == 1){
    if(sType == kPPDATA)
      leg->AddEntry(hist_p, "pp", "p");
    else if(sType == kHIDATA)
      leg->AddEntry(hist_p, "PbPb", "p");
    else
      leg->AddEntry(hist_p, "P+H", "f");
  }

  return;
}


void makeJtVarPanel(const char* filePbPbName, const char* fileTagPbPb, const char* outName, Int_t setNum, const char* jtVarIn, const char* filePPName = "", const char* fileTagPP = "", Bool_t isHighPtTrk = false, const char* filePbPbMCName = "", const char* fileTagPbPbMC = "")
{
  TFile* histPbPbFile_p = new TFile(filePbPbName, "READ");
  TFile* histPPFile_p;
  TFile* histPbPbMCFile_p;
  if(strcmp(fileTagPP, "") != 0)
    histPPFile_p = new TFile(filePPName, "READ");

  if(strcmp(fileTagPbPbMC, "") != 0)
    histPbPbMCFile_p = new TFile(filePbPbMCName, "READ");

  const char* cent[4] = {"50100", "3050", "1030", "010"};
  const char* cent2[4] = {"50-100%", "30-50%", "10-30%", "0-10%"};
  const char* cuts[4] = {"p_{T,1}>120 GeV/c", "p_{T,2}>50 GeV/c", "#Delta#phi_{1,2}>5#pi/6", "|#eta_{1}|,|#eta_{2}|<1.6"};

  TCanvas* jtVarPanel_p = new TCanvas(Form("%s%sPanel_%s", algType[setNum], jtVarIn, fileTagPbPb), Form("%s%sPanel_%s", algType[setNum], jtVarIn, fileTagPbPb), 300*4, 1*350);
  jtVarPanel_p->Divide(4, 1, 0, 0);
  std::cout << "FourPanel Init" << std::endl;

  TLegend* leg_p = new TLegend(0.60, 0.40, 0.95, 0.60);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSizePixels(28);
  leg_p->SetBorderSize(0);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  for(Int_t centIter = 0; centIter < 4; centIter++){
    if(strcmp(fileTagPbPbMC, "") != 0)
      addJtHistToPanel(histPbPbMCFile_p, fileTagPbPbMC, jtVarPanel_p, setNum, jtVarIn, leg_p, centIter+1, cent[centIter], "E1 HIST", kHIMC);


    addJtHistToPanel(histPbPbFile_p, fileTagPbPb, jtVarPanel_p, setNum, jtVarIn, leg_p, centIter+1, cent[centIter], "E1 SAME", kHIDATA);
    
    if(strcmp(fileTagPP, "") != 0)
      addJtHistToPanel(histPPFile_p, fileTagPP, jtVarPanel_p, setNum, jtVarIn, leg_p, centIter+1, "PP", "E1 SAME", kPPDATA);

    //edit here

    label_p->DrawLatex(.60, .30, cent2[centIter]);
    label_p->DrawLatex(.40, .80, cuts[centIter]);
  }

  jtVarPanel_p->cd(1);
  leg_p->Draw("SAME");

  jtVarPanel_p->cd(4);
  if(isHighPtTrk)
    label_p->DrawLatex(.40, .6, "p_{T}^{trk}>8 GeV/c");


  TFile* outFile_p = new TFile(outName, "UPDATE");
  jtVarPanel_p->Write();
  claverCanvasSaving(jtVarPanel_p, Form("pdfDir/%s%sPanel_%s", algType[setNum], jtVarIn, fileTagPbPb), "pdf");
  outFile_p->Close();

  delete outFile_p;
  delete jtVarPanel_p;

  if(strcmp(fileTagPP, "") != 0)
    delete histPPFile_p;

  delete histPbPbFile_p;

  return;
}


void makeDiJetPlots_Jet(const char* filePbPbName, const char* fileTagPbPb, const char* outName, Bool_t montecarlo = false, const char* filePPName = "", const char* fileTagPP = "", Bool_t isPercent = false, Bool_t isHighPtTrk = false, const char* filePbPbMCName = "", const char* fileTagPbPbMC = "")
{
  TH1::SetDefaultSumw2();

  Int_t jetAlgMax = 2;
  
  if(montecarlo)
    jetAlgMax = 3;
  

  for(Int_t algIter = 1; algIter < jetAlgMax; algIter++){
    makeJtVarPanel(filePbPbName, fileTagPbPb, outName, algIter, "Aj", filePPName, fileTagPP, isHighPtTrk, filePbPbMCName, fileTagPbPbMC);
  }

  return;
}

