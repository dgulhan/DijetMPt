#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TString.h"

#include <iostream>
#include <fstream>
#include <string>

#include "TMath.h"
#include "TLatex.h"
#include "TDatime.h"
#include "TLegend.h"

const std::string algType[7] = {"Pu3Calo", "Pu4Calo", "Pu5Calo", "Vs2Calo", "Vs3Calo", "Vs4Calo", "Vs5Calo"};
const std::string algTypePP[7] = {"Vs3Calo", "Vs4Calo", "Vs5Calo", "Vs2Calo", "Vs3Calo", "Vs4Calo", "Vs5Calo"};
const Int_t totJtPtCut = 50;

const std::string histType1[7] = {"LeadPt", "SubleadPt", "ThirdPt", "LeadEta", "SubleadEta", "DelPhi12", "AJ"};
const std::string histType2[7] = {"p_{T,1}", "p_{T,2}", "p_{T,3}", "#eta_{1}", "#eta_{2}", "#Delta#phi_{1,2}", "A_{J}"};
const std::string etaCut[7] = {"|#eta_{1,2}|<2.0", "|#eta_{1,2}|<2.0", "|#eta_{1,2,3}|<2.0", "|#eta_{1,2}|<2.0", "|#eta_{1,2}|<2.0", "|#eta_{1,2}|<2.0", "|#eta_{1,2}|<1.6"};
std::string phiCut[7] = {"", "", "p_{T,3}>30 GeV/c", "", "", "", "#Delta#phi>5#pi/6"};


const Float_t histYMax[7] = {.99, .99, .99, .0999, .0999, .9, .2499};
const Float_t histYMin[7] = {.0000002, .0000002, .0000002, .0001, .0001, .00002, .0001};

const Float_t histXMax[7] = {699.5, 699.5, 399.5, 1.9999, 1.9999, 3.14, .9999};
const Float_t histXMin[7] = {119.5, 49.5, 29.5, -1.9999, -1.9999, 0.0001, .0001};

const Bool_t histLog[7] = {true, true, true, false, false, true, false};

const Float_t legXLow[7] = {.65, .65, .65, .60, .60, .30, .60};
const Float_t legXHi[7] = {.85, .85, .85, .85, .85, .55, .85};
const Float_t legYLow[7] = {.65, .65, .65, .55, .55, .55, .55};
const Float_t legYHi[7] = {.85, .85, .85, .75, .75, .75, .75};

const Int_t YPanelNum[7] = {1, 1, 1, 1, 1, 1, 2};
const Float_t XOff[7] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.5};
const Float_t YOff[7] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 2.5};

void handsomeTH1( TH1 *a=0, Int_t col =1, Float_t size=1, Int_t markerstyle=20)
{
  a->SetMarkerColor(col);
  a->SetMarkerSize(size);
  a->SetMarkerStyle(markerstyle);
  a->SetLineColor(col);
}

void niceTH1(TH1F* uglyTH1, float max , float min, float ndivX, float ndivY)
{
  handsomeTH1(uglyTH1);
  uglyTH1->SetMaximum(max);
  uglyTH1->SetMinimum(min);
  uglyTH1->GetXaxis()->SetNdivisions(ndivX);
  uglyTH1->GetYaxis()->SetNdivisions(ndivY);
  return;
}

void claverCanvasSaving(TCanvas* c, TString s,TString format="gif") {
  TDatime* date = new TDatime();
  c->SaveAs(Form("%s_%d.%s",s.Data(),date->GetDate(), format.Data()));
}


void drawMeanLine(Float_t histMean, Int_t lineColor, Float_t lineMin, Float_t lineMax, Int_t lineStyle = 2)
{
  TLine* oneLine_p = new TLine(histMean, lineMin, histMean, lineMax);
  oneLine_p->SetLineColor(lineColor);
  oneLine_p->SetLineStyle(lineStyle);

  oneLine_p->Draw("Same");
}


void SetTitleLabel(TH1F* inHist_p, Int_t titleSize, Int_t labelSize, std::string xTitle, Float_t xOff, std::string yTitle, Float_t yOff)
{
  inHist_p->GetXaxis()->SetTitleOffset(xOff);
  inHist_p->GetYaxis()->SetTitleOffset(yOff);
  inHist_p->GetXaxis()->CenterTitle();
  inHist_p->GetYaxis()->CenterTitle();

  inHist_p->GetXaxis()->SetTitleFont(43);
  inHist_p->GetYaxis()->SetTitleFont(43);
  inHist_p->GetXaxis()->SetLabelFont(43);
  inHist_p->GetYaxis()->SetLabelFont(43);

  inHist_p->GetXaxis()->SetTitleColor(1);

  inHist_p->GetXaxis()->SetTitleSize(titleSize);
  inHist_p->GetYaxis()->SetTitleSize(titleSize);
  inHist_p->GetXaxis()->SetLabelSize(labelSize);
  inHist_p->GetYaxis()->SetLabelSize(labelSize);

  inHist_p->SetXTitle(xTitle.c_str());
  inHist_p->SetYTitle(yTitle.c_str());
}

void plotDijetAJ(const std::string histFileName, const std::string histFileTag, const std::string dataFileName, const std::string dataFileTag, const std::string ppFileName, const std::string ppFileTag, const std::string ppMCFileName, const std::string ppMCFileTag, Int_t setNum, const Int_t histNum)
{
  const Int_t titleSize = 26;
  const Int_t labelSize = 22;

  const Float_t labelPosX[5] = {.50, .35, .35, .35, .35};
  
  const std::string centString[5] = {"50100", "3050", "2030", "1020", "010"};
  const std::string centString2[5] = {"50-100%", "30-50%", "20-30%", "10-20%", "0-10%"};
  const std::string cutString[5] = {"p_{T,1}>120 GeV/c", "p_{T,2}>50 GeV/c", etaCut[histNum], phiCut[histNum], algType[setNum].c_str()};

  TFile* ppMCFile_p = new TFile(ppMCFileName.c_str(), "UPDATE");
  TH1F* getPPMC_p;

  getPPMC_p = (TH1F*)ppMCFile_p->Get(Form("%s%s_PP_%s_h", algTypePP[setNum].c_str(), histType1[histNum].c_str(), ppMCFileTag.c_str()));
  niceTH1(getPPMC_p, histYMax[histNum], histYMin[histNum], 505, 505);

  TFile* ppFile_p = new TFile(ppFileName.c_str(), "UPDATE");
  TH1F* getPP_p;

  getPP_p = (TH1F*)ppFile_p->Get(Form("%s%s_PP_%s_h", algTypePP[setNum].c_str(), histType1[histNum].c_str(), ppFileTag.c_str()));
  niceTH1(getPP_p, histYMax[histNum], histYMin[histNum], 505, 505);


  TFile* dataFile_p = new TFile(dataFileName.c_str(), "UPDATE");
  TH1F* getData_p[5];

  for(Int_t iter = 0; iter < 5; iter++){
    getData_p[iter] = (TH1F*)dataFile_p->Get(Form("%s%s_%s_%s_h", algType[setNum].c_str(), histType1[histNum].c_str(), centString[iter].c_str(), dataFileTag.c_str()));
    niceTH1(getData_p[iter], histYMax[histNum], histYMin[histNum], 505, 505);
  }

  TFile* histFile_p = new TFile(histFileName.c_str(), "UPDATE");
  TH1F* getHist_p[5];

  for(Int_t iter = 0; iter < 5; iter++){
    getHist_p[iter] = (TH1F*)histFile_p->Get(Form("%s%s_%s_%s_h", algType[setNum].c_str(), histType1[histNum].c_str(), centString[iter].c_str(), histFileTag.c_str()));
    niceTH1(getHist_p[iter], histYMax[histNum], histYMin[histNum], 505, 505);
  }

  TCanvas* plotCanv_p = new TCanvas(Form("%s%s_%s_c", algType[setNum].c_str(), histType1[histNum].c_str(), histFileTag.c_str()), Form("%s%s_%s_c", algType[setNum].c_str(), histType1[histNum].c_str(), histFileTag.c_str()), 5*300, YPanelNum[histNum]*350);
  plotCanv_p->Divide(5, YPanelNum[histNum], 0, 0);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(23);

  for(Int_t iter = 0; iter < 5; iter++){
    plotCanv_p->cd(iter+1);
    if(histLog[histNum]) gPad->SetLogy();
    SetTitleLabel(getHist_p[iter], titleSize, labelSize, histType2[histNum], XOff[histNum], "Event Fraction", YOff[histNum]);
    getHist_p[iter]->SetFillColor(11);
    getHist_p[iter]->GetXaxis()->SetLimits(histXMin[histNum], histXMax[histNum]);
    getHist_p[iter]->DrawCopy("HIST E1");
    getHist_p[iter]->DrawCopy("E1 SAME");
    if(histNum == 6) drawMeanLine(getHist_p[iter]->GetMean(), 1, histYMin[histNum], histYMax[histNum]);
    getData_p[iter]->SetMarkerColor(kRed);
    getData_p[iter]->GetXaxis()->SetLimits(histXMin[histNum], histXMax[histNum]);
    getData_p[iter]->DrawCopy("E1 SAME");
    if(histNum == 6) drawMeanLine(getData_p[iter]->GetMean(), kRed, histYMin[histNum], histYMax[histNum]);
    getPP_p->SetMarkerColor(kBlue);
    getPP_p->GetXaxis()->SetLimits(histXMin[histNum], histXMax[histNum]);
    getPP_p->DrawCopy("E1 SAME");
    if(histNum == 6) drawMeanLine(getPP_p->GetMean(), kBlue, histYMin[histNum], histYMax[histNum]);
    getPPMC_p->SetMarkerStyle(24);
    getPPMC_p->GetXaxis()->SetLimits(histXMin[histNum], histXMax[histNum]);
    getPPMC_p->DrawCopy("E1 SAME");

    label_p->DrawLatex(labelPosX[iter], .92, centString2[iter].c_str());
    label_p->DrawLatex(labelPosX[iter], .84, cutString[iter].c_str());
  }

  TLegend* leg_p = new TLegend(legXLow[histNum], legYLow[histNum], legXHi[histNum], legYHi[histNum]);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSizePixels(16);
  leg_p->SetBorderSize(0);

  leg_p->AddEntry(getHist_p[0], "PYT+HYD", "F P");
  leg_p->AddEntry(getData_p[0], "PbPb", "P");
  leg_p->AddEntry(getPP_p, "PP", "P");
  leg_p->AddEntry(getPPMC_p, "PYTHIA", "P");

  plotCanv_p->cd(1);
  leg_p->Draw("SAME");


  if(histNum == 6){
    TH1::SetDefaultSumw2();
    for(Int_t iter = 0; iter < 5; iter++){
      plotCanv_p->cd(10-iter);

      niceTH1(getData_p[4 - iter], 3.9999, 0.0001, 505, 505);
      SetTitleLabel(getData_p[4 - iter], titleSize, labelSize, histType2[histNum], XOff[histNum], "R_{cp}", YOff[histNum]);

      getData_p[4-iter]->SetMarkerColor(kRed);
      getData_p[4 - iter]->Divide(getData_p[0]);
      getData_p[4 - iter]->DrawCopy("E1");

      niceTH1(getHist_p[4 - iter], 3.9999, 0.0001, 505, 505);
      getHist_p[4 - iter]->Divide(getHist_p[0]);
      getHist_p[4 - iter]->DrawCopy("E1 SAME");
    }

  }

  plotCanv_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(plotCanv_p, Form("pdfDir/%s%s_%s", algType[setNum].c_str(), histType1[histNum].c_str(), histFileTag.c_str()), "pdf");


  delete leg_p;
  delete label_p;
  delete plotCanv_p;

  histFile_p->Close();
  delete histFile_p;

  dataFile_p->Close();
  delete dataFile_p;

  return;
}


void makeDijetPlots_Jet(const std::string histFileName, const std::string histFileTag, const std::string dataFileName, const std::string dataFileTag, const std::string ppFileName, const std::string ppFileTag, const std::string ppMCFileName, const std::string ppMCFileTag)
{
  TH1::SetDefaultSumw2();

  for(Int_t iter = 0; iter < 7; iter++){
    for(Int_t histIter = 0; histIter < 7; histIter++){
      plotDijetAJ(histFileName, histFileTag, dataFileName, dataFileTag, ppFileName, ppFileTag, ppMCFileName, ppMCFileTag, iter, histIter);
    }
  }

}
