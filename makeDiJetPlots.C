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


void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns, const Int_t rows, const Float_t leftOffset, const Float_t bottomOffset, const Float_t leftMargin, const Float_t bottomMargin, const Float_t edge, const char* CNCR = ""){
  if(canv==0){
    Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];

  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];
  Float_t PadWidth = (1.0-leftOffset)/((1.0/(1.0-leftMargin)) + (1.0/(1.0-edge))+(Float_t)columns-2.0);
  Float_t PadHeight = (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) + (1.0/(1.0-edge))+(Float_t)rows-2.0);
  Xlow[0] = leftOffset;
  Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
  Xup[columns-1] = 1;
  Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

  Yup[0] = 1;
  Ylow[0] = 1.0-PadHeight/(1.0-edge);
  Ylow[rows-1] = bottomOffset;
  Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

  for(Int_t i=1; i<columns-1; i++){
    Xlow[i] = Xup[0] + (i-1)*PadWidth;
    Xup[i] = Xup[0] + (i)*PadWidth;
  }
  Int_t ct = 0;
  for(Int_t i=rows-2; i>0; i--){
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }

  TString padName;
  for(Int_t i=0; i<columns; i++){
    for(Int_t j=0; j<rows; j++){
      canv->cd();
      padName = Form("p_%d_%d",i,j);

      if(i==0 && j==1){
        pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j],Xup[i]*0.890,Yup[j]*0.935);
      }
      else if(i==0 && j==0){
        if(strcmp(CNCR, "") ==0)
          pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j]*0.933,Xup[i],Yup[j]);
        else
          pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j]*0.934,Xup[i],Yup[j]*.9985);
      }
      else if(i==1 && j==1){
        if(strcmp(CNCR, "") == 0)
          pad[i][j] = new TPad(padName.Data(),padName.Data(),0.815*Xlow[i],Ylow[j],Xup[i],Yup[j]*.995);
        else
          pad[i][j] = new TPad(padName.Data(),padName.Data(),0.80*Xlow[i],Ylow[j],Xup[i],Yup[j]*.995);
      }
      else if(j == 0){
        pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j],Xup[i],Yup[j]*.998);
      }
      else pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j],Xup[i],Yup[j]*.995);

      if(i==0){
        if(j == 0)
          pad[i][j]->SetLeftMargin(leftMargin*1.2);
       else
          pad[i][j]->SetLeftMargin(leftMargin);
      }
      else if(i==1 && j==1){
        if(strcmp(CNCR, "") == 0)
          pad[i][j]->SetLeftMargin(PadWidth);
        else
          pad[i][j]->SetLeftMargin(PadWidth*.65);
      }
      else pad[i][j]->SetLeftMargin(0);

      if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
      else pad[i][j]->SetRightMargin(0);

      if(j==0){
        if(i==0)pad[i][j]->SetTopMargin(edge);
        else pad[i][j]->SetTopMargin(edge);
      }
      else pad[i][j]->SetTopMargin(0);

      if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
      else if(i==0 && j==0) pad[i][j]->SetBottomMargin(0.17*PadHeight);
      else pad[i][j]->SetBottomMargin(0);


      pad[i][j]->Draw();
      pad[i][j]->cd();
      pad[i][j]->SetNumber(columns*j+i+1);

    }
  }
  pad[0][0]->cd();

  return;
}


Double_t sumYForPTStack(Double_t dIn = 0, Double_t comp1 = 0, Double_t comp2 = 0, Double_t comp3 = 0, Double_t comp4 = 0)
{
  Double_t dOut = dIn;

  if(sameSign(comp1, dOut))
    dOut += comp1;

  if(sameSign(comp2, dOut))
    dOut += comp2;

  if(sameSign(comp3, dOut))
    dOut += comp3;

  if(sameSign(comp4, dOut))
    dOut += comp4;

  return dOut;
}


void makeHistForPtStack(TH1F* h_p[6], Int_t pos = 4, const char* Tight = "", const char* CNCR = "", Bool_t isPercent = false)
{
  Int_t nBins = 4;

  if(!strcmp(CNCR, "R") || !strcmp(CNCR, "RD") || !strcmp(CNCR, "RU"))
    nBins = 10;
  else if(strcmp(Tight, "") != 0)
    nBins = 8;

  for(Int_t iter = 0; iter < nBins; iter++){
    h_p[0]->SetBinContent(iter + 1, sumYForPTStack(h_p[0]->GetBinContent(iter+1), h_p[1]->GetBinContent(iter+1), h_p[2]->GetBinContent(iter+1), h_p[3]->GetBinContent(iter+1), h_p[4]->GetBinContent(iter+1)));

    h_p[1]->SetBinContent(iter + 1, sumYForPTStack(h_p[1]->GetBinContent(iter+1), h_p[2]->GetBinContent(iter+1), h_p[3]->GetBinContent(iter+1), h_p[4]->GetBinContent(iter+1)));

    h_p[2]->SetBinContent(iter + 1, sumYForPTStack(h_p[2]->GetBinContent(iter+1), h_p[3]->GetBinContent(iter+1), h_p[4]->GetBinContent(iter+1)));

    h_p[3]->SetBinContent(iter + 1, sumYForPTStack(h_p[3]->GetBinContent(iter+1), h_p[4]->GetBinContent(iter+1)));
  }

  const char* xTitle;

  if(!strcmp(CNCR, "R") || !strcmp(CNCR, "RD") || !strcmp(CNCR, "RU"))
    xTitle = "#DeltaR";
  else{
    if(!isPercent)
      xTitle = "A_{J}";
    else
      xTitle = "A_{J}%";
  }

  h_p[4]->SetXTitle(xTitle);
  h_p[3]->SetXTitle(xTitle);
  h_p[2]->SetXTitle(xTitle);
  h_p[1]->SetXTitle(xTitle);
  h_p[0]->SetXTitle(xTitle);
  h_p[5]->SetXTitle(xTitle);

  if(pos == 1){
    h_p[4]->SetYTitle("   <#slash{p}_{T}^{||}> (GeV/c)");
    h_p[3]->SetYTitle("   <#slash{p}_{T}^{||}> (GeV/c)");
    h_p[2]->SetYTitle("   <#slash{p}_{T}^{||}> (GeV/c)");
    h_p[1]->SetYTitle("   <#slash{p}_{T}^{||}> (GeV/c)");
    h_p[0]->SetYTitle("   <#slash{p}_{T}^{||}> (GeV/c)");
    h_p[5]->SetYTitle("   <#slash{p}_{T}^{||}> (GeV/c)");
  }

  return;
}

void drawHistToPTStack(TH1F* drawHist_p, Int_t color, const char* drawOpt, Bool_t isSub = false, const char* CNCR = "", Bool_t isPercent = false)
{
  if(isSub){
    if(!strcmp(CNCR, "R") || !strcmp(CNCR, "RD") || !strcmp(CNCR, "RU"))
      niceTH1(drawHist_p, 4.999, -10., 505, 503);
    else{
      if(!isPercent)
	niceTH1(drawHist_p, 59.999, -60., 505, 406);
      else
	niceTH1(drawHist_p, 59.999, -60., 504, 406);
    }
  }

  drawHist_p->SetFillColor(color);
  drawHist_p->SetMarkerStyle(6);
  drawHist_p->SetMarkerSize(.5);
  drawHist_p->DrawCopy(drawOpt);
  drawHist_p->DrawCopy("E1 SAME");

  return;
}

void drawFullStack(TH1F* h_p[6], Int_t color, Int_t style, TLegend* leg_p = 0, Bool_t isSub = false, const char* CNCR = "", Bool_t isPercent = false)
{
  drawHistToPTStack(h_p[0], kBlue - 9, "E1 HIST", isSub, CNCR, isPercent);
  drawHistToPTStack(h_p[1], kYellow - 9, "E1 HIST SAME", isSub, CNCR, isPercent);
  drawHistToPTStack(h_p[2], kOrange + 1, "E1 HIST SAME", isSub, CNCR, isPercent);
  drawHistToPTStack(h_p[3], kGreen + 3, "E1 HIST SAME", isSub, CNCR, isPercent);
  drawHistToPTStack(h_p[4], kRed + 1, "E1 HIST SAME", isSub, CNCR, isPercent);

  if(isSub){
    if(!strcmp(CNCR, "R") || !strcmp(CNCR, "RD") || !strcmp(CNCR, "RU"))
      niceTH1(h_p[5], 4.999, -10., 505, 403);
    else{
      if(!isPercent)
	niceTH1(h_p[5], 59.999, -60., 505, 406);
      else 
	niceTH1(h_p[5], 59.999, -60., 504, 406);
    }
  }

  h_p[5]->SetFillColor(color);
  h_p[5]->SetMarkerStyle(style);

  h_p[5]->DrawCopy("SAME E1");

  if(leg_p != 0){
    leg_p->AddEntry(h_p[0], "0.5 - 1.0", "F");
    leg_p->AddEntry(h_p[1], "1.0 - 2.0", "F");
    leg_p->AddEntry(h_p[2], "2.0 - 4.0", "F");
    leg_p->AddEntry(h_p[3], "4.0 - 8.0", "F");
    leg_p->AddEntry(h_p[4], "8.0 - 300.0", "F");
  }

  return;
}


void makeSysError(Float_t sysArr[4], TH1F* hist_p)
{
  for(Int_t iter = 0; iter < hist_p->GetNbinsX(); iter++){
    Float_t yVal = hist_p->GetBinContent(iter+1);
    Float_t sys = sysArr[iter];
    TLine* l = new TLine(hist_p->GetBinLowEdge(iter+1) + .01, yVal - TMath::Sqrt(sys*sys), hist_p->GetBinLowEdge(iter+2) - .01, yVal - TMath::Sqrt(sys*sys));
    l->SetLineColor(1);
    l->Draw();
    l->DrawLine(hist_p->GetBinLowEdge(iter+1) + .01, yVal - TMath::Sqrt(sys*sys), hist_p->GetBinLowEdge(iter+1) + .01, yVal - TMath::Sqrt(sys*sys) + 2);
    l->DrawLine(hist_p->GetBinLowEdge(iter + 2) - .01, yVal - TMath::Sqrt(sys*sys), hist_p->GetBinLowEdge(iter + 2) - .01, yVal - TMath::Sqrt(sys*sys) + 2);
    l->DrawLine(hist_p->GetBinLowEdge(iter + 1) + .01, yVal + TMath::Sqrt(sys*sys), hist_p->GetBinLowEdge(iter + 2) -.01, yVal + TMath::Sqrt(sys*sys));

    l->DrawLine(hist_p->GetBinLowEdge(iter + 1) + .01, yVal + TMath::Sqrt(sys*sys) - 2, hist_p->GetBinLowEdge(iter + 1) + .01, yVal + TMath::Sqrt(sys*sys));

    l->DrawLine(hist_p->GetBinLowEdge(iter + 2) - .01, yVal + TMath::Sqrt(sys*sys) - 2, hist_p->GetBinLowEdge(iter + 2) - .01, yVal + TMath::Sqrt(sys*sys));
  }

  return;
}


 void makeImbPtStack(const char* filePbPbName, const char* fileTagPbPb, const char* outName, const char* gorr, Int_t setNum, const char* Corr = "", const char* CNCR = "", Bool_t montecarlo = false, const char* filePPName = "", const char* fileTagPP = "", const char* Tight = "", Bool_t isPercent = false)
{
  TFile* histPbPbFile_p = new TFile(filePbPbName, "READ");
  TFile* histPPFile_p;
  if(strcmp(fileTagPP, "") != 0)
    histPPFile_p = new TFile(filePPName, "READ");

  const char* mcLabel[4] = {"PYTHIA", "PYTHIA + HYDJET", "PYTHIA + HYDJET", "(P + H) - P"};
  const char* dataLabel[4] = {"pp 5.3 pb^{-1}", "PbPb 150 #mub^{-1}", "PbPb", "PbPb - pp"};

  const char* overLabel[4];
  Float_t overCoord[4] = {.84, .76, .90, .82};

  if(strcmp(CNCR, "") != 0){
    for(Int_t coordIter = 0; coordIter < 4; coordIter++){
      overCoord[coordIter] += .04;
    }
  }

  for(Int_t iter = 0; iter < 4; iter++){
    if(montecarlo)
      overLabel[iter] = mcLabel[iter];
    else
      overLabel[iter] = dataLabel[iter];
  }

  //Grab hists for stack

  TH1F* hist1_p[6];
  TH1F* hist2_p[6];
  TH1F* hist3_p[6];
  TH1F* hist4_p[6];
  TH1F* histPP_p[6];

  const char* FPT[6] = {"0_1", "1_2", "2_4", "4_8", "8_100", "F"};

  for(Int_t histIter = 0; histIter < 6; histIter++){
    if(strcmp(fileTagPP, "") != 0)
      histPP_p[histIter] = (TH1F*)histPPFile_p->Get(Form("%s%sImbProjA%s%s%s%s_PP_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPP));

    if(!strcmp(CNCR, "")){
      hist1_p[histIter] = (TH1F*)histPbPbFile_p->Get(Form("%s%sImbProjA%s%s%s%s_50100_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPbPb));
      hist2_p[histIter] = (TH1F*)histPbPbFile_p->Get(Form("%s%sImbProjA%s%s%s%s_3050_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPbPb));
      hist3_p[histIter] = (TH1F*)histPbPbFile_p->Get(Form("%s%sImbProjA%s%s%s%s_1030_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPbPb));
      hist4_p[histIter] = (TH1F*)histPbPbFile_p->Get(Form("%s%sImbProjA%s%s%s%s_010_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPbPb));
    }
    else{
      hist1_p[histIter] = (TH1F*)histPbPbFile_p->Get(Form("%s%sImbProjA%s%s%s%s_30100_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPbPb));
      hist2_p[histIter] = (TH1F*)histPbPbFile_p->Get(Form("%s%sImbProjA%s%s%s%s_030_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPbPb));
    }
  }

  //Draw first PP panel

  makeHistForPtStack(histPP_p, 1, Tight, CNCR, isPercent);

  TCanvas* profPanel_p;

  if(!strcmp(CNCR, "")){
    profPanel_p = new TCanvas(Form("%s%sImb%s%s%sPTStackPP_%s_c", gorr, algType[setNum], CNCR, Corr, Tight, fileTagPbPb), Form("%s%sImb%s%s%sPTStackPP_%s_c", gorr, algType[setNum], CNCR, Corr, Tight, fileTagPbPb), 5*300, 700);
    makeMultiPanelCanvas(profPanel_p, 5, 2, 0.0, 0.0, 0.2, 0.2, 0.05);
    std::cout << "FivePanel Init" << std::endl;
  }
  else{
    profPanel_p = new TCanvas(Form("%s%sImb%s%s%sPTStackPP_%s_c", gorr, algType[setNum], CNCR, Corr, Tight, fileTagPbPb), Form("%s%sImb%s%s%sPTStackPP_%s_c", gorr, algType[setNum], CNCR, Corr, Tight, fileTagPbPb), 300*3, 700);
    makeMultiPanelCanvas(profPanel_p, 3, 2, 0.0, 0.0, 0.2, 0.2, 0.01, CNCR);
    std::cout << "ThreePanel Init" << std::endl;
  }

  //Make legend

  TLegend* legA_p = new TLegend(0.25, 0.18, 0.99, 0.88);
  TLegend* legB_p = new TLegend(0.25, 0.1, 0.55, 0.40);

  legA_p->SetFillColor(0);
  legA_p->SetFillStyle(0);
  legA_p->SetTextFont(43);
  legA_p->SetTextSizePixels(28);
  legA_p->SetBorderSize(0);

  legB_p->SetFillColor(0);
  legB_p->SetFillStyle(0);
  legB_p->SetTextFont(43);
  legB_p->SetTextSizePixels(28);
  legB_p->SetBorderSize(0);

  profPanel_p->cd(1);

  histPP_p[0]->GetYaxis()->SetTitleOffset(2.2);

  drawFullStack(histPP_p, 0, 25, legA_p, false, CNCR, isPercent);

  Float_t sysAPP[4] = {2.2, 3.3, 4.4, 5.5};
  if(!montecarlo && !strcmp(CNCR, "") && !strcmp(Tight, "")) makeSysError(sysAPP, histPP_p[5]);

  TLine* zeroLine_p;
  if(!isPercent)  zeroLine_p = new TLine(0., 0., 0.5, 0.);
  else zeroLine_p = new TLine(0., 0., 100.00, 0.);

  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(1);
  zeroLine_p->Draw();

  TLatex* label1_p = new TLatex();
  label1_p->SetNDC();
  label1_p->SetTextFont(43);
  label1_p->SetTextSizePixels(28);

  label1_p->DrawLatex(.28, overCoord[0], Form("%s", overLabel[0]));
  label1_p->DrawLatex(.28, overCoord[1], "CMS Preliminary");

  legB_p->Draw("SAME");

  profPanel_p->cd(2);

  makeHistForPtStack(hist1_p, 2, Tight, CNCR, isPercent);

  drawFullStack(hist1_p, 0, 28, 0, false, CNCR, isPercent);

  Float_t sysA50100[4] = {4.3, 4.7, 5.2, 5.8};
  if(!montecarlo && !strcmp("", CNCR) && !strcmp(Tight, "")) makeSysError(sysA50100, hist1_p[5]);

  zeroLine_p->Draw();

  label1_p->DrawLatex(.05, overCoord[0], Form("%s", overLabel[1]));
  if(strcmp("", CNCR) == 0)
    label1_p->DrawLatex(.05, overCoord[1], "50-100%");
  else
    label1_p->DrawLatex(.05, overCoord[1], "30-100%");

  label1_p->DrawLatex(.05, .05, "#sqrt{s_{NN}} = 2.76 TeV");

  if(!strcmp(CNCR, ""))
    profPanel_p->cd(6);
  else
    profPanel_p->cd(4);

  TH1F* histDum_p = new TH1F("histDum_p", "histDum_p", 10, 0, 1);
  histDum_p->SetMarkerStyle(28);
  legA_p->AddEntry(histDum_p, "> 0.5", "p");
  legB_p->AddEntry(histPP_p[5], "pp", "p");
  legB_p->AddEntry(histDum_p, "PbPb", "p");
  legA_p->Draw("SAME");

  label1_p->DrawLatex(.30, .92, "p_{T}^{trk} (|#eta|<2.4)");

  profPanel_p->cd(3);

  //Draw second PbPb hist

  makeHistForPtStack(hist2_p, 3, Tight, CNCR, isPercent);

  drawFullStack(hist2_p, 0, 28, 0, false, CNCR, isPercent);

  Float_t sysA3050[4] = {3.8, 4.4, 5.8, 6.5};
  if(!montecarlo && !strcmp("", CNCR) && !strcmp(Tight, "")) makeSysError(sysA3050, hist2_p[5]);

  zeroLine_p->Draw();

  label1_p->DrawLatex(.05, overCoord[0], Form("%s", overLabel[2]));
  if(!strcmp("", CNCR))
    label1_p->DrawLatex(.05, overCoord[1], "30-50%");
  else
    label1_p->DrawLatex(.05, overCoord[1], "0-30%");

  //Draw third and fourth PbPb panels, if applicable

  if(!strcmp(CNCR, "")){
    profPanel_p->cd(4);

    makeHistForPtStack(hist3_p, 4, Tight, CNCR, isPercent);

    drawFullStack(hist3_p, 0, 28, 0, false, CNCR, isPercent);

    Float_t sysA1030[4] = {2.3, 2.9, 3.7, 4.3};
    if(!montecarlo && !strcmp("", CNCR) && !strcmp(Tight, "")) makeSysError(sysA1030, hist3_p[5]);

    zeroLine_p->Draw();

    label1_p->DrawLatex(.05, overCoord[0], Form("%s", overLabel[2]));
    label1_p->DrawLatex(.05, overCoord[1], "10-30%");

    profPanel_p->cd(5);

    makeHistForPtStack(hist4_p, 5, Tight, CNCR, isPercent);

    drawFullStack(hist4_p, 0, 28, 0, false, CNCR, isPercent);

    Float_t sysA010[4] = {3.1, 3.5, 4.2, 5.5};
    if(!montecarlo && !strcmp("", CNCR) && !strcmp(Tight, "")) makeSysError(sysA010, hist4_p[5]);

    zeroLine_p->Draw();

    label1_p->DrawLatex(.05, overCoord[0], Form("%s", overLabel[2]));
    label1_p->DrawLatex(.05, overCoord[1], "0-10%");
  }

  for(Int_t histIter = 0; histIter < 6; histIter++){
    if(strcmp(fileTagPP, "") != 0)
      histPP_p[histIter] = (TH1F*)histPPFile_p->Get(Form("%s%sImbProjA%s%s%s%s_PP_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPP));

    if(!strcmp(CNCR, "")){
      hist1_p[histIter] = (TH1F*)histPbPbFile_p->Get(Form("%s%sImbProjA%s%s%s%s_50100_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPbPb));
      hist2_p[histIter] = (TH1F*)histPbPbFile_p->Get(Form("%s%sImbProjA%s%s%s%s_3050_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPbPb));
      hist3_p[histIter] = (TH1F*)histPbPbFile_p->Get(Form("%s%sImbProjA%s%s%s%s_1030_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPbPb));
      hist4_p[histIter] = (TH1F*)histPbPbFile_p->Get(Form("%s%sImbProjA%s%s%s%s_010_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPbPb));

      hist3_p[histIter]->Add(histPP_p[histIter], -1);
      hist4_p[histIter]->Add(histPP_p[histIter], -1);
    }
    else{
      hist1_p[histIter] = (TH1F*)histPbPbFile_p->Get(Form("%s%sImbProjA%s%s%s%s_30100_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPbPb));
      hist2_p[histIter] = (TH1F*)histPbPbFile_p->Get(Form("%s%sImbProjA%s%s%s%s_030_%s_h", gorr, algType[setNum], CNCR, FPT[histIter], Corr, Tight, fileTagPbPb));
    }

    hist1_p[histIter]->Add(histPP_p[histIter], -1);
    hist2_p[histIter]->Add(histPP_p[histIter], -1);
  }

  Int_t panels;
  Int_t ppStart;
  const char* ppChar[2];

  if(!strcmp(CNCR, "")){
    panels = 10;
    ppStart = 7;
    ppChar[0] = "50-100%";
    ppChar[1] = "30-50%";
  }
  else{
    panels = 6;
    ppStart = 5;
    ppChar[0] = "30-100%";
    ppChar[1] = "0-30%";
  }

  if(strcmp(CNCR, "") != 0){
    for(Int_t coordIter = 0; coordIter < 4; coordIter++){
      overCoord[coordIter] -= .04;
    }
  }

  profPanel_p->cd(ppStart);
  makeHistForPtStack(hist1_p, ppStart, Tight, CNCR, isPercent);
  drawFullStack(hist1_p, 0, 24, 0, true, CNCR, isPercent);
  label1_p->DrawLatex(.22, overCoord[2], overLabel[3]);
  label1_p->DrawLatex(.22, overCoord[3], ppChar[0]);
  label1_p->DrawLatex(.24, .38, "p_{T,1}>120 GeV/c");
  label1_p->DrawLatex(.24, .28, "p_{T,2}>50 GeV/c");

  zeroLine_p->Draw();

  profPanel_p->cd(ppStart+1);
  makeHistForPtStack(hist2_p, ppStart+1, Tight, CNCR, isPercent);
  drawFullStack(hist2_p, 0, 24, 0, true, CNCR, isPercent);
  label1_p->DrawLatex(.05, overCoord[2], overLabel[3]);
  label1_p->DrawLatex(.05, overCoord[3], ppChar[1]);

  if(!strcmp(CNCR, "R") || !strcmp(CNCR, "RD") || !strcmp(CNCR, "RU"))
    label1_p->DrawLatex(.05, .38, "|#eta|_{1},|#eta|_{2}<0.5");
  else
    label1_p->DrawLatex(.05, .38, "|#eta|_{1},|#eta|_{2}<1.6");

  label1_p->DrawLatex(.05, .28, "#Delta#phi_{1,2}>5#pi/6");

  zeroLine_p->Draw();

  legB_p->AddEntry(hist2_p[5], "PbPb - pp", "p");

  if(!strcmp(CNCR, "")){
    profPanel_p->cd(9);
    makeHistForPtStack(hist3_p, 9, Tight, CNCR, isPercent);
    drawFullStack(hist3_p, 0, 24, 0, true, CNCR, isPercent);
    label1_p->DrawLatex(.05, overCoord[2], overLabel[3]);
    label1_p->DrawLatex(.05, overCoord[3], "10-30%");
    label1_p->DrawLatex(.05, .38, "anti-k_{T} Calo R = 0.3");

    zeroLine_p->Draw();

    profPanel_p->cd(10);
    makeHistForPtStack(hist4_p, 10, Tight, CNCR, isPercent);
    drawFullStack(hist4_p, 0, 24, 0, true, CNCR, isPercent);
    label1_p->DrawLatex(.05, overCoord[2], overLabel[3]);
    label1_p->DrawLatex(.05, overCoord[3], "0-10%");

    zeroLine_p->Draw();
  }

  for(Int_t panelIter = 0; panelIter < panels; panelIter++){
    profPanel_p->cd(panelIter+1)->RedrawAxis();
  }

  TFile* outFile_p = new TFile(outName, "UPDATE");
  profPanel_p->Write();
  claverCanvasSaving(profPanel_p, Form("pdfDir/%s%sImbProjA%s%s%sPTStack_%s", gorr, algType[setNum], CNCR, Corr, Tight, fileTagPbPb), "pdf");
  outFile_p->Close();

  delete outFile_p;
  delete label1_p;
  delete zeroLine_p;
  delete legB_p;
  delete legA_p;
  delete profPanel_p;

  histPbPbFile_p->Close();

  delete histPbPbFile_p;

  return;
}


void makeDiJetPlots(const char* filePbPbName, const char* fileTagPbPb, const char* outName, Bool_t montecarlo = false, const char* filePPName = "", const char* fileTagPP = "", Bool_t isPercent = false)
{
  TH1::SetDefaultSumw2();

  Int_t jetAlgMax = 2;
  
  if(montecarlo)
    jetAlgMax = 3;
  
  const char* corr[2] = {"", "Corr"};
  const char* CNCR[6] = {"", "C", "NC", "R", "RD", "RU"};
  const char* Tight[2] = {"", "Tight"};

  for(Int_t algIter = 0; algIter < jetAlgMax; algIter++){
    for(Int_t tightIter = 0; tightIter < 2; tightIter++){
      for(Int_t corrIter = 0; corrIter < 2; corrIter++){
	for(Int_t CNCRIter = 0; CNCRIter < 6; CNCRIter++){
	  if((CNCRIter == 3 || CNCRIter == 4 || CNCRIter == 5) && tightIter == 1)
	    continue;

	  makeImbPtStack(filePbPbName, fileTagPbPb, outName, "r", algIter, corr[corrIter], CNCR[CNCRIter], montecarlo, filePPName, fileTagPP, Tight[tightIter], isPercent);

	  if(montecarlo && corrIter > 0)
	    makeImbPtStack(filePbPbName, fileTagPbPb, outName, "g", algIter, corr[corrIter], CNCR[CNCRIter], montecarlo, filePPName, fileTagPP, Tight[tightIter], isPercent);
	}
      }
    }
  }

  return;
}

