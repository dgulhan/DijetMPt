//=============================================                                                
// Author: Chris McGinn                                                                         
//                                                                                            
// DiJet Plotter                                                              
//                                                                                            
//=============================================                                               

#include "commonUtility.h"
#include "TDatime.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"

TFile* histFile_p = 0;

TFile* plotFile_p = 0;

const char* algType[3] = {"PuCalo", "VsCalo", "T"};

Double_t quadSum(Double_t one, Double_t two)
{
  Double_t err = TMath::Sqrt(one*one + two*two);
  return err;
}


void drawPatch(float x1, float y1, float x2, float y2){
  TLegend *t1=new TLegend(x1,y1,x2,y2);
  t1->SetFillColor(kWhite);
  t1->SetBorderSize(0);
  t1->SetFillStyle(1001);
  t1->Draw("");
}


void niceTH1(TH1F* uglyTH1, float max , float min, float ndivX, float ndivY)
{
  handsomeTH1N(uglyTH1);
  uglyTH1->SetMaximum(max);
  uglyTH1->SetMinimum(min);
  uglyTH1->SetNdivisions(ndivX);
  uglyTH1->SetNdivisions(ndivY, "Y");
}


void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns, const Int_t rows, const Float_t leftOffset, const Float_t bottomOffset, const Float_t leftMargin, const Float_t bottomMargin, const Float_t edge, const char* CNC = ""){
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
        if(strcmp(CNC, "") ==0)
          pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j]*0.933,Xup[i],Yup[j]);
        else
          pad[i][j] = new TPad(padName.Data(),padName.Data(),Xlow[i],Ylow[j]*0.934,Xup[i],Yup[j]*.9985);
      }
      else if(i==1 && j==1){
        if(strcmp(CNC, "") == 0)
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
        if(strcmp(CNC, "") == 0)
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


void makeHistForPtStack(TH1F* h_p[6], Int_t pos = 4, const char* CNC = "")
{
  for(Int_t iter = 0; iter < 4; iter++){
    h_p[3]->SetBinContent(iter + 1, sumYForPTStack(h_p[3]->GetBinContent(iter+1), h_p[4]->GetBinContent(iter+1)));

    h_p[2]->SetBinContent(iter + 1, sumYForPTStack(h_p[2]->GetBinContent(iter+1), h_p[3]->GetBinContent(iter+1), h_p[4]->GetBinContent(iter+1)));

    h_p[1]->SetBinContent(iter + 1, sumYForPTStack(h_p[1]->GetBinContent(iter+1), h_p[2]->GetBinContent(iter+1), h_p[3]->GetBinContent(iter+1), h_p[4]->GetBinContent(iter+1)));

    h_p[0]->SetBinContent(iter + 1, sumYForPTStack(h_p[0]->GetBinContent(iter+1), h_p[1]->GetBinContent(iter+1), h_p[2]->GetBinContent(iter+1), h_p[3]->GetBinContent(iter+1), h_p[4]->GetBinContent(iter+1)));
  }

  h_p[4]->SetXTitle("A_{J}");
  h_p[3]->SetXTitle("A_{J}");
  h_p[2]->SetXTitle("A_{J}");
  h_p[1]->SetXTitle("A_{J}");
  h_p[0]->SetXTitle("A_{J}");
  h_p[5]->SetXTitle("A_{J}");

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

void drawHistToPTStack(TH1F* drawHist_p, Int_t color, const char* drawOpt)
{
  drawHist_p->SetFillColor(color);
  drawHist_p->SetMarkerStyle(6);
  drawHist_p->SetMarkerSize(.5);
  drawHist_p->DrawCopy(drawOpt);
  drawHist_p->DrawCopy("E1 SAME");
}

void drawFullStack(TH1F* h_p[6], Int_t color, Int_t style, TLegend* leg_p = 0)
{
  drawHistToPTStack(h_p[0], kBlue - 9, "E1 HIST");
  drawHistToPTStack(h_p[1], kYellow - 9, "E1 HIST SAME");
  drawHistToPTStack(h_p[2], kOrange + 1, "E1 HIST SAME");
  drawHistToPTStack(h_p[3], kGreen + 3, "E1 HIST SAME");
  drawHistToPTStack(h_p[4], kRed + 1, "E1 HIST SAME");

  h_p[5]->SetFillColor(color);
  h_p[5]->SetMarkerStyle(style);

  h_p[5]->DrawCopy("SAME E1");

  if(!leg_p){
    leg_p->AddEntry(h_p[0], "0.5 - 1.0", "F");
    leg_p->AddEntry(h_p[1], "1.0 - 2.0", "F");
    leg_p->AddEntry(h_p[2], "2.0 - 4.0", "F");
    leg_p->AddEntry(h_p[3], "4.0 - 8.0", "F");
    leg_p->AddEntry(h_p[4], "8.0 - 300.0", "F");
  }

  return;
}


void makeImbAsymmPtStack(const char* filePbPbName, const char* fileTagPbPb, const char* outName, const char* gorr, Int_t setNum, const char* Corr = "", const char* CNC = "", Bool_t montecarlo = false)
{
  TFile histPbPbFile_p = new TFile(filePbPbName, "READ");

  const char* mcLabel[4] = {"PYTHIA", "PYTHIA + HYDJET", "PYTHIA + HYDJET", "(P + H) - P"};
  const char* dataLabel[4] = {"pp 5.3 pb^{-1}", "PbPb 150 #mub^{-1}", "PbPb", "PbPb - pp"};

  const char* overLabel[4];
  Float_t overCoord[4] = {.84, .76, .90, .82};

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
  TH1F* histpp_p[6];

  Float_t binArrayX[5] = {.00, .11, .22, .33, .499};

  const char* FPT[6] = {"0_1", "1_2", "2_4", "4_8", "8_100", "F"};

  for(Int_t histIter = 0; histIter < 6; histIter++){
    if(!strcmp(CNC, "")){
      hist1_p[histIter] = (TH1F*)histPbPbFile_p->Get("%s%sImbAsymmProjA%s%s%s_50100_%s_h", gorr, algType[setNum], CNC, FPT[histIter], Corr, fileTag);
      hist2_p[histIter] = (TH1F*)histPbPbFile_p->Get("%s%sImbAsymmProjA%s%s%s_3050_%s_h", gorr, algType[setNum], CNC, FPT[histIter], Corr, fileTag);
      hist3_p[histIter] = (TH1F*)histPbPbFile_p->Get("%s%sImbAsymmProjA%s%s%s_1030_%s_h", gorr, algType[setNum], CNC, FPT[histIter], Corr, fileTag);
      hist4_p[histIter] = (TH1F*)histPbPbFile_p->Get("%s%sImbAsymmProjA%s%s%s_010_%s_h", gorr, algType[setNum], CNC, FPT[histIter], Corr, fileTag);
    }
    else{
      hist1_p[histIter] = (TH1F*)histPbPbFile_p->Get("%s%sImbAsymmProjA%s%s%s_30100_%s_h", gorr, algType[setNum], CNC, FPT[histIter], Corr, fileTag);
      hist2_p[histIter] = (TH1F*)histPbPbFile_p->Get("%s%sImbAsymmProjA%s%s%s_030_%s_h", gorr, algType[setNum], CNC, FPT[histIter], Corr, fileTag);
      niceTH1(hist3_p[histIter], 59.999, -59.999, 505, 406);
      niceTH1(hist4_p[histIter], 59.999, -59.999, 505, 406);
    }

    niceTH1(hist1_p[histIter], 59.999, -59.999, 505, 406);
    niceTH1(hist2_p[histIter], 59.999, -59.999, 505, 406);
  }

  //Draw first PbPb panel

  makeHistForPtStack(hist1_p, 2, CNC);

  TCanvas* profPanel_p;

  if(!strcmp(CNC, "")){
    profPanel_p = new TCanvas(Form("%s%sImbAsymm%s%s%s%sPTStackPP_%s_%s_c", gorr, algType[setNum], Tight, perpProj, CNC, Corr, GLN, fileTag1), Form("%s%sImbAsymm%s%s%s%sPTStackPP_%s_%s_c", gorr, algType[setNum], Tight, perpProj, CNC, Corr, GLN, fileTag1), 5*300, 700);
    makeMultiPanelCanvas(profPanel_p, 5, 2, 0.0, 0.0, 0.2, 0.2, 0.05);
    std::cout << "FivePanel Init" << std::endl;
  }
  else{
    profPanel_p = new TCanvas(Form("%s%sImbAsymm%s%s%s%sPTStackPP_%s_%s_c", gorr, algType[setNum], Tight, perpProj, CNC, Corr, GLN, fileTag1), Form("%s%sImbAsymm%s%s%s%sPTStackPP_%s_%s_c", gorr, algType[setNum], Tight, perpProj, CNC, Corr, GLN, fileTag1), 300*3, 700);
    makeMultiPanelCanvas(profPanel_p, 3, 2, 0.0, 0.0, 0.2, 0.2, 0.01, CNC);
    std::cout << "ThreePanel Init" << std::endl;
  }

  //Make legend

  TLegend* legA_p = new TLegend(0.25, 0.18, 0.99, 0.88);
  TLegend* legB_p = new TLegend(0.45, 0.01, 0.99, 0.80);

  legA_p->SetFillColor(0);
  legA_p->SetTextFont(43);
  legA_p->SetTextSizePixels(28);
  legA_p->SetBorderSize(0);

  legB_p->SetFillColor(0);
  legB_p->SetTextFont(43);
  legB_p->SetTextSizePixels(22);
  legB_p->SetBorderSize(0);

  profPanel_p->cd(2);

  drawFullStack(hist1_p, 0, 28, legA_p);

  TLine* zeroLine_p = new TLine(0., 0., 0.5, 0.);
  zeroLine_p->SetLineColor(1);
  zeroLine_p->SetLineStyle(1);
  zeroLine_p->Draw();

  TLatex* label1_p = new TLatex();
  label1_p->SetNDC();
  label1_p->SetTextFont(43);
  label1_p->SetTextSizePixels(28);

  label1_p->DrawLatex(.05, overCoord[0], Form("%s", overLabel[1]));
  if(strcmp("", CNC) == 0)
    label1_p->DrawLatex(.05, overCoord[1], "50-100%");
  else
    label1_p->DrawLatex(.05, overCoord[1], "30-100%");

  label1_p->DrawLatex(.05, .05, "#sqrt{s_{NN}} = 2.76 TeV");

  profPanel_p->cd(3);

  //Draw second PbPb hist

  makeHistForPtStack(hist2_p, 3, CNC);

  drawFullStack(hist2_p, 0, 28);

  zeroLine_p->Draw();

  label1_p->DrawLatex(.05, overCoord[0], Form("%s", overLabel[2]));
  if(!strcmp("", CNC))
    label1_p->DrawLatex(.05, overCoord[1], "30-50%");
  else
    label1_p->DrawLatex(.05, overCoord[1], "0-30%");

  //Draw third and fourth PbPb panels, if applicable

  if(!strcmp(CNC, "")){
    profPanel_p->cd(4);

    makeHistForPtStack(hist3_p, 4, CNC);

    hist3_p[5]->SetMarkerStyle(28);
    hist3_p[5]->SetFillColor(0);

    drawFullStack(hist3_p, 0, 28);

    zeroLine_p->Draw();

    label1_p->DrawLatex(.05, overCoord[0], Form("%s", overLabel[2]));
    label1_p->DrawLatex(.05, overCoord[1], "10-30%");

    profPanel_p->cd(5);

    makeHistForPtStack(hist4_p, 5, CNC);

    hist4_p[5]->SetMarkerStyle(28);
    hist4_p[5]->SetFillColor(0);

    drawFullStack(hist4_p, 0, 28);

    zeroLine_p->Draw();

    label1_p->DrawLatex(.05, overCoord[0], Form("%s", overLabel[2]));
    label1_p->DrawLatex(.05, overCoord[1], "0-10%");
  }

  TFile* outFile_p = new TFile(outName, "RECREATE");
  profPanel_p->Write();
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


void makeDiJetPlots(const char* filePbPbName, const char* fileTagPbPb, const char* outName, Bool_t montecarlo = false)
{
  TH1::SetDefaultSumw2();

  Int_t jetAlgMax = 2;
  
  if(montecarlo)
    jetAlgMax = 3;
  
  const char* corr[2] = {"", "Corr"};
  const char* CNC[3] = {"", "C", "NC"};
  
  for(Int_t algIter = 0; algIter < jetAlgMax; algIter++){
    for(Int_t corrIter = 0; corrIter < 2; corrIter++){
      for(Int_t CNCIter = 0; CNCIter < 3; CNCIter++){
	makeImbAsymmPtStack(filePbPbName, fileTagPbPb, outName, "r", setNum, corr[corrIter], CNC[CNCIter], montecarlo);
	
	if(montecarlo && corrIter > 0)
	  makeImbAsymmPtStack(filePbPbName, fileTagPbPb, outName, "g", setNum, corr[corrIter], CNC[CNCIter], montecarlo);
      }
    }
  }
  
  return;
}

