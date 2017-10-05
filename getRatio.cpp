#include <TSystem.h>
#include <TMath.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TPaveStats.h>
#include <TTree.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <iostream>
#include <string>
#include "JpsiFunc.h"

void getRatio(){

  string date="1122";
  string dataset="Data";
  string vername="PromptReco_5TeV_Double_All";
  string refname="PromptReco_5TeV_Double_All";
  const int Ntrig = 19;
  string trigname[Ntrig]={
    "HLT_PAL1DoubleMuOpen_v1",
    "HLT_PAL1DoubleMuOpen_OS_v1",
    "HLT_PAL1DoubleMuOpen_SS_v1",
    "HLT_PAL1DoubleMu0_v1",
    "HLT_PAL1DoubleMu0_MGT1_v1",
    "HLT_PAL1DoubleMu0_HighQ_v1",
    "HLT_PAL2DoubleMu0_v1",
    "HLT_PAL3DoubleMu0_v1",
    "HLT_PAL3DoubleMu0_HIon_v1",
    "HLT_PAL1DoubleMu10_v1",
    "HLT_PAL2DoubleMu10_v1",
    "HLT_PAL3DoubleMu10_v1",
    "HLT_PAL2Mu12_v1",
    "HLT_PAL2Mu15_v1",
    "HLT_PAL3Mu3_v1",
    "HLT_PAL3Mu5_v3",
    "HLT_PAL3Mu7_v1",
    "HLT_PAL3Mu12_v1",
    "HLT_PAL3Mu15_v1"
  };

  f1->cd();
  double *xValuesA[i] = divide_nume_p[i]_by_deno_p[i]->GetX();
  double *yValuesA[i] = divide_nume_p[i]_by_deno_p[i]->GetY();
  f2->cd();
  double *xValuesB[i] = divide_nume_p[i]_by_deno_p[i]->GetX();
  double *yValuesB[i] = divide_nume_p[i]_by_deno_p[i]->GetY();

  string save[]={"MuonL3_Pt", "MuonL3_Eta"};
  string file1[]= {
    form("figs/%s/%s/%s/h_pt_efficiency_%s.root",date.c_str(),dataset.c_str(),vername.c_str(),trigname[i].c_str())
  };
  string file2[]= {
    form("figs/%s/%s/%s/h_pt_efficiency_%s.root",date.c_str(),dataset.c_str(),refname.c_str(),trigname[i].c_str())
  };

  //version histo name
  string h1name[]={ "eff14p", "eff14e"};
  //reference histo name
  string h2name[]={ "eff14p", "eff14e"};

  TLine *lptx = new TLine(3.,0, 3.,1);
  TLine *lpty = new TLine(0,1, 30,1);
  TLine *letay = new TLine(-2.4,1,2.4,1);

  lptx->SetLineStyle(2);
  lpty->SetLineStyle(2);
  letay->SetLineStyle(2);
  lptx->SetLineWidth(2);
  lpty->SetLineWidth(2);
  letay->SetLineWidth(2);
  lptx->SetLineColor(kRed+2);

  for (int i=0; i<Ntrig; i++){

    TFile *f1 = new TFile(file1[i].c_str()); TH1F *h1 = (TH1F*)f1->Get(h1name[i].c_str()); TH1F *h3=(TH1F*)h1->Clone("h3"); 
    TFile *f2 = new TFile(file2[i].c_str()); TH1F *h2 = (TH1F*)f2->Get(h2name[i].c_str());

    TCanvas *c1 = new TCanvas(save[i].c_str(),save[i].c_str(),800,400); 
    gStyle->SetOptStat(0);
    c1->Divide(2,1); 
    c1->cd(1); h1->Draw(); 
    if(i==0) {lptx->Draw("sames"); lpty->Draw("sames");}
    if(i==1) letay->Draw("sames");
    h1->SetTitle(Form("%s",save[i].c_str())); 
    h1->SetMarkerColor(kRed+2); 
    h1->SetLineColor(kRed+2); //h1->SetMarkerStyle(8); 
    //h2->Draw("sames"); h2->SetMarkerColor(kBlue+2); h2->SetLineColor(kBlue+2); h2->SetMarkerStyle(33);
    //h3->Draw("sames"); h3->SetMarkerColor(kBlack); h3->SetLineColor(kBlack); h3->SetTitle(Form("%s",save[i].c_str()));//h3->SetMarkerStyle(25); h3->SetMarkerSize(1.2);
    c1->cd(1); h2->Draw("sames"); h2->SetTitle(Form("%s_Online_vs_Offline",save[i].c_str())); h2->SetMarkerStyle(33); h2->SetMarkerColor(kBlue+2); h2->SetLineColor(kBlue+2);
    h2->GetYaxis()->SetRangeUser(0,1.6);
    c1->cd(2); h3->Divide(h2); h3->Draw("pe"); 
    if(i==0) {lptx->Draw("sames"); lpty->Draw("sames");}
    if(i==1) letay->Draw("sames");
    h3->SetTitle("Ratio Offline/Online"); 
    h3->SetYTitle("Ratio"); h3->SetMarkerColor(kBlack); h3->SetLineColor(kBlack); 
    h3->GetYaxis()->SetRangeUser(0.8,1.2);
    c1->SaveAs(Form("test_%s.png",save[i].c_str()));
  }
}
