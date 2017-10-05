#include <TSystem.h>
#include <TMath.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
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

void gethist_dist(){

  //TH1::SetDefaultSumw2();
  gROOT->Macro("~/JpsiStyle.C");

  double ptarr[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 40};
  const int Nptarr = sizeof(ptarr)/sizeof(double);
  double raparr[] = {-2.4,-2.2,-2.,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.,2.2,2.4};
  const int Nraparr = sizeof(raparr)/sizeof(double);

  TFile *f1 = new TFile("OniaTrees/run285684_All.root");
  TFile *f2 = new TFile("OniaTrees/run285726_All.root");
  TFile *f3 = new TFile("OniaTrees/run285759_All.root");
  TTree *hlttree1 = (TTree*)f1->Get("hltbitanalysis/HltTree");
  TTree *hlttree2 = (TTree*)f2->Get("hltbitanalysis/HltTree");
  TTree *hlttree3 = (TTree*)f3->Get("hltbitanalysis/HltTree");

  /////////////////hltanalyzer/HltTree/////////////////
  TH1D *f1h1 = new TH1D("f1h1","f1h1",Nraparr-1,raparr);
  TH1D *f1h2 = new TH1D("f1h2","f1h2",Nraparr-1,raparr);
  TH1D *f1h3 = new TH1D("f1h3","f1h3",Nraparr-1,raparr);
  TH1D *f1h4 = new TH1D("f1h4","f1h4",Nraparr-1,raparr);
  TH1D *f1h5 = new TH1D("f1h5","f1h5",Nraparr-1,raparr);
  TH1D *f2h1 = new TH1D("f2h1","f2h1",Nraparr-1,raparr);
  TH1D *f2h2 = new TH1D("f2h2","f2h2",Nraparr-1,raparr);
  TH1D *f2h3 = new TH1D("f2h3","f2h3",Nraparr-1,raparr);
  TH1D *f2h4 = new TH1D("f2h4","f2h4",Nraparr-1,raparr);
  TH1D *f2h5 = new TH1D("f2h5","f2h5",Nraparr-1,raparr);
  TH1D *f3h1 = new TH1D("f3h1","f3h1",Nraparr-1,raparr);
  TH1D *f3h2 = new TH1D("f3h2","f3h2",Nraparr-1,raparr);
  TH1D *f3h3 = new TH1D("f3h3","f3h3",Nraparr-1,raparr);
  TH1D *f3h4 = new TH1D("f3h4","f3h4",Nraparr-1,raparr);
  TH1D *f3h5 = new TH1D("f3h5","f3h5",Nraparr-1,raparr);

     //BX:1-270
     //BX:1162-1431
     //BX:1871-2323
     //BX:2323-2673

  hlttree1->Draw("L1Stage2MuonPt>>f1h5","","pe");
  hlttree2->Draw("L1Stage2MuonPt>>f2h5","","pe");
  hlttree3->Draw("L1Stage2MuonPt>>f3h5","","pe");
  hlttree1->Draw("L1Stage2MuonEta>>f1h1",
      "(HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1>0)&&(L1Stage2MuonBx==0)&&(Bx>1 && Bx<270)","pe");
  hlttree1->Draw("L1Stage2MuonEta>>f1h2",
      "(HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1>0)&&(L1Stage2MuonBx==0)&&(Bx>1162 && Bx<1431)","pe");
  hlttree1->Draw("L1Stage2MuonEta>>f1h3",
      "(HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1>0)&&(L1Stage2MuonBx==0)&&(Bx>1871 && Bx<2323)","pe");
  hlttree1->Draw("L1Stage2MuonEta>>f1h4",
      "(HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1>0)&&(L1Stage2MuonBx==0)&&(Bx>2323 && Bx<2673)","pe");
  hlttree2->Draw("L1Stage2MuonEta>>f2h1",
      "(HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1>0)&&(L1Stage2MuonBx==0)&&(Bx>1 && Bx<270)","pe");
  hlttree2->Draw("L1Stage2MuonEta>>f2h2",
      "(HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1>0)&&(L1Stage2MuonBx==0)&&(Bx>1162 && Bx<1431)","pe");
  hlttree2->Draw("L1Stage2MuonEta>>f2h3",
      "(HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1>0)&&(L1Stage2MuonBx==0)&&(Bx>1871 && Bx<2323)","pe");
  hlttree2->Draw("L1Stage2MuonEta>>f2h4",
      "(HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1>0)&&(L1Stage2MuonBx==0)&&(Bx>2323 && Bx<2673)","pe");
  hlttree3->Draw("L1Stage2MuonEta>>f3h1",
      "(HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1>0)&&(L1Stage2MuonBx==0)&&(Bx>1 && Bx<270)","pe");
  hlttree3->Draw("L1Stage2MuonEta>>f3h2",
      "(HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1>0)&&(L1Stage2MuonBx==0)&&(Bx>1162 && Bx<1431)","pe");
  hlttree3->Draw("L1Stage2MuonEta>>f3h3",
      "(HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1>0)&&(L1Stage2MuonBx==0)&&(Bx>1871 && Bx<2323)","pe");
  hlttree3->Draw("L1Stage2MuonEta>>f3h4",
      "(HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1>0)&&(L1Stage2MuonBx==0)&&(Bx>2323 && Bx<2673)","pe");

  TCanvas *c1 = new TCanvas("c1","c1",800,400);
  TCanvas *c2 = new TCanvas("c2","c2",800,400);
  TCanvas *c3 = new TCanvas("c3","c3",800,400);
  gStyle->SetOptStat(0);
  gStyle->SetFillColor(0);

  TLegend *leg = new TLegend(0.175,0.675,0.55,0.925);
  leg->AddEntry(f1h1,"Bx>1 && Bx<270","lp");
  leg->AddEntry(f1h2,"Bx>1162 && Bx<1431","lp");
  leg->AddEntry(f1h3,"Bx>1871 && Bx<2323","lp");
  leg->AddEntry(f1h4,"Bx>2323 && Bx<2673","lp");
  SetLegendStyle(leg);
  leg->SetFillStyle(0);

  SetHistStyle(f1h1,0,1,0.,0.35);
  SetHistStyle(f1h2,1,2,0.,0.35);
  SetHistStyle(f1h3,2,3,0.,0.35);
  SetHistStyle(f1h4,3,4,0.,0.35);
  SetHistStyle(f2h1,0,1,0.,0.35);
  SetHistStyle(f2h2,1,2,0.,0.35);
  SetHistStyle(f2h3,2,3,0.,0.35);
  SetHistStyle(f2h4,3,4,0.,0.35);
  SetHistStyle(f3h1,0,1,0.,0.35);
  SetHistStyle(f3h2,1,2,0.,0.35);
  SetHistStyle(f3h3,2,3,0.,0.35);
  SetHistStyle(f3h4,3,4,0.,0.35);

  c1->cd();
  f1h1->DrawNormalized("pe");
  f1h2->DrawNormalized("sames");
  f1h3->DrawNormalized("sames");
  f1h4->DrawNormalized("sames");
  leg->Draw("sames");
  f1h4->SetXTitle("p_{T}(GeV/c)");
  f1h4->GetYaxis()->SetRangeUser(0,0.7);
  gPad->Update(); c1->Update();
  c1->SaveAs("test1.png");

  c2->cd();
  f2h1->DrawNormalized("pe");
  f2h2->DrawNormalized("sames");
  f2h3->DrawNormalized("sames");
  f2h4->DrawNormalized("sames");
  leg->Draw("sames");
  f2h4->SetXTitle("p_{T}(GeV/c)");
  f2h4->GetYaxis()->SetRangeUser(0,0.7);
  gPad->Update(); c2->Update();
  c2->SaveAs("test2.png");

  c3->cd();
  f3h1->DrawNormalized("pe");
  f3h2->DrawNormalized("sames");
  f3h3->DrawNormalized("sames");
  f3h4->DrawNormalized("sames");
  leg->Draw("sames");
  f3h4->SetXTitle("p_{T}(GeV/c)");
  f3h4->GetYaxis()->SetRangeUser(0,0.7);
  gPad->Update(); c3->Update();
  c3->SaveAs("test3.png");

}
