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
#include <TGraphAsymmErrors.h>
#include <iostream>
#include <string>
#include "JpsiFunc.h"

// use tnp true and track true for mutrk triggers only. Usetnp false and isTrk false for SingleMu triggers except track triggers. Use tnp true and is trk false for DoubleMu triggers only.
bool usetnp = true;
bool isTrk = false;

void getCorrectedEffErr(const int nbins, TH1D *hrec, TH1D *hgen, TH1D *heff) {
  for (int a=0; a<nbins; a++) {
    double genInt = hgen->GetBinContent(a+1);
    double genErr = hgen->GetBinError(a+1);
    double recInt = hrec->GetBinContent(a+1);
    double recErr = hrec->GetBinError(a+1);
    double eff = recInt / genInt;

    double tmpErrGen1 = TMath::Power(eff,2) / TMath::Power(genInt,2);
    double tmpErrRec1 = TMath::Power(recErr,2);
    double tmpErr1 = tmpErrGen1 * tmpErrRec1;

    double tmpErrGen2 = TMath::Power(1-eff,2) / TMath::Power(genInt,2);
    double tmpErrRec2 = TMath::Abs(TMath::Power(genErr,2) - TMath::Power(recErr,2));
    double tmpErr2 = tmpErrGen2 * tmpErrRec2;
    double effErr = TMath::Sqrt(tmpErr1 + tmpErr2);

    if (genInt == 0) {
      heff->SetBinContent(a+1, 0);
      heff->SetBinError(a+1, 0);
    } else {
      heff->SetBinContent(a+1, eff);
      heff->SetBinError(a+1, effErr);
    }
  }
}

bool AcceptanceCut (TLorentzVector* Muon)
{
  return ((fabs(Muon->Eta())<1.2&&Muon->Pt()>=3.3)||(1.2<=fabs(Muon->Eta())&&fabs(Muon->Eta())<2.1&&Muon->Pt()>=(4.0-1.1*fabs(Muon->Eta())))||(2.1<=fabs(Muon->Eta())&&fabs(Muon->Eta())<2.4&&Muon->Pt()>=1.3));
}; 

//main
void getEff1117_Santona(){


  TChain *fcha = new TChain("hionia/myTree");

  //MC
  // For 2017 pp ref run, Prompt JPsi MC sample (test)
  //fcha->Add("/afs/cern.ch/user/t/twang/public/ForAndre/OniaForest.root");

  // For 2017 pp ref run, JPsi Gun MC
  //fcha->Add("./InputFiles/Pythia8_JPsiGun_pp_2017pp502_Onia_20171031_merge.root");
  // For 2017 pp ref run, Mu Gun MC
  //fcha->Add("./InputFiles/Pythia8_MuGun_pp_2017pp502_Onia_20171031_merge.root");
  // For 2017 pp ref run, Pr JPsi MC (larger)
  //fcha->Add("./InputFiles/Pythia8_PrJPsi_pp_2017pp502_Onia_20171031_merge.root");
  // For 2017 pp ref run, NonPr JPsi MC
  //fcha->Add("./InputFiles/Pythia8_NonPrJPsi_pp_2017pp502_Onia_20171031_merge.root");
  // For 2017 pp ref run, Z MC
  //fcha->Add("./InputFiles/Pythia8_Zm10m10_pp_2017pp502_Onia_20171031_merge.root");

  // 11/17 onia trees
  fcha->Add("/eos/cms/store/group/phys_heavyions/jaebeom/OniaTree_PPref_SingleMuTnP/SingleMuonTnP/OniaTree_PPref_SingleMuTnP/171116_031052/0000/OniaForest_321.root");
  fcha->Add("/eos/cms/store/group/phys_heavyions/jaebeom/OniaTree_PPref_SingleMuTnP/SingleMuonTnP/OniaTree_PPref_SingleMuTnP/171116_031052/0000/OniaForest_322.root");
  fcha->Add("/eos/cms/store/group/phys_heavyions/jaebeom/OniaTree_PPref_SingleMuTnP/SingleMuonTnP/OniaTree_PPref_SingleMuTnP/171116_031052/0000/OniaForest_323.root");




  double ptmin = 0;
  double ptmax = 50;
  double etamin = -2.4;
  double etamax = 2.4;
  double massmin = 3;
  double massmax = 3.2;  // tight mass cut
  int minpTtrigfordeno = 18; // L3Mu3
  int L3SingleMuTrig = 18;

  bool isZ = false;
  if(isZ==1){
   massmin = 70;
   massmax = 110;
   minpTtrigfordeno = 13;
  }

  string date="111817";
  string dataset="Data";
  //string dataset="MC";

  //Santona // For 2017 pp ref run
  //previous smaller sample //string vername="PromptJPsiMC";
  //string vername="JPsiGun";
  //string vername="MuGun";
  //string vername="PrJPsi";
  //string vername="NonPrJPsi";
  //string vername="Zm10m10";

  string vername="OniaForest_ppRefRunTest";

  //define Trigger
  const int Ntrig = 30;  // was 27 // 18 for 2015
  string trigname[Ntrig]={ //was Ntrig+1
  // 2017 pp Ref Data. 30 total
  //Double
  "HLT_HIL1DoubleMuOpen_v1", // 0
  "HLT_HIL1DoubleMuOpen_OS_v1", // 1
  "HLT_HIL1DoubleMuOpen_SS_v1", // 2
  "HLT_HIL1DoubleMu0_v1", //3
  "HLT_HIL1DoubleMu0_HighQ_v1", //4
  "HLT_HIL1DoubleMu10_v1", //5
  "HLT_HIL2DoubleMu0_v1", // 6
  "HLT_HIL2DoubleMu10_v1", //7
  "HLT_HIL3DoubleMu0_v1", //8
  "HLT_HIL3DoubleMu10_v1", //9
  // Single
  "HLT_HIL1Mu12_v1", //10
  "HLT_HIL1Mu16_v1", //11
  "HLT_HIL2Mu3_NHitQ10_v1", //12
  "HLT_HIL2Mu5_NHitQ10_v1", //13
  "HLT_HIL2Mu7_v1", //14
  "HLT_HIL2Mu12_v1", //15
  "HLT_HIL2Mu15_v1", //16
  "HLT_HIL2Mu20_v1", //17
  "HLT_HIL3Mu3_v1", //18
  "HLT_HIL3Mu3_NHitQ10_v1", //19
  "HLT_HIL3Mu5_v1", //20
  "HLT_HIL3Mu5_NHitQ10_v1", //21
  "HLT_HIL3Mu7_v1", //22
  "HLT_HIL3Mu12_v1", //23
  "HLT_HIL3Mu15_v1", //24
  "HLT_HIL3Mu20_v1", //25
  "HLT_HIL3Mu3_Track1_Jpsi_v1", //26
  "HLT_HIL3Mu5_Track1_Jpsi_v1", //27
  "HLT_HIL3Mu3_Track1_v1", //28
  "HLT_HIL3Mu5_Track1_v1" // 29
  };

/*    // For 2017 (Latest) 33, 30 needed (was 28) total
    //Double
    "HLT_HIL1DoubleMuOpen_v1",//0
    "HLT_HIL1DoubleMuOpen_OS_v1",//1
    "HLT_HIL1DoubleMuOpen_SS_v1",//2
    "HLT_HIL1DoubleMu0_v1",//3
    "HLT_HIL1DoubleMu0_HighQ_v1",//4
    "HLT_HIL1DoubleMu10_v1",//5
    "HLT_HIL2DoubleMu0_v1",//6
    "HLT_HIL2DoubleMu10_v1",//7
    "HLT_HIL3DoubleMu0_v1",//8
    "HLT_HIL3DoubleMu10_v1",//9
    //Single
    "HLT_HIL1Mu3_v1",//10 // not needed
    "HLT_HIL1Mu5_v1",//11 //not needed
    "HLT_HIL1Mu7_v1",//12 //not needed
    "HLT_HIL1Mu12_v1",//13
    "HLT_HIL1Mu16_v1",//14
    "HLT_HIL2Mu3_NHitQ10_v1",//15 // new N Hit requirement
    "HLT_HIL2Mu5_NHitQ10_v1",//16 // new N hit requirement
    "HLT_HIL2Mu7_v1",//17
    "HLT_HIL2Mu12_v1",//18
    "HLT_HIL2Mu15_v1",//19
    "HLT_HIL2Mu20_v1",//20
    "HLT_HIL3Mu3_v1",//21
    "HLT_HIL3Mu5_v1",//22
    "HLT_HIL3Mu7_v1",//23
    "HLT_HIL3Mu12_v1",//24
    "HLT_HIL3Mu15_v1",//25
    "HLT_HIL3Mu20_v1",//26
    "HLT_HIL3Mu3_NHitQ10_v1",//27 //new
    "HLT_HIL3Mu5_NHitQ10_v1",//28 //new
    "HLT_HIL3Mu3_Track1_Jpsi_v1",//29
    "HLT_HIL3Mu5_Track1_Jpsi_v1",//30 //new
    "HLT_HIL3Mu3_Track1_v1",//31 //new
    "HLT_HIL3Mu5_Track1_v1" //32 //new
  };
// */



  //Set Branch
//  Int_t           Centrality;
  UInt_t          runNb;
  //Int_t           Gen_mu_size; //for MC
  Int_t           Reco_mu_size;
  Int_t           Reco_QQ_size;
  Int_t           trigPrescale[1000];
  //TClonesArray    *Gen_mu_4mom; //for MC
  TClonesArray    *Reco_mu_4mom;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_QQ_mupl_4mom;
  TClonesArray    *Reco_QQ_mumi_4mom;
  ULong64_t       Reco_mu_trig[1000];
  ULong64_t       Reco_QQ_trig[1000];
  ULong64_t       Reco_QQ_mupl_trig[1000];
  ULong64_t       Reco_QQ_mumi_trig[1000];
  Int_t           Reco_QQ_sign[1000];
  ULong64_t       HLTriggers;
  Bool_t          Reco_mu_isGoodMuon[1000];
  Int_t           Reco_mu_nPixWMea[1000];
  Int_t           Reco_mu_nTrkWMea[1000];
  Float_t         Reco_mu_dxy[1000];
  Float_t         Reco_mu_dz[1000];
  ULong64_t       Reco_QQ_VtxProb[1000];
  Bool_t          Reco_QQ_mumi_isGoodMuon[1000];
  Int_t           Reco_QQ_mumi_nPixWMea[1000];
  Int_t           Reco_QQ_mumi_nTrkWMea[1000];
  Float_t         Reco_QQ_mumi_dxy[1000];
  Float_t         Reco_QQ_mumi_dz[1000];
  Bool_t          Reco_QQ_mupl_isGoodMuon[1000];
  Int_t           Reco_QQ_mupl_nPixWMea[1000];
  Int_t           Reco_QQ_mupl_nTrkWMea[1000];
  Float_t         Reco_QQ_mupl_dxy[1000];
  Float_t         Reco_QQ_mupl_dz[1000];

//  TBranch        *b_Centrality;
  TBranch        *b_runNb;
  TBranch        *b_trigPrescale;
  //TBranch        *b_Gen_mu_size; //for MC
  TBranch        *b_Reco_mu_size;
  TBranch        *b_Reco_QQ_size;
  //TBranch        *b_Gen_mu_4mom; //for MC
  TBranch        *b_Reco_mu_4mom;
  TBranch        *b_Reco_QQ_4mom;
  TBranch        *b_Reco_QQ_VtxProb;
  TBranch        *b_Reco_QQ_mupl_4mom;
  TBranch        *b_Reco_QQ_mumi_4mom;
  TBranch        *b_Reco_mu_trig;
  TBranch        *b_Reco_QQ_trig;
  TBranch        *b_Reco_QQ_sign;
  TBranch        *b_HLTriggers;
  TBranch        *b_Reco_mu_isGoodMuon;
  TBranch        *b_Reco_mu_nPixWMea;
  TBranch        *b_Reco_mu_nTrkWMea;
  TBranch        *b_Reco_mu_dxy;
  TBranch        *b_Reco_mu_dz;
  TBranch        *b_Reco_QQ_mupl_trig;
  TBranch        *b_Reco_QQ_mumi_trig;
  TBranch        *b_Reco_QQ_mumi_isGoodMuon;
  TBranch        *b_Reco_QQ_mumi_nPixWMea;
  TBranch        *b_Reco_QQ_mumi_nTrkWMea;
  TBranch        *b_Reco_QQ_mumi_dxy;
  TBranch        *b_Reco_QQ_mumi_dz;
  TBranch        *b_Reco_QQ_mupl_isGoodMuon;
  TBranch        *b_Reco_QQ_mupl_nPixWMea;
  TBranch        *b_Reco_QQ_mupl_nTrkWMea;
  TBranch        *b_Reco_QQ_mupl_dxy;
  TBranch        *b_Reco_QQ_mupl_dz;

  //  Gen_mu_4mom = 0; //for MC
  Reco_mu_4mom = 0;
  Reco_QQ_4mom = 0;
  Reco_QQ_mupl_4mom = 0;
  Reco_QQ_mumi_4mom = 0;

//  fcha->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  fcha->SetBranchAddress("runNb", &runNb, &b_runNb);
  fcha->SetBranchAddress("trigPrescale", &trigPrescale, &b_trigPrescale);
  //  fcha->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size); //for MC
  //  fcha->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom); //for MC
  fcha->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  fcha->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  fcha->SetBranchAddress("Reco_QQ_sign", &Reco_QQ_sign, &b_Reco_QQ_sign);
  fcha->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  fcha->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  fcha->SetBranchAddress("Reco_QQ_VtxProb", &Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  fcha->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
  fcha->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
  fcha->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  fcha->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  fcha->SetBranchAddress("Reco_QQ_mupl_trig", Reco_QQ_mupl_trig, &b_Reco_QQ_mupl_trig);
  fcha->SetBranchAddress("Reco_QQ_mumi_trig", Reco_QQ_mumi_trig, &b_Reco_QQ_mumi_trig);
  fcha->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  fcha->SetBranchAddress("Reco_mu_isGoodMuon", Reco_mu_isGoodMuon, &b_Reco_mu_isGoodMuon);
  fcha->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  fcha->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  fcha->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  fcha->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  fcha->SetBranchAddress("Reco_QQ_mumi_isGoodMuon", Reco_QQ_mumi_isGoodMuon, &b_Reco_QQ_mumi_isGoodMuon);
  fcha->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
  fcha->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
  fcha->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
  fcha->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
  fcha->SetBranchAddress("Reco_QQ_mupl_isGoodMuon", Reco_QQ_mupl_isGoodMuon, &b_Reco_QQ_mupl_isGoodMuon);
  fcha->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
  fcha->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
  fcha->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
  fcha->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);

  Int_t nevent = fcha->GetEntries();
  cout<<"nevent: "<<nevent<<endl;

  double ptarr[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 40, 50};
  const int Nptarr = sizeof(ptarr)/sizeof(double);

  double raparr[] = {-2.4,-2.2,-2.,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.,2.2,2.4};
  const int Nraparr = sizeof(raparr)/sizeof(double);

  double phiarr[] = {-3.1,-3.0,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.2,-2.,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,
                      0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1};
  const int Nphiarr = sizeof(phiarr)/sizeof(double);

//  double centarr[] = {0, 20, 40, 60, 80, 100};
//  const int Ncentarr = sizeof(centarr)/sizeof(double);

  double legmin[]={0.875,0.850,0.825,0.800,0.775,0.750,0.725,0.700,0.675};
  const int Nlegmin = sizeof(legmin)/sizeof(double);
  double legmax[]={0.850,0.825,0.800,0.775,0.750,0.725,0.700,0.675,0.650};
  const int Nlegmax = sizeof(legmax)/sizeof(double);

  //Define histograms
  TH1D *test_p[Ntrig];
  TH1D *nume_p[Ntrig];
  TH1D *nume_e[Ntrig];
  TH1D *nume_phi[Ntrig];
  TH1D *deno_p[Ntrig];
  TH1D *deno_e[Ntrig];
  TH1D *deno_phi[Ntrig];
  TGraphAsymmErrors *eff_p[Ntrig];
  TGraphAsymmErrors *eff_e[Ntrig];
  //make histrograms
  for(int i=0; i<Ntrig; i++){
    test_p[i] = new TH1D(Form("test_p%d",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    nume_p[i] = new TH1D(Form("nume_p%d",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    nume_e[i] = new TH1D(Form("nume_e%d",i),";#eta;Events",Nraparr-1,raparr);
    nume_phi[i] = new TH1D(Form("nume_phi%d",i),";#phi;Events",Nphiarr-1,phiarr);
    deno_p[i] = new TH1D(Form("deno_p%d",i),";p_{T}(GeV/c);Events",Nptarr-1,ptarr);
    deno_e[i] = new TH1D(Form("deno_e%d",i),";#eta;Events",Nraparr-1,raparr);
    deno_phi[i] = new TH1D(Form("deno_phi%d",i),";#phi;Events",Nphiarr-1,phiarr);
  }

  for(int i=0; i<nevent; i++){ // i<nevent
   fcha->GetEvent(i);
   if(i%500000==0){cout<<">>>>> EVENT "<<i<<" / "<<fcha->GetEntries()<<" ("<<(int)(100.*i/fcha->GetEntries())<<"%)"<<endl;}

   if(isTrk==true){
   //Track triggers (Done TnP style)
      for(int j=0; j<Reco_QQ_size; j++){
        TLorentzVector *recoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(j);
        TLorentzVector *recoQQpl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(j);
        TLorentzVector *recoQQmi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(j);

        //SoftMuon Cut for track trigger - Should they be different from soft muon cut for doubleMu trigs?
        Bool_t Condmi = true;
	Bool_t Condpl = true;
        Condmi = Condmi && (Reco_QQ_mumi_isGoodMuon[j]==1);
        Condmi = Condmi && (Reco_QQ_mumi_nTrkWMea[j] > 5);
        Condmi = Condmi && (Reco_QQ_mumi_nPixWMea[j] > 0);
        Condmi = Condmi && (fabs(Reco_QQ_mumi_dxy[j]) < 0.3);
        Condmi = Condmi && (fabs(Reco_QQ_mumi_dz[j]) < 20.);
        Condpl = Condpl && (Reco_QQ_mupl_isGoodMuon[j]==1);
        Condpl = Condpl && (Reco_QQ_mupl_nTrkWMea[j] > 5);
        Condpl = Condpl && (Reco_QQ_mupl_nPixWMea[j] > 0);
        Condpl = Condpl && (fabs(Reco_QQ_mupl_dxy[j]) < 0.3);
        Condpl = Condpl && (fabs(Reco_QQ_mupl_dz[j]) < 20.);

        //for Fill->Reco_QQ_singleMuon with SingleMu track trigger Matching
        for(int k=0; k<Ntrig; k++){
          if( //Cond&&//soft muon cut                
              Reco_QQ_sign[j]==0&&//opposite sign  
              recoQQ4mom->M()>massmin && recoQQ4mom->M()<massmax&&//mass window
              AcceptanceCut(recoQQpl4mom)&&AcceptanceCut(recoQQmi4mom)&&//Acceptance cut
	      (k>=26)  // only simgleMu track triggers
            ) // Track trigger conditions not including soft muon cuts
          {
	    if(k==26||k==28){L3SingleMuTrig = 18;} // L3Mu3
	    else if(k==27||k==29){L3SingleMuTrig = 20;} // L3Mu5
            if(
                ((Reco_QQ_mupl_trig[j]&((ULong64_t)pow(2, L3SingleMuTrig)))==((ULong64_t)pow(2, L3SingleMuTrig)))&& // Tag muon matched to L3 SingleMu3
		Condpl
                //recoQQpl4mom->Pt()>4
              )//select tag (mupl)
            {
              // Inside probe (mumi)
	      // Denominator: probe belongs to event which fired L3Mu3
	      if(
	         ((HLTriggers&((ULong64_t)pow(2, L3SingleMuTrig)))==((ULong64_t)pow(2, L3SingleMuTrig)))&& // Event fired L3 SingleMu3
		 Condmi
	        )
	      {
	         //Fill deno
                 deno_p[k]->Fill(recoQQmi4mom->Pt());
                 if((k==26||k==28)&&recoQQmi4mom->Pt()>3) deno_e[k]->Fill(recoQQmi4mom->Eta());
                 if((k==27||k==29)&&recoQQmi4mom->Pt()>5) deno_e[k]->Fill(recoQQmi4mom->Eta()); // for 5
	      } // if deno

	      if(
                  ((HLTriggers&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k)))&& // Event fired track trigger
		  Condmi
		)
              {
		 //Fill nume
                 nume_p[k]->Fill(recoQQmi4mom->Pt());
                 if((k==26||k==28)&&recoQQmi4mom->Pt()>3) nume_e[k]->Fill(recoQQmi4mom->Eta());
                 if((k==27||k==29)&&recoQQmi4mom->Pt()>5) nume_e[k]->Fill(recoQQmi4mom->Eta());
	      } //if nume
	    } // for tag mupl

            if(
                ((Reco_QQ_mumi_trig[j]&((ULong64_t)pow(2, L3SingleMuTrig)))==((ULong64_t)pow(2, L3SingleMuTrig)))&& // Tag muon matched to L3 SingleMu3
		Condmi
                //recoQQmi4mom->Pt()>4
              )//select tag (mumi)
            {
              // Inside probe (mupl)
              // Denominator: probe belongs to event which fired L3Mu3
              if(
                 ((HLTriggers&((ULong64_t)pow(2, L3SingleMuTrig)))==((ULong64_t)pow(2, L3SingleMuTrig)))&& // Event fired L3 SingleMu3
		 Condpl
                )
              {
                 //Fill deno
                 deno_p[k]->Fill(recoQQpl4mom->Pt());
                 if((k==26||k==28)&&recoQQpl4mom->Pt()>3) deno_e[k]->Fill(recoQQpl4mom->Eta());
                 if((k==27||k==29)&&recoQQpl4mom->Pt()>5) deno_e[k]->Fill(recoQQpl4mom->Eta()); // for 5
	      } //if deno

              if(
                  ((HLTriggers&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k)))&&
		  Condpl
		)
              {  
		 //Fill nume
                 nume_p[k]->Fill(recoQQpl4mom->Pt());
                 if((k==26||k==28)&&recoQQpl4mom->Pt()>3) nume_e[k]->Fill(recoQQpl4mom->Eta());
                 if((k==27||k==29)&&recoQQpl4mom->Pt()>5) nume_e[k]->Fill(recoQQpl4mom->Eta());
              } //if nume
            } //for tag mumi
          } // if Track trigger conditions
        } // For Ntrig - looping over triggers	
      } // For Filling with Reco_QQ_SingleMuTrkTrigger
   } // Track triggers (isTrk true)

   else{ // Not track
    if(usetnp==true){
      //for DoubleMu trigger - ALL
      for(int j=0; j<Reco_QQ_size; j++){
        TLorentzVector *recoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(j);
        TLorentzVector *recoQQpl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(j);
        TLorentzVector *recoQQmi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(j);
        //SoftMuon Cut for DoubleMu trigger
        Bool_t Cond = true;
        Cond = Cond && (Reco_QQ_mumi_isGoodMuon[j]==1);
        Cond = Cond && (Reco_QQ_mumi_nTrkWMea[j] > 5);
        Cond = Cond && (Reco_QQ_mumi_nPixWMea[j] > 0);
        Cond = Cond && (fabs(Reco_QQ_mumi_dxy[j]) < 0.3);
        Cond = Cond && (fabs(Reco_QQ_mumi_dz[j]) < 20.);
        Cond = Cond && (Reco_QQ_mupl_isGoodMuon[j]==1);
        Cond = Cond && (Reco_QQ_mupl_nTrkWMea[j] > 5);
        Cond = Cond && (Reco_QQ_mupl_nPixWMea[j] > 0);
        Cond = Cond && (fabs(Reco_QQ_mupl_dxy[j]) < 0.3);
        Cond = Cond && (fabs(Reco_QQ_mupl_dz[j]) < 20.);

        //for Fill->RecoQQ with DoubleMuon trigger Matching, except DoubleMu10
        for(int k=0; k<Ntrig; k++){
          if( Cond&&//soft muon cut                
              Reco_QQ_sign[j]==0&&//opposite sign  
              recoQQ4mom->M()>massmin && recoQQ4mom->M()<massmax&&//mass window
              AcceptanceCut(recoQQpl4mom)&&AcceptanceCut(recoQQmi4mom)&&//Acceptance cut
	      k<=9&&
              k!=5&&k!=7&&k!=9//except DoubleMu10 trigs
            ) // DoubleMu trigger conditions, except DoubleMu10 and SingleMu trigs
          {
            if(
                ((Reco_QQ_mupl_trig[j]&((ULong64_t)pow(2, minpTtrigfordeno)))==((ULong64_t)pow(2, minpTtrigfordeno)))&& 
                ((HLTriggers&((ULong64_t)pow(2, minpTtrigfordeno)))==((ULong64_t)pow(2, minpTtrigfordeno)))&& 
                recoQQpl4mom->Pt()>4
              )//select tag (mupl)
            {
              // Inside probe (mumi)
              deno_p[k]->Fill(recoQQmi4mom->Pt());
              //if(k<=9)
              //{
                deno_e[k]->Fill(recoQQmi4mom->Eta());//denominator for DoubleMu trigs
		//if(recoQQmi4mom->Pt()>10){deno_e[k]->Fill(recoQQmi4mom->Eta());} //Denominator for DoubleMu trigs
              //}

              if(
                  ((Reco_QQ_mumi_trig[j]&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k))) && 
                  ((HLTriggers&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k)))) 
              {
                nume_p[k]->Fill(recoQQmi4mom->Pt());
                //DoubleMu trig
	        //Santona
		//if(k<=9){
                  nume_e[k]->Fill(recoQQmi4mom->Eta());
		  //if(recoQQmi4mom->Pt()>10){nume_e[k]->Fill(recoQQmi4mom->Eta());}
                //}
              }  // Mumi probe numerator
            } // Mupl Tag 

            //Reco_QQ_mumi
            else if(
                ((Reco_QQ_mumi_trig[j]&((ULong64_t)pow(2, minpTtrigfordeno)))==((ULong64_t)pow(2, minpTtrigfordeno))) && 
                ((HLTriggers&((ULong64_t)pow(2, minpTtrigfordeno)))==((ULong64_t)pow(2, minpTtrigfordeno)))&&  
                recoQQmi4mom->Pt()>4
              ) //select tag (mumi)
            {
              deno_p[k]->Fill(recoQQpl4mom->Pt());
              //if(k<=9){
                deno_e[k]->Fill(recoQQpl4mom->Eta());
              //}

              if(
                  ((Reco_QQ_mupl_trig[j]&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k))) && 
                  ((HLTriggers&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k)))
                ) 
              {
                nume_p[k]->Fill(recoQQpl4mom->Pt());
                //DoubleMu trig
	        //if(k<=9){
                  nume_e[k]->Fill(recoQQpl4mom->Eta());
                //}
             } // mupl probe numerator
            }// mumi tag  
          } // If DoubleMu Trigger condition,  EXCEPT DOUBLEMU10
        }//for Fill - DoubleMu triggers except DoubleMu10

        //ONLY DoubleMu10 for enabling pT>10 cut
        for(int k=0; k<Ntrig; k++){
          if( Cond&&//soft muon cut                
              Reco_QQ_sign[j]==0&&//opposite sign  
              recoQQ4mom->M()>massmin && recoQQ4mom->M()<massmax&&
              AcceptanceCut(recoQQpl4mom)&&AcceptanceCut(recoQQmi4mom)&&//Acceptance cut
              (k==5 || k==7 || k==9)
            )
          {
            if(
                ((Reco_QQ_mupl_trig[j]&((ULong64_t)pow(2, minpTtrigfordeno)))==((ULong64_t)pow(2, minpTtrigfordeno))) && 
                ((HLTriggers&((ULong64_t)pow(2, minpTtrigfordeno)))==((ULong64_t)pow(2, minpTtrigfordeno)))&&  
                recoQQpl4mom->Pt()>4
              ) // Select tag muon (mupl) 
            {
              deno_p[k]->Fill(recoQQmi4mom->Pt());
	      if(recoQQmi4mom->Pt()>10){deno_e[k]->Fill(recoQQmi4mom->Eta());}  // Santona: only pT>10
              if(
                  ((Reco_QQ_mumi_trig[j]&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k))) && 
                  ((HLTriggers&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k)))) {
                nume_p[k]->Fill(recoQQmi4mom->Pt()); //,trigPrescale[k]);
                if(recoQQmi4mom->Pt()>10){nume_e[k]->Fill(recoQQmi4mom->Eta());}
              } //mumi probe numerator
            } // mupl tag 
            else if(
                ((Reco_QQ_mumi_trig[j]&((ULong64_t)pow(2, minpTtrigfordeno)))==((ULong64_t)pow(2, minpTtrigfordeno))) && 
                ((HLTriggers&((ULong64_t)pow(2, minpTtrigfordeno)))==((ULong64_t)pow(2, minpTtrigfordeno)))&& //L3Mu3
                recoQQmi4mom->Pt()>4
              ) // Select tag muom (mumi)
            {
              deno_p[k]->Fill(recoQQpl4mom->Pt());
	      if(recoQQpl4mom->Pt()>10){deno_e[k]->Fill(recoQQpl4mom->Eta());}
              if(
                  ((Reco_QQ_mupl_trig[j]&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k))) && 
                  ((HLTriggers&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k)))
                ) 
              {
                nume_p[k]->Fill(recoQQpl4mom->Pt()); //,trigPrescale[k]);
                if(recoQQpl4mom->Pt()>10){nume_e[k]->Fill(recoQQpl4mom->Eta());}
              }//mupl probe numerator
            }//mumi tag
          }//if DoubleMu10 trigger conditions (need to fix mass - should be same mass range as other DoubleMu triggers)
        }//for Fill for DoubleMu10 ONLY
      }//for DoubleMu triggers - ALL
    }//for usetnp True

    //Not using  TnP start:
    // Santona: USE THIS FOR SIMGLEMU TRIGGERS
    else{
      if(i==10||i==100||i==1000){cout<<"@@@@@@@@@@@@@@@@@@@@@@ IT'S NOT! TNP @@@@@@@@@@@@@@@@@@@@@"<<endl;}
      //for SingleMu trigger. Uses Reco_mu
      for(int j=0; j<Reco_mu_size; j++){
        TLorentzVector *recomu4mom = (TLorentzVector*) Reco_mu_4mom->At(j);
        //SoftMuon Cut
        Bool_t SingleCond = true;
        SingleCond = SingleCond && (Reco_mu_isGoodMuon[j]==1);
        SingleCond = SingleCond && (Reco_mu_nTrkWMea[j] > 5);
        SingleCond = SingleCond && (Reco_mu_nPixWMea[j] > 0);
        SingleCond = SingleCond && (fabs(Reco_mu_dxy[j]) < 0.3);
        SingleCond = SingleCond && (fabs(Reco_mu_dz[j]) < 20.);
        //for Fill->Reco_mu_4mom
        for(int k=0; k<Ntrig; k++){
          //SingleMu trigger selection // no trigger matching requirement is placed for denominator since it's not tag and probe style
	  if (k>=10 &&
              AcceptanceCut(recomu4mom) &&
              SingleCond
             )
          {
            deno_p[k]->Fill(recomu4mom->Pt());
            if((k==12||k==18||k==19)&&recomu4mom->Pt()>3) deno_e[k]->Fill(recomu4mom->Eta());
            if((k==13||k==20||k==21)&&recomu4mom->Pt()>5) deno_e[k]->Fill(recomu4mom->Eta());
            if((k==14||k==22)&&recomu4mom->Pt()>7) deno_e[k]->Fill(recomu4mom->Eta());
            if((k==10||k==15||k==23)&&recomu4mom->Pt()>12) deno_e[k]->Fill(recomu4mom->Eta());
            if((k==16||k==24)&&recomu4mom->Pt()>15) deno_e[k]->Fill(recomu4mom->Eta());
            if(k==11&&recomu4mom->Pt()>16) deno_e[k]->Fill(recomu4mom->Eta());
            if((k==17||k==25)&&recomu4mom->Pt()>20) deno_e[k]->Fill(recomu4mom->Eta());

            if(
                ((Reco_mu_trig[j]&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k))) &&
                ((HLTriggers&((ULong64_t)pow(2, k)))==((ULong64_t)pow(2, k)))
              )
            {
              nume_p[k]->Fill(recomu4mom->Pt());
              if((k==12||k==18||k==19)&&recomu4mom->Pt()>3) nume_e[k]->Fill(recomu4mom->Eta());
              if((k==13||k==20||k==21)&&recomu4mom->Pt()>5) nume_e[k]->Fill(recomu4mom->Eta());
              if((k==14||k==22)&&recomu4mom->Pt()>7) nume_e[k]->Fill(recomu4mom->Eta());
              if((k==10||k==15||k==23)&&recomu4mom->Pt()>12) nume_e[k]->Fill(recomu4mom->Eta());
              if((k==16||k==24)&&recomu4mom->Pt()>15) nume_e[k]->Fill(recomu4mom->Eta());
              if(k==11&&recomu4mom->Pt()>16) nume_e[k]->Fill(recomu4mom->Eta());
              if((k==17||k==25)&&recomu4mom->Pt()>20) nume_e[k]->Fill(recomu4mom->Eta());
            }//if nume (trigger matching)
          }//if SingleMu trigger selection (All muons, no trigger matching required)
        }//for Fill - SingleMu triggers
      }//for SingleMu triggers
    }//NOT use TnP
   }//Not Track
  }//for nevent



  gStyle->SetOptStat(1);
  gStyle->SetFillColor(0);
  gSystem->mkdir(Form("figs/%s/%s/%s",date.c_str(),dataset.c_str(),vername.c_str()),1);

  TLine *lptx3 = new TLine(3.,0, 3.,1);
  TLine *lptx5 = new TLine(5.,0, 5.,1);
  TLine *lptx7 = new TLine(7.,0, 7.,1);
  TLine *lptx10 = new TLine(10.,0, 10.,1);
  TLine *lptx12 = new TLine(12.,0, 12.,1);
  TLine *lptx15 = new TLine(15.,0, 15.,1);
  TLine *lptx16 = new TLine(16.,0, 16.,1);
  TLine *lptx20 = new TLine(20.,0, 20.,1);
  TLine *lpty[Ntrig];
  TLine *letay[Ntrig];
  lptx3->SetLineStyle(2);   lptx3->SetLineWidth(2);  lptx3->SetLineColor(kMagenta+2);
  lptx5->SetLineStyle(2);   lptx5->SetLineWidth(2);  lptx5->SetLineColor(kRed+2);
  lptx7->SetLineStyle(2);   lptx7->SetLineWidth(2);  lptx7->SetLineColor(kGreen+2);
  lptx10->SetLineStyle(2); lptx10->SetLineWidth(2); lptx10->SetLineColor(kRed+2);
  lptx12->SetLineStyle(2); lptx12->SetLineWidth(2); lptx12->SetLineColor(kBlue+2);
  lptx15->SetLineStyle(2); lptx15->SetLineWidth(2); lptx15->SetLineColor(kMagenta+2);
  lptx16->SetLineStyle(2); lptx16->SetLineWidth(2); lptx16->SetLineColor(kBlue+2);
  lptx20->SetLineStyle(2); lptx20->SetLineWidth(2); lptx20->SetLineColor(kRed+2);

  TCanvas *c1;
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  TCanvas *c3 = new TCanvas("c3","c3",1200,600);
  c3->Divide(2,1);
  TCanvas *c4 = new TCanvas("c4","c4",1200,600);
  c4->Divide(2,1);
  TCanvas *c5 = new TCanvas("c5","c5",1200,600);
  c5->Divide(2,1);
  TCanvas *c6 = new TCanvas("c6","c6",1200,600);
  c6->Divide(2,1);
  TCanvas *c7 = new TCanvas("c7","c7",1200,600);
  c7->Divide(2,1);
  TCanvas *c8 = new TCanvas("c8","c8",1200,600);
  c8->Divide(2,1);
  TCanvas *c9 = new TCanvas("c9","c9",1200,600);
  c9->Divide(2,1);
  TCanvas *c10 = new TCanvas("c10","c10",1200,600);
  c10->Divide(2,1);
  TCanvas *c11 = new TCanvas("c11","c11",1200,600);
  c11->Divide(2,1);
  TCanvas *test = new TCanvas("test","test",1200,600);
  test->Divide(2,1);

  TLegend *leg[Ntrig];
  TLegend *leg1[Ntrig];
  TLegend *leg2[Ntrig];
  TLegend *leg3[Ntrig];
  TLegend *leg4[Ntrig];
  TLegend *leg5[Ntrig];
  TLegend *leg6[Ntrig];
  TLegend *leg7[Ntrig];
  TLegend *leg8[Ntrig];
  TLegend *leg9[Ntrig];
  TLegend *leg10[Ntrig];

  //Set histograms & Draw histograms
  for(int i=0; i<Ntrig; i++){
    cout<<"Enter histograms : "<<i<<endl;
    eff_p[i] = new TGraphAsymmErrors(nume_p[i],deno_p[i]);
    eff_e[i] = new TGraphAsymmErrors(nume_e[i],deno_e[i]);
    eff_p[i]->GetXaxis()->SetLimits(ptmin,ptmax);
    eff_e[i]->GetXaxis()->SetLimits(-2.4,2.4);

    SetHistStyle(eff_p[i],i,i,0,1.7);
    SetHistStyle(eff_e[i],i,i,0,1.7);
    cout<<"*************************"<<endl;
    cout<<"*         eff_p:        *"<<eff_p[i]<<endl;
    cout<<"*        trigname:      *"<<trigname[i].c_str()<<endl;
    cout<<"*       numerator:      *"<<nume_p[i]->GetEntries()<<endl;
    cout<<"*      denominator:     *"<<deno_p[i]->GetEntries()<<endl;
    cout<<"*************************"<<endl;
    for(int a=0; a<eff_e[i]->GetN(); a++){
      double x;
      double y;
      eff_e[i]->GetPoint(a,x,y);
      //cout<<"\nGetPoint eta distribution.\n"<<endl;
      cout<<"\t"<<"eta bin number: "<<a<<" x: "<<x<<" y: "<<y<<endl;
    }
    leg[i] = new TLegend(0.175,0.850,0.55,0.875);
    leg[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
    SetLegendStyle(leg[i]);
    leg[i]->SetFillStyle(0);
    lpty[i] = new TLine(ptmin,1, ptmax,1);
    lpty[i]->SetLineStyle(2);
    lpty[i]->SetLineWidth(2);
    letay[i] = new TLine(etamin,1, etamax,1);
    letay[i]->SetLineStyle(2);
    letay[i]->SetLineWidth(2);

    // Use TnP - for plotting
    if(usetnp==true){
      c1 = new TCanvas("c1","c1",1200,600);
      c1->Divide(2,1);
      c1->cd(1);
      eff_p[i]->Draw("ap");
      leg[i]->Draw("sames");//cout<<"legend name: "<<Form("%s",trigname[i].c_str())<<endl;
      eff_p[i]->SaveAs(
          Form("figs/%s/%s/%s/h_pt_efficiency_%s.root",date.c_str(),dataset.c_str(),vername.c_str(),trigname[i].c_str())
          );
      //Set line for pt TurnOn
      if(i==5 || i==7 || i==9) {lptx10->Draw("sames");} 
      lpty[i]->Draw("sames");
      c1->cd(2);
      eff_e[i]->Draw("ap");
      leg[i]->Draw("sames");
      eff_p[i]->GetHistogram()->GetXaxis()->SetTitle("p_{T}(GeV/c)");
      eff_e[i]->GetHistogram()->GetXaxis()->SetTitle("#eta");
      eff_p[i]->GetHistogram()->GetYaxis()->SetTitle("Efficiency");
      eff_e[i]->GetHistogram()->GetYaxis()->SetTitle("Efficiency");
      eff_e[i]->SaveAs(
          Form("figs/%s/%s/%s/h_eta_efficiency_%s.root",date.c_str(),dataset.c_str(),vername.c_str(),trigname[i].c_str())
          );
      letay[i]->Draw("sames");
      c1->SaveAs(
          Form("figs/%s/%s/%s/efficiency_%s.png",date.c_str(),dataset.c_str(),vername.c_str(),trigname[i].c_str())
          );

      gStyle->SetOptStat(1);
      test->cd(1); deno_phi[12]->Draw("pe"); test->cd(2); nume_phi[12]->Draw("pe");
      test->SaveAs("test.png");

      //for turnOn all double mu trig
      if(i<=9){
        cout<<"DoubleMu triggers : "<<i<<endl;
        c2->cd(1); 
        if(i==0) eff_p[i]->Draw("alp"); 
        eff_p[i]->Draw("lp"); 
        lpty[i]->Draw("sames");
        leg1[0] = new TLegend(0.15,0.850,0.50,0.875);
        leg1[1] = new TLegend(0.15,0.825,0.50,0.850);
        leg1[2] = new TLegend(0.15,0.800,0.50,0.825);
        leg1[3] = new TLegend(0.15,0.775,0.50,0.800);
        leg1[4] = new TLegend(0.15,0.750,0.50,0.775);

        leg1[5] = new TLegend(0.55,0.850,0.90,0.875);
        leg1[6] = new TLegend(0.55,0.825,0.90,0.850);
        leg1[7] = new TLegend(0.55,0.800,0.90,0.825);
        leg1[8] = new TLegend(0.55,0.775,0.90,0.800);
        leg1[9] = new TLegend(0.55,0.750,0.90,0.775);
        leg1[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
        SetLegendStyle(leg1[i]);
        leg1[i]->SetFillStyle(0);
        leg1[i]->Draw("sames");
        c2->cd(2); 
        if(i==0) eff_e[0]->Draw("ap"); 
        eff_e[i]->Draw("p");
        letay[i]->Draw("sames");
        leg1[i]->Draw("sames");
        c2->SaveAs(
            Form("figs/%s/%s/%s/efficiency_DoubleMuTurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
            );
      }
      // Changing to just OS and SS (double mu0) (check)
      if(i==1||i==2){  
        c3->cd(1); 
        if(i==1) eff_p[i]->Draw("alp");
        eff_p[i]->Draw("lp"); 
        lpty[i]->Draw("sames");
        leg2[1] = new TLegend(0.175,0.850,0.70,0.875);
        leg2[2] = new TLegend(0.175,0.825,0.70,0.850);
        leg2[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
        SetLegendStyle(leg2[i]);
        leg2[i]->SetFillStyle(0);
        leg2[i]->Draw("sames");
        c3->cd(2); 
        if(i==1) eff_e[i]->Draw("ap"); 
        eff_e[i]->Draw("p"); 
        letay[i]->Draw("sames");
        leg2[i]->Draw("sames");
        c3->SaveAs(
            Form("figs/%s/%s/%s/efficiency_DoubleL1Mu_Open_OS_SS_TurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
            );
      }
      //for turnOn all L1  double mu Open_0_HighQ trig (check)
      if(i==0||i==3||i==4){ 
        c4->cd(1); 
        if(i==0) eff_p[i]->Draw("alp"); 
        eff_p[i]->Draw("lp"); 
        lpty[i]->Draw("sames");
        leg3[0] = new TLegend(0.175,0.850,0.70,0.875);
        leg3[3] = new TLegend(0.175,0.825,0.70,0.850);
        leg3[4] = new TLegend(0.175,0.800,0.70,0.825);
        leg3[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
        SetLegendStyle(leg3[i]);
        leg3[i]->SetFillStyle(0);
        leg3[i]->Draw("sames");
        c4->cd(2); 
        if(i==0) eff_e[0]->Draw("ap"); 
        eff_e[i]->Draw("p");
        letay[i]->Draw("sames");
        leg3[i]->Draw("sames");
        c4->SaveAs(
            Form("figs/%s/%s/%s/efficiency_DoubleL1Mu_Open_0_HighQ_TurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
            );
      }
      //for turnOn all double mu L1L2L3_0 (check)
      if(i==3||i==6||i==8){  
        c5->cd(1); 
        if(i==3) eff_p[i]->Draw("alp"); 
        eff_p[i]->Draw("lp"); 
        lpty[i]->Draw("sames");
        leg4[3] = new TLegend(0.175,0.850,0.70,0.875);
        leg4[6] = new TLegend(0.175,0.825,0.70,0.850);
        leg4[8] = new TLegend(0.175,0.800,0.70,0.825);
        leg4[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
        SetLegendStyle(leg4[i]);
        leg4[i]->SetFillStyle(0);
        leg4[i]->Draw("sames");
        c5->cd(2); 
        if(i==3) eff_e[i]->Draw("ap"); 
        eff_e[i]->Draw("p");
        letay[i]->Draw("sames");
        leg4[i]->Draw("sames");
        c5->SaveAs(
            Form("figs/%s/%s/%s/efficiency_DoubleL1L2L3Mu_0_TurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
            );
      }
      //for turnOn all double mu L1L2L3_10 (check - although only needed for Z)
      if(i==5||i==7||i==9){ 
        c6->cd(1);
        if(i==5) eff_p[i]->Draw("alp"); 
        eff_p[i]->Draw("lp"); 
        lptx10->Draw("sames");
        lpty[i]->Draw("sames");
        leg5[5] = new TLegend(0.175,0.850,0.70,0.875);
        leg5[7] = new TLegend(0.175,0.825,0.70,0.850);
        leg5[9] = new TLegend(0.175,0.800,0.70,0.825);
        leg5[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
        SetLegendStyle(leg5[i]);
        leg5[i]->SetFillStyle(0);
        leg5[i]->Draw("sames");
        c6->cd(2); 
        if(i==5) eff_e[i]->Draw("ap"); 
        eff_e[i]->Draw("p");
        letay[i]->Draw("sames");
        leg5[i]->Draw("sames");
        c6->SaveAs(
            Form("figs/%s/%s/%s/efficiency_DoubleL1L2L3Mu_10_TurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
            );
      }

/*      //for turnOn single mu trig L1Mu12,L2Mu12,L3Mu12 (check)
      if(i==13||i==18||i==24){
        c7->cd(1); 
        if(i==13) eff_p[i]->Draw("alp"); 
        eff_p[i]->Draw("lp"); 
        lptx12->Draw("sames");
        lpty[i]->Draw("sames");
        leg6[13] = new TLegend(0.175,0.850,0.70,0.875);
        leg6[18] = new TLegend(0.175,0.825,0.70,0.850);
	leg6[24] = new TLegend(0.175,0.800,0.70,0.825);
        leg6[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
        SetLegendStyle(leg6[i]);
        leg6[i]->SetFillStyle(0);
        leg6[i]->Draw("sames");
        c7->cd(2); 
        if(i==13) eff_e[i]->Draw("ap"); 
        eff_e[i]->Draw("p");
        letay[i]->Draw("sames");
        leg6[i]->Draw("sames");
        c7->SaveAs(
            Form("figs/%s/%s/%s/efficiency_SingleL1L2L3Mu12TurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
            );
      }
      //for turnOn single mu trig L2Mu15,L3Mu15
      if(i==19||i==25){ 
        c8->cd(1); 
        if(i==19) eff_p[i]->Draw("alp"); 
        eff_p[i]->Draw("lp"); 
        lptx15->Draw("sames");
        lpty[i]->Draw("sames");
        leg7[19] = new TLegend(0.175,0.850,0.70,0.875);
        leg7[25] = new TLegend(0.175,0.825,0.70,0.850);
        leg7[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
        SetLegendStyle(leg7[i]);
        leg7[i]->SetFillStyle(0);
        leg7[i]->Draw("sames");
        c8->cd(2); 
        if(i==19) eff_e[i]->Draw("ap"); 
        if(i==25) eff_e[i]->Draw("p");
        letay[i]->Draw("sames");
        leg7[i]->Draw("sames");
        c8->SaveAs(
            Form("figs/%s/%s/%s/efficiency_SingleL2L3Mu15TurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
            );
      }
      //for turnOn single mu trig L1Mu16
      if(i==14){
        c9->cd(1); 
        if(i==14) eff_p[i]->Draw("alp"); 
        eff_p[i]->Draw("lp"); 
        lptx16->Draw("sames");
        lpty[i]->Draw("sames");
        leg8[14] = new TLegend(0.175,0.850,0.70,0.875);
        leg8[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
        SetLegendStyle(leg8[i]);
        leg8[i]->SetFillStyle(0);
        leg8[i]->Draw("sames");
        c9->cd(2); 
        if(i==14) eff_e[i]->Draw("ap"); 
        letay[i]->Draw("sames");
        leg8[i]->Draw("sames");
        c9->SaveAs(
            Form("figs/%s/%s/%s/efficiency_SingleL1Mu16TurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
            );   
      }    
      //for turnOn single mu trig L2Mu20,L3Mu20
      if(i==20||i==26){
        c10->cd(1); 
        if(i==20) eff_p[i]->Draw("alp");
        eff_p[i]->Draw("lp");
        lptx20->Draw("sames");
        lpty[i]->Draw("sames");
        leg9[20] = new TLegend(0.175,0.850,0.70,0.875); 
        leg9[26] = new TLegend(0.175,0.825,0.70,0.850);
        leg9[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
        SetLegendStyle(leg9[i]);
        leg9[i]->SetFillStyle(0);
        leg9[i]->Draw("sames");
        c10->cd(2); 
        if(i==20) eff_e[i]->Draw("ap"); 
        if(i==26) eff_e[i]->Draw("p");
        letay[i]->Draw("sames");
        leg9[i]->Draw("sames");
        c10->SaveAs(
            Form("figs/%s/%s/%s/efficiency_SingleL2L3Mu20TurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
            );   
      }    

      //for turnOn single mu trig 3,5,7,12,15,16,20 //I'm including all L1,L2,L3
      if(i==12||i==11||i==12||i==13||i==14||i==15||i==16||i==17||i==18||i==19||i==20||i==21||i==22||i==23||i==24||i==25||i==26){ //Santona //was 14,15,16,17,18 //This one is confusing. Why was L2 left out?
        c11->cd(1); 
        if(i==10) eff_p[i]->Draw("alp"); 
        eff_p[i]->Draw("lp"); 
        lptx3->Draw("sames"); lptx5->Draw("sames"); lptx7->Draw("sames"); lptx12->Draw("sames"); lptx15->Draw("sames"); lptx16->Draw("sames"); lptx20->Draw("sames"); 
        lpty[i]->Draw("sames");
        leg10[10] = new TLegend(0.15,0.850,0.50,0.875);
        leg10[11] = new TLegend(0.15,0.825,0.50,0.850);
        leg10[12] = new TLegend(0.15,0.800,0.50,0.825);
        leg10[13] = new TLegend(0.15,0.775,0.50,0.800);
        leg10[14] = new TLegend(0.15,0.750,0.50,0.775);
        leg10[15] = new TLegend(0.15,0.725,0.50,0.750);
        leg10[16] = new TLegend(0.15,0.700,0.50,0.725);
        leg10[17] = new TLegend(0.15,0.675,0.50,0.700);
        leg10[18] = new TLegend(0.15,0.650,0.50,0.675);

        leg10[19] = new TLegend(0.55,0.850,0.90,0.875);
        leg10[20] = new TLegend(0.55,0.825,0.90,0.850);
        leg10[21] = new TLegend(0.55,0.800,0.90,0.825);
        leg10[22] = new TLegend(0.55,0.775,0.90,0.800);
        leg10[23] = new TLegend(0.55,0.750,0.90,0.775);
        leg10[24] = new TLegend(0.55,0.725,0.90,0.750);
        leg10[25] = new TLegend(0.55,0.700,0.90,0.725);
        leg10[26] = new TLegend(0.55,0.675,0.90,0.700);
        leg10[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
        SetLegendStyle(leg10[i]);
        leg10[i]->SetFillStyle(0);
        leg10[i]->Draw("sames");
        c11->cd(2); 
        if(i==10) eff_e[i]->Draw("ap"); 
        eff_e[i]->Draw("p");
        letay[i]->Draw("sames");
        leg10[i]->Draw("sames");
        c11->SaveAs(
            Form("figs/%s/%s/%s/efficiency_SingleMuTurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
            );
      }
// */

      if(isTrk==true){
      //for L3 singleMu track triggers 3,5 (check)
        if(i==26||i==27||i==28||i==29){
          c8->cd(1);
          if(i==26) eff_p[i]->Draw("alp");
          eff_p[i]->Draw("lp");
          //lptx3->Draw("sames"); lptx5->Draw("sames");
          lpty[i]->Draw("sames");
          leg9[26] = new TLegend(0.175,0.850,0.70,0.875);
          leg9[27] = new TLegend(0.175,0.825,0.70,0.850);
          leg9[28] = new TLegend(0.175,0.800,0.70,0.825);
          leg9[29] = new TLegend(0.175,0.775,0.70,0.800);
          leg9[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
          SetLegendStyle(leg9[i]);
          leg9[i]->SetFillStyle(0);
          leg9[i]->Draw("sames");
          c8->cd(2);
          if(i==26) eff_e[i]->Draw("ap");
          eff_e[i]->Draw("p");
          letay[i]->Draw("sames");
          leg9[i]->Draw("sames");
          c8->SaveAs(
                Form("figs/%s/%s/%s/efficiency_SingleL3MuTrack.png",date.c_str(),dataset.c_str(),vername.c_str())
              );
        }
      } //Track triggers plot
    }//use TnP, for DoubleMu triggers and singleMu track triggers - For plotting

    // Not using TnP - for SingleMu trigs
    else{
      if(i>=10){
        c1 = new TCanvas("c1","c1",1200,600);
        c1->Divide(2,1);
        c1->cd(1);
        if(i==10) eff_p[i]->Draw("alp"); 
        eff_p[i]->Draw("lp");
        leg[i]->Draw("sames");//cout<<"legend name: "<<Form("%s",trigname[i].c_str())<<endl;
        eff_p[i]->SaveAs(
            Form("figs/%s/%s/%s/h_pt_efficiency_%s.root",date.c_str(),dataset.c_str(),vername.c_str(),trigname[i].c_str())
            );
        //Set line for pt TurnOn
        if(i==12||i==18||i==19) {lptx3->Draw("sames");} // track triggers 26 and 28
        if(i==13||i==20||i==21) {lptx5->Draw("sames");}  // Track triggers 27 and 29
        if(i==14||i==22) {lptx7->Draw("sames");}
        if(i==10 || i==15 || i==23) {lptx12->Draw("sames");}
        if(i==16 || i==24) {lptx15->Draw("sames");}
        if(i==11) {lptx16->Draw("sames");}
        if(i==17 || i==25) {lptx20->Draw("sames");}
        lpty[i]->Draw("sames");
        c1->cd(2);
        if(i==10) eff_e[i]->Draw("ap"); 
        eff_e[i]->Draw("p");
        leg[i]->Draw("sames");
        eff_p[i]->GetHistogram()->GetXaxis()->SetTitle("p_{T}(GeV/c)");
        eff_e[i]->GetHistogram()->GetXaxis()->SetTitle("#eta");
        eff_p[i]->GetHistogram()->GetYaxis()->SetTitle("Efficiency");
        eff_e[i]->GetHistogram()->GetYaxis()->SetTitle("Efficiency");
        eff_e[i]->SaveAs(
            Form("figs/%s/%s/%s/h_eta_efficiency_%s.root",date.c_str(),dataset.c_str(),vername.c_str(),trigname[i].c_str())
            );
        letay[i]->Draw("sames");
        c1->SaveAs(
            Form("figs/%s/%s/%s/efficiency_%s.png",date.c_str(),dataset.c_str(),vername.c_str(),trigname[i].c_str())
            );

        //for turnOn single mu trig L1Mu12,16 //was 3,5,7,12,16 (ALL) (check)
        if(i==10||i==11){
          c4->cd(1);
          if(i==10) eff_p[i]->Draw("alp");
          eff_p[i]->Draw("lp");
          lptx12->Draw("sames"); lptx16->Draw("sames");
          lpty[i]->Draw("sames");
          leg5[10] = new TLegend(0.175,0.850,0.70,0.875);
          leg5[11] = new TLegend(0.175,0.825,0.70,0.850);
          leg5[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
          SetLegendStyle(leg5[i]);
          leg5[i]->SetFillStyle(0);
          leg5[i]->Draw("sames");
          c4->cd(2);
          if(i==10) eff_e[i]->Draw("ap");
          eff_e[i]->Draw("p");
          letay[i]->Draw("sames");
          leg5[i]->Draw("sames");
          c4->SaveAs(
              Form("figs/%s/%s/%s/efficiency_SingleL1MuALLTurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
              );
        }

        //for turnOn single mu trig L2Mu3,5,7,12,15,20 (ALL) (check)
        if(i==12||i==13||i==14||i==15||i==16||i==17){ 
          c5->cd(1);
          if(i==12) eff_p[i]->Draw("alp");
          eff_p[i]->Draw("lp");
          lptx3->Draw("sames"); lptx5->Draw("sames"); lptx7->Draw("sames"); lptx12->Draw("sames"); lptx15->Draw("sames"); lptx20->Draw("sames");
          lpty[i]->Draw("sames");
          leg6[12] = new TLegend(0.175,0.850,0.70,0.875);
          leg6[13] = new TLegend(0.175,0.825,0.70,0.850);
          leg6[14] = new TLegend(0.175,0.800,0.70,0.825);
          leg6[15] = new TLegend(0.175,0.775,0.70,0.800);
          leg6[16] = new TLegend(0.175,0.750,0.70,0.775);
          leg6[17] = new TLegend(0.175,0.725,0.70,0.750);
          leg6[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
          SetLegendStyle(leg6[i]);
          leg6[i]->SetFillStyle(0);
          leg6[i]->Draw("sames");
          c5->cd(2);
          if(i==12) eff_e[i]->Draw("ap");
          eff_e[i]->Draw("p");
          letay[i]->Draw("sames");
          leg6[i]->Draw("sames");
          c5->SaveAs(
              Form("figs/%s/%s/%s/efficiency_SingleL2MuALLTurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
              );
        }

        //for turnOn single mu trig L3Mu3,5,7,12,15,20 (ALL except track) (check)
        if(i==18||i==19||i==20||i==21||i==22||i==23||i==24||i==25){ 
          c2->cd(1);
          if(i==18) eff_p[i]->Draw("alp");
          eff_p[i]->Draw("lp");
          lptx3->Draw("sames"); lptx5->Draw("sames"); lptx7->Draw("sames"); lptx12->Draw("sames"); lptx15->Draw("sames"); lptx20->Draw("sames");
          lpty[i]->Draw("sames");
          leg3[18] = new TLegend(0.175,0.850,0.70,0.875);
          leg3[19] = new TLegend(0.175,0.825,0.70,0.850);
          leg3[20] = new TLegend(0.175,0.800,0.70,0.825);
          leg3[21] = new TLegend(0.175,0.775,0.70,0.800);
          leg3[22] = new TLegend(0.175,0.750,0.70,0.775);
          leg3[23] = new TLegend(0.175,0.725,0.70,0.750);
          leg3[24] = new TLegend(0.175,0.700,0.70,0.725);
          leg3[25] = new TLegend(0.175,0.675,0.70,0.700);
          leg3[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
          SetLegendStyle(leg3[i]);
          leg3[i]->SetFillStyle(0);
          leg3[i]->Draw("sames");
          c2->cd(2);
          if(i==18) eff_e[i]->Draw("ap");
          eff_e[i]->Draw("p");
          letay[i]->Draw("sames");
          leg3[i]->Draw("sames");
          c2->SaveAs(
              Form("figs/%s/%s/%s/efficiency_SingleL3MuALLTurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
              );
        }

        //for turnOn single mu trig L1Mu12,L2Mu12,L3Mu12 (check)
        if(i==10||i==15||i==23){ 
          c7->cd(1); 
          if(i==10) eff_p[i]->Draw("alp"); 
          eff_p[i]->Draw("lp"); 
          lptx12->Draw("sames");
          lpty[i]->Draw("sames");
          leg6[10] = new TLegend(0.175,0.850,0.70,0.875);
          leg6[15] = new TLegend(0.175,0.825,0.70,0.850);
	  leg6[23] = new TLegend(0.175,0.800,0.70,0.825);
          leg6[i]->AddEntry(eff_p[i],Form("%s",trigname[i].c_str()),"lp");
          SetLegendStyle(leg6[i]);
          leg6[i]->SetFillStyle(0);
          leg6[i]->Draw("sames");
          c7->cd(2); 
          if(i==10) eff_e[i]->Draw("ap"); 
          eff_e[i]->Draw("p");
          letay[i]->Draw("sames");
          leg6[i]->Draw("sames");
          c7->SaveAs(
              Form("figs/%s/%s/%s/efficiency_SingleL1L2L3Mu12TurnOn.png",date.c_str(),dataset.c_str(),vername.c_str())
              );
        }

      } // If SingleMu trig - Plotting
    }//it does not use TnP - Plotting
  }//for Ntrig - Plotting
}//main END
