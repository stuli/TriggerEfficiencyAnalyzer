#include <TGraphAsymmErrors.h>
#include <iostream>
#include <string>
#include "JpsiFunc.h"

void getTriggers(){

  //TChain *fcha = new TChain("hionia");

  //fcha->Add("/eos/cms/store/group/phys_heavyions/jaebeom/OniaTree_PPref_SingleMuTnP/SingleMuonTnP/OniaTree_PPref_SingleMuTnP/171116_031052/0000/OniaForest_321.root");

  TFile* file = new TFile("/eos/cms/store/group/phys_heavyions/jaebeom/OniaTree_PPref_SingleMuTnP/SingleMuonTnP/OniaTree_PPref_SingleMuTnP/171116_031052/0000/OniaForest_321.root");

  file.ls();

  //TH1D* hmyHisto;
  //hTriggers = 
  //file->GetObject("hStats", hmyHisto);
  //hmyHisto = file->Get("hStats");
  //for (int i=2; i<hStats->GetNbinsX()/2+1; i++) cout << i-2 << " " << hStats->GetXaxis()->GetBinLabel(i) << endl;

  TH1F * hStats = (TH1F*)file.Get("hionia/hStats");

  //hmyHisto = file->Get("hionia/hStats");
  //hStats->Draw();
  //hmyHisto->Draw();
  //cout << hmyHisto->GetNbinsX() << endl;

  
  for (int i=2; i<hStats->GetNbinsX()/2+1; i++) cout << i-2 << " " << hStats->GetXaxis()->GetBinLabel(i) << endl;


}
