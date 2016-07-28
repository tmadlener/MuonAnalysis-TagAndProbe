//C/C++
#include <iostream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>
//Root
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TIterator.h"
#include "TObject.h"
#include "TKey.h"
#include "TText.h"
#include "TPaveText.h"
#include "TString.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TGaxis.h"

void ratio(){

  //steering variable: 1 = eta, 2 = abseta, 3 = pt, 4 = vertex
  // id = Tight, Loose, Soft
  // separation = "", _seagulls, _cowboys
  int scenario = 3;
  std::string id = "Loose";
  std::string title = "Loose ID";
  std::string separation = "";
  std::string absetaBin;
  if(scenario == 3) absetaBin = "1";
  std::string add = "_2012_highestBinsMerged"; //_2012_highestBinsMerged

  gROOT->SetStyle("Plain");
  //gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetPadLeftMargin(2.);
  //gStyle->SetPadBottomMargin(1.3);
  //gStyle->SetTitleYOffset(1.1);

  double x1 = 0, x2 = 0, y1 = 0, y2 = 0,
    t1 = 0.43, t2 = 0.8, t3 = 0.63, t4 = 0.9,
    lx = 0.65, hy2 = 1.05, hy1 = 0.;
  std::string xtitle, values;
  std::string scen;
  std::vector <double> bins;

  // binning
  //eta
  const std::vector<double> bins1 = {-2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1};

  //abseta
  const std::vector<double> bins2 = {0., 0.9, 1.2, 2.1};

  //pt
  const std::vector<double> bins31 = {2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 20.0};
  const std::vector<double> bins3 = {2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0};

  //vertex
  const std::vector<double> bins4 = {0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 22.5, 24.5, 26.5, 28.5, 30.5};


  switch(scenario){
  case 1:
    std::cout << "processing eta" << std::endl;
    xtitle = "#eta";
    scen = "eta";
    values = "p_{T} > 8 GeV/c";
    bins = bins1;
    x1 = -2.2; // -2.2
    x2 = 2.2; // 2.2
    y1 = 0.9; // 0.5
    y2 = 1.1; // 1.3
    hy1 = 0.75;
    hy2 = 1.1;
    break;
  case 2:
    std::cout << "processing abseta" << std::endl;
    xtitle = "|#eta|";
    scen = "plateau_abseta";
    values = "p_{T} > 8 GeV/c";
    bins = bins2;
    x1 = 0;
    x2 = 2.2;
    y1 = 0.3;
    y2 = 1.3;
    break;
  case 3:
    std::cout << "processing pt" << std::endl;
    xtitle = "p_{T} [GeV/c]";
    scen = "pt_abseta";
    if(absetaBin == "0"){
      values = "|#eta| < 0.9";
      bins = bins31;
      y2 = 1.2;
      if(id == "Tight") y2 = 1.4;
    }else if(absetaBin == "1"){
      values = "0.9 < |#eta| < 1.2";
      bins = bins3;
      y2 = 1.7;
    }else if(absetaBin == "2"){
      values = "1.2 < |#eta| < 2.1";
      bins = bins3;
      y2 = 1.2;
    }
    x1 = 1;
    x2 = 21;
    y1 = 0.9;
    t2 = 0.8;
    hy2 = 1.1;
    hy1 = 0.;
    t1 = lx;
    t3 = 0.9;
    t2 = 0.4;
    t4 = 0.5;
    break;
  case 4:
    std::cout << "processing vertex" << std::endl;
    xtitle = "number of vertices";
    scen = "vtx";
    values = "|#eta| < 2.1, p_{T} > 8 GeV/c";
    bins = bins4;
    x1 = 0.;
    x2 = 31.;
    y1 = 0.9;
    y2 = 1.1;
    hy1 = 0.75;
    hy2 = 1.1;
    lx = 0.56;
    break;
  }

  std::stringstream file;
  file << "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/MuonID_" << id.c_str() << "2012_" << scen.c_str() << absetaBin.c_str() << separation.c_str() << add.c_str() << ".root"; //_2012
  TFile *f = TFile::Open(file.str().c_str());

  TGraphAsymmErrors *g_data = (TGraphAsymmErrors *) f->Get("DATA");
  TGraphAsymmErrors *g_mc = (TGraphAsymmErrors *) f->Get("MC");
  TGraphAsymmErrors *g_ratio = (TGraphAsymmErrors *) f->Get("RATIO");

  TCanvas *c = new TCanvas("c", "c", 500, 500);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0.05);
  pad1->SetTickx();
  pad1->SetTicky();
  pad1->Draw();
  pad1->cd();

  TH1F *hFrame1 = pad1->DrawFrame(x1, hy1, x2, hy2);
  hFrame1->SetYTitle("#epsilon");
  if(scenario == 3)  hFrame1->GetYaxis()->SetTitleOffset(1.3);
  else  hFrame1->GetYaxis()->SetTitleOffset(1.6);
  hFrame1->GetXaxis()->SetLabelOffset(999);
  hFrame1->GetXaxis()->SetLabelSize(0);
  hFrame1->GetXaxis()->SetLabelFont(63); //font in pixels
  hFrame1->GetXaxis()->SetLabelSize(16); //in pixels
  hFrame1->GetYaxis()->SetLabelFont(63); //font in pixels
  hFrame1->GetYaxis()->SetLabelSize(16); //in pixels
  hFrame1->GetYaxis()->SetTitleFont(69); //font in pixels
  hFrame1->GetYaxis()->SetTitleSize(16); //in pixels
  hFrame1->GetYaxis()->SetNdivisions(510);
  hFrame1->GetYaxis()->SetDecimals();

  TLegend *l1 = new TLegend(lx, 0.13, 0.9, 0.35, values.c_str());
  l1->SetTextSize(0.05);
  //l1->SetBorderSize(0);
  l1->SetFillColor(kWhite);
  l1->SetShadowColor(0);

  TLatex *latex = new TLatex();
  latex->SetNDC(kTRUE);
  latex->SetTextSize(0.05);
  double left = 0.12, top = 0.92, /*bottom = 0.23,*/ right = 0.8;
  latex->DrawLatex(left,top,"CMS preliminary             Run 2012");

  TLatex *latex1 = new TLatex();
  latex1->SetNDC(kTRUE);
  latex1->SetTextSize(0.05);
  latex1->DrawLatex(right,top,"#sqrt{s} = 8 TeV");

  TLine *line = new TLine();
  line->SetLineStyle(3);

  for(unsigned int iBins = 0; iBins < bins.size(); iBins++){
    line->DrawLine(bins[iBins], hy1, bins[iBins], hy2);
  }
  pad1->SetGridy();

  g_data->Draw("P SAME");
  g_mc->SetMarkerStyle(20);
  l1->AddEntry(g_data, "Data", "PL");

  g_mc->Draw("P SAME");
  g_mc->SetMarkerColor(kRed);
  g_mc->SetLineColor(kRed);
  g_mc->SetMarkerStyle(22);
  l1->AddEntry(g_mc, "MC", "PL");
  l1->Draw();

  //TLatex *latex2 = new TLatex();
  //latex2->SetNDC(kTRUE);
  //latex2->SetTextSize(0.04);
  //left = 0.18;
  //if(scenario == 3) left = 0.25;
  //latex2->DrawLatex(left,bottom,values.c_str());

  TPaveText *text = new TPaveText(t1,t2,t3,t4,"NDC");
  text->SetFillColor(kWhite);
  text->SetTextSize(0.05);
  text->AddText(title.c_str());
  text->Draw();
  c->cd();

  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.28);
  pad2->SetTickx();
  pad2->SetTicky();
  pad2->Draw();
  pad2->cd();

  TH1F *hFrame2 = pad2->DrawFrame(x1, y1, x2, y2);
  hFrame2->SetXTitle(xtitle.c_str());
  hFrame2->SetYTitle("Data/MC");
  if(scenario == 3) hFrame2->GetYaxis()->SetTitleOffset(1.3);
  else hFrame2->GetYaxis()->SetTitleOffset(1.6);
  hFrame2->GetXaxis()->SetTitleOffset(3.5);
  hFrame2->GetXaxis()->SetLabelFont(63); //font in pixels
  hFrame2->GetXaxis()->SetLabelSize(16); //in pixels
  hFrame2->GetYaxis()->SetLabelFont(63); //font in pixels
  hFrame2->GetXaxis()->SetTitleFont(63); //font in pixels
  hFrame2->GetXaxis()->SetTitleSize(16); //in pixels
  hFrame2->GetYaxis()->SetLabelSize(16); //in pixels
  hFrame2->GetYaxis()->SetTitleFont(69); //font in pixels
  hFrame2->GetYaxis()->SetTitleSize(16); //in pixels
  hFrame2->GetYaxis()->SetNdivisions(5,5,0);
  hFrame2->GetYaxis()->SetDecimals();

  g_ratio->Draw("P SAME");
  g_ratio->SetMarkerColor(kBlue);
  g_ratio->SetLineColor(kBlue);
  g_ratio->SetMarkerStyle(21);

  for(unsigned int iBins = 0; iBins < bins.size(); iBins++){
    line->DrawLine(bins[iBins], y1, bins[iBins], y2);
  }
  pad2->SetGridy();

  c->cd();
  std::stringstream name;
  name << "Figures/Approval/" << id.c_str() << "ID_" << scen.c_str() << absetaBin.c_str() << separation.c_str() << "_2012.pdf";
  c->SaveAs(name.str().c_str());
}

#ifndef __CINT__
int main (int argc, char* const argv[])
{
  ratio();

  return 0;
}
#endif
