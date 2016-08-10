////////////////////////////////
//
//written by Ilse Kraetschmer
//last updated on 19th August 2013
//
////////////////////////////////
//change for combination of tracks

#include "common_utils.h"

//C/C++
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <limits>

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
#include "TMath.h"

//////////////////////////////////////////////
//MAIN FUNCTION
//////////////////////////////////////////////
void createRootFiles(){
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);

  // const std::vector<std::string> effName = {"Loose2012", "Soft2012", "newSoft2012", "Tight2012"};
  const std::vector<std::string> effName = {"Loose2016"};
  int iEff = 0;
  const std::vector<std::string> scenario = {"_eta", "_vtx", "_plateau_abseta"};
  int iScen = 0;

  //input files
  std::string datafile = "/afs/hephy.at/work/t/tmadlener/CMSSW_8_0_12/src/data_rootfiles/TnP_MuonID__data_all__" + effName[iEff] + scenario[iScen] + ".root";
  std::string mcfile = "/afs/hephy.at/work/t/tmadlener/CMSSW_8_0_12/src/mc_rootfiles/TnP_MuonID__signal_mc__" + effName[iEff] + scenario[iScen] + ".root";

  //output file
  std::string outputfile = "/afs/hephy.at/work/t/tmadlener/CMSSW_8_0_12/src/outputfiles/MuonID_" + effName[iEff] + scenario[iScen] + ".root";
  TFile *output = new TFile(outputfile.c_str(),"RECREATE");

  // Name of samples: data and MC
  const std::vector<std::string> effSampleName = {"DATA", "MC"};
  //const std::string effSampleName[] = {"MC NO TRIGGER", "MC"};
  const size_t nEffSample = effSampleName.size();

  //Declare bins according to efficiency
  //eta
  const std::vector<double> etabins = {-2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1};

  //vtx
  const std::vector<double> vtxbins = {0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 22.5, 24.5, 26.5, 28.5, 30.5};

  //plateau_abseta
  const std::vector<double> absetabins = {0., 0.9, 1.2, 2.1, 2.4}; //abseta


  std::vector <double> bins1;
  if(iScen == 0) bins1 = etabins;
  else if(iScen == 1) bins1 = vtxbins;
  else if(iScen == 2) bins1 = absetabins;
  else{
    std::cout << "This scenario is not possible!" << std::endl;
    return;}

  const auto nBins1 = bins1.size();
  size_t nBins2 = 2;

  // structure to store values
  storage values[nEffSample][nBins1][nBins2];

  // initialize storage
  for(size_t iEffSample = 0; iEffSample < nEffSample; iEffSample++){
    for (size_t iBins1 = 0; iBins1 < nBins1; iBins1++){
      for (size_t iBins2 = 0; iBins2 < nBins2; iBins2++){
        values[iEffSample][iBins1][iBins2].null();
      } //iBins1
    } //iBins2
  } //iEffSample


  for(size_t iEffSample = 0; iEffSample < nEffSample; iEffSample++){

    // open input files
    TFile *file;
    if(iEffSample==0){
      file = open(datafile.c_str());
      if(!file) return;
      std:: cout << "Sucessfully opened data file" << std::endl;
    }
    else{
      file = open(mcfile.c_str());
      if(!file) return;
      std:: cout << "Sucessfully opened MC file"<< std::endl;
    }

    //Jump to tpTree
    TDirectory* dir_tpTree=cd(file,"tpTree");
    if (!dir_tpTree) return;

    //Jump to next directory
    std::stringstream directory;
    directory << effName[iEff] << scenario[iScen];
    if (iScen == 1) {
      directory << "_Mu7p5_Track2_Jpsi"; // for nVertices this is part of the name of the TDirectory
      // TODO: check how this can be implemented cleaner
    }
    TDirectory* dir_run=cd(dir_tpTree, directory.str());
    if (!dir_run) return;

    //Jump to fit directory
    TDirectory* dir_fit_eff = cd(dir_run,"fit_eff_plots");
    if (!dir_fit_eff) return;
    //std:: cout << "Found fit directory" << std::endl;

    //create TGraphAsymmErrors
    TGraphAsymmErrors *graph = new TGraphAsymmErrors();

    std::string plot;
    if(iScen == 0)
      plot = "eta_PLOT_Mu7p5_Track2_Jpsi_TK_pass_&_tag_Mu7p5_Track2_Jpsi_MU_pass";
    else if(iScen == 1)
      plot = "tag_nVertices_PLOT_Mu7p5_Track2_Jpsi_TK_pass_&_tag_Mu7p5_Track2_Jpsi_MU_pass";
    else if(iScen == 2){
      //if(iEffSample==1)
      plot = "abseta_PLOT_Mu7p5_Track2_Jpsi_TK_pass_&_tag_Mu7p5_Track2_Jpsi_MU_pass";
      // plot = "abseta_PLOT_Mu7_Track7_Jpsi_TK_pass_&_tag_Mu7_Track7_Jpsi_MU_pass";
      //else
      //plot = "abseta_PLOT";
    }

    for(size_t iBins2 = 0; iBins2 < nBins2-1; iBins2++){

      std::stringstream name;
      name << effSampleName[iEffSample];
      std::cout << name.str() << std::endl;

      std::cout << plot << std::endl;
      graph->SetName(name.str().c_str());
      graph->SetTitle(plot.c_str());
      TCanvas *c = dynamic_cast<TCanvas*>(dir_fit_eff->Get(plot.c_str()));
      if (!c) std::cout << "No " << plot.c_str() << " in " << effSampleName[iEffSample] << std::endl;
      TGraphAsymmErrors *get_plot = dynamic_cast<TGraphAsymmErrors*>(c->FindObject("hxy_fit_eff"));
      //check if there are values in PLOT
      int N = get_plot->GetN();
      if (N == 0) continue;

      int points = 0;
      for(int iBins1 = 0; iBins1 < N; iBins1++){
        //get values from plot
        double x = 0, y = 0;
        if(get_plot->GetPoint(iBins1, x, y) != iBins1) {
          std::cout << "Error while getting point " << iBins1 << " from graph" << std::endl;
        }
        double err_high = get_plot->GetErrorYhigh(iBins1);
        double err_low = get_plot->GetErrorYlow(iBins1);
        double var_high = get_plot->GetErrorXhigh(iBins1);
        double var_low = get_plot->GetErrorXlow(iBins1);

        //if(iEffSample == 1 && err_high > 0.02) {err_high = err_low; std::cout << "changed high error" << std::endl;}
        if(y + err_high > 1) err_high = 1-y;
        //if(err_low > 0.5) err_low = 0.025;

        //store values
        values[iEffSample][iBins1][iBins2].setEff(y, err_low, err_high);
        values[iEffSample][iBins1][iBins2].setVar(x, var_low, var_high);

        //fill TGraphsAsymmErrors
        if (!isZero(values[iEffSample][iBins1][iBins2].eff) && values[iEffSample][iBins1][iBins2].eff > 0.3){
          graph->SetPoint(points, values[iEffSample][iBins1][iBins2].var, values[iEffSample][iBins1][iBins2].eff);
          graph->SetPointError(points,
                               values[iEffSample][iBins1][iBins2].var_low, values[iEffSample][iBins1][iBins2].var_high,
                               values[iEffSample][iBins1][iBins2].eff_low, values[iEffSample][iBins1][iBins2].eff_high);
          points++;
          std::cout << effName[iEff] << " " << effSampleName[iEffSample] << " etabin = "  << iBins1
                    << " mean = " << values[iEffSample][iBins1][iBins2].var << " eff = " << values[iEffSample][iBins1][iBins2].eff
                    << " low = " << values[iEffSample][iBins1][iBins2].eff_low << " high = " << values[iEffSample][iBins1][iBins2].eff_high
                    << std::endl;
        }

      } //iBins1

      graph->SetMarkerStyle(8);
      output->cd();
      graph->Write();

    } //iBins2

  } // iEffSample

    // ratio: DATA/MC
  TGraphAsymmErrors *ratio = new TGraphAsymmErrors();

  for(size_t iBins2 = 0; iBins2 < nBins2-1; iBins2++){

    std::stringstream name1;
    name1 << "RATIO";
    std::cout << name1.str().c_str() << std::endl;

    int points = 0;
    for(size_t iBins1 = 0; iBins1 < nBins1-1; iBins1++){

      // compute ratio
      double eff_ratio = values[0][iBins1][iBins2].eff / values[1][iBins1][iBins2].eff;
      double err_low_ratio = eff_ratio*
        TMath::Sqrt(TMath::Power(values[0][iBins1][iBins2].eff_low/values[0][iBins1][iBins2].eff,2)
                    + TMath::Power(values[1][iBins1][iBins2].eff_low/values[1][iBins1][iBins2].eff,2));
      double err_high_ratio = eff_ratio*
        TMath::Sqrt(TMath::Power(values[0][iBins1][iBins2].eff_high/values[0][iBins1][iBins2].eff,2)
                    + TMath::Power(values[1][iBins1][iBins2].eff_high/values[1][iBins1][iBins2].eff,2));

      double mean_var = (values[0][iBins1][iBins2].var + values[1][iBins1][iBins2].var)/2;

      double x_high = 0, x_low = 0;
      for (size_t i = 0; i < nBins1-1; i++){
        if ((bins1[i] <= mean_var) && (bins1[i+1] >= mean_var)){
          //std::cout << "found bin: " << bins1[i] << " " << bins1[i+1] << std::endl;
          x_low = mean_var - bins1[i];
          x_high = - mean_var + bins1[i+1];
        }
      }

      //std::cout << x_low << " " << x_high << " " << mean_var << std::endl;
      //fill TGraphsAsymmErrors
      if (!TMath::IsNaN(eff_ratio)){
        ratio->SetPoint(points, mean_var, eff_ratio);
        ratio->SetPointError(points, TMath::Abs(x_low), TMath::Abs(x_high), err_low_ratio, err_high_ratio);
        points++;
        //std::cout << effName[iEff] << " etabin = "  << iBins1 << " mean = " << mean_var
        //          << " eff = " << eff_ratio << " low = " << err_low_ratio << " high = " << err_high_ratio
        //          << std::endl;
      }

      std::cout.precision(4);
      std::cout << " & " << bins1[iBins1] << " $ < |eta| < $ " << bins1[iBins1+1] << " & $"
                << std::fixed << values[0][iBins1][iBins2].eff << "^{+"
                << values[0][iBins1][iBins2].eff_high << "}_{-"
                << values[0][iBins1][iBins2].eff_low << "}$ & $"
                << values[1][iBins1][iBins2].eff << "^{+"
                << values[1][iBins1][iBins2].eff_high << "}_{-"
                << values[1][iBins1][iBins2].eff_low << "}$ & $"
                // << std::fixed << values[1][iBins1][iBins2].eff << "^{+"
                // << values[1][iBins1][iBins2].eff_high << "}_{-"
                // << values[1][iBins1][iBins2].eff_low << "}$ & $"
                // << values[0][iBins1][iBins2].eff << "^{+"
                // << values[0][iBins1][iBins2].eff_high << "}_{-"
                // << values[0][iBins1][iBins2].eff_low << "}$ & $"
                << eff_ratio << "^{+"
                << err_high_ratio << "}_{-"
                << err_low_ratio << "}$\\" << std::endl;


    } //iBins1

    ratio->SetMarkerStyle(8);
    ratio->SetName(name1.str().c_str());
    ratio->SetTitle("Ratio Data/MC");
    output->cd();
    ratio->Write();

  } //iBins2

} //void

#ifndef __CINT__
int main(int argc, char* const argv[])
{
  createRootFiles();

  return 0;
}

#endif
