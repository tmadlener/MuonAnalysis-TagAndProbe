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
void createPtRootFiles(const std::string& datafile, const std::string& mcfile, const std::string& effName,
                       const std::string& outputfile = "createPtRootFiles_out.root")
{
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);

  //output file
  TFile *output = new TFile(outputfile.c_str(),"RECREATE");

  // Name of samples: data and MC
  const std::vector<std::string> effSampleName = {"DATA", "MC"};
  //const std::string effSampleName[] = {"MC", "MC"};
  const auto nEffSample = effSampleName.size();
  // Name of trigger: Mu5_Track2, Mu7_Track7
  // const std::vector<std::string> trackName = {"Mu5_Track2", "Mu7_Track7"};
  const std::vector<std::string> trackName = {"Mu7p5_Track2"};
  const auto nTrack = trackName.size();

  //Declare bins according to efficiency
  //pt_abseta
  const std::vector<double> pTBins = {2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 40.0};
  const std::vector<double> etaBins = {0,2.4}; // does not matter since we always look at only one abseta bin

  const auto nPTBins = pTBins.size() - 1;
  const auto nEtaBins = etaBins.size() - 1;

  // structure to store values
  storage values[nEffSample][nTrack][nPTBins][nEtaBins];

  // initialize storage
  for(size_t iEffSample = 0; iEffSample < nEffSample; iEffSample++){
    for(size_t iTrack = 0; iTrack < nTrack; iTrack++){
      for (size_t iPTBins = 0; iPTBins < nPTBins; iPTBins++){
        for (size_t iEtaBins = 0; iEtaBins < nEtaBins; iEtaBins++){
          values[iEffSample][iTrack][iPTBins][iEtaBins].null();
        } //iPTBins
      } //iEtaBins
    }//iTrack
  } //iEffSample


  for(size_t iEffSample = 0; iEffSample < nEffSample; iEffSample++){

    //create TGraphAsymmErrors to store graph
    TGraphAsymmErrors *graph = new TGraphAsymmErrors();

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
    for(size_t iTrack = 0; iTrack < nTrack; iTrack++){

      std::stringstream directory;
      directory << effName << "_pt_abseta_" << trackName[iTrack] << "_Jpsi";
      //directory << "Soft_pt_abseta";
      TDirectory* dir_run=cd(dir_tpTree, directory.str().c_str());
      if (!dir_run) return;

      //Jump to fit directory
      TDirectory* dir_fit_eff = cd(dir_run,"fit_eff_plots");
      if (!dir_fit_eff) return;
      //std:: cout << "Found fit directory" << std::endl;

      std::string plot = "pt_PLOT_" + trackName[iTrack] + "_Jpsi_TK_pass_&_tag_" + trackName[iTrack] + "_Jpsi_MU_pass";

      for(size_t iEtaBins = 0; iEtaBins < nEtaBins; iEtaBins++){

        std::stringstream name;
        name << effSampleName[iEffSample];
        std::cout << name.str().c_str() << std::endl;

        std::cout << plot.c_str() << std::endl;
        graph->SetName(name.str().c_str());
        graph->SetTitle(plot.c_str());
        TCanvas *c = dynamic_cast<TCanvas*>(dir_fit_eff->Get(plot.c_str()));
        if (!c) std::cout << "No " << plot.c_str() << " in " << effSampleName[iEffSample] << std::endl;
        TGraphAsymmErrors *get_plot = dynamic_cast<TGraphAsymmErrors*>(c->FindObject("hxy_fit_eff"));
        //check if there are values in PLOT
        int N = get_plot->GetN();
        if (N == 0) continue;

        for(int iPTBins = 0; iPTBins < N; iPTBins++){

          //get values from plot
          double x = 0, y = 0;
          if (get_plot->GetPoint(iPTBins, x, y) != iPTBins) {
            std::cout << "Error while getting point " << iPTBins << " from graph" << std::endl;
          }
          double err_high = get_plot->GetErrorYhigh(iPTBins);
          double err_low = get_plot->GetErrorYlow(iPTBins);
          double var_high = get_plot->GetErrorXhigh(iPTBins);
          double var_low = get_plot->GetErrorXlow(iPTBins);
          //if(err_high > 0.05) {err_high = err_low; std::cout << "changed high error" << std::endl;}
          //if(err_low > 0.05) {err_low = err_high; std::cout << "changed low error" << std::endl;}
          //if(iEffSample==1 && err_high > 0.05) {err_high = err_low; std::cout << "changed MC high error" << std::endl;}
          if(y + err_high > 1) err_high = 1 - y;

          //store values
          for(size_t s = 0; s < nPTBins; s++){
            if(x > pTBins[s] && x < pTBins[s+1]){
              values[iEffSample][iTrack][s][iEtaBins].setEff(y, err_low, err_high);
              values[iEffSample][iTrack][s][iEtaBins].setVar(x, var_low, var_high);
              std::cout << s << ". " << trackName[iTrack] << ": eff = " << values[iEffSample][iTrack][s][iEtaBins].eff
                        << " pt = " << values[iEffSample][iTrack][s][iEtaBins].var << std::endl;
              break;
            }
          }//s

        } //iPTBins
      } //iEtaBins
    } // iTrack


    for(size_t iEtaBins = 0; iEtaBins < nEtaBins; iEtaBins++){

      int points = 0;

      for(size_t iPTBins = 0; iPTBins < nPTBins; iPTBins++){

        //fill TGraphsAsymmErrors
        //fill with Mu5_Track2 for pt < 9 GeV - former 7 GeV
        if(/*pTBins[iPTBins] < 8 &&*/ !isZero(values[iEffSample][0][iPTBins][iEtaBins].eff)){
          graph->SetPoint(points, values[iEffSample][0][iPTBins][iEtaBins].var, values[iEffSample][0][iPTBins][iEtaBins].eff);
          graph->SetPointError(points,
                               values[iEffSample][0][iPTBins][iEtaBins].var_low, values[iEffSample][0][iPTBins][iEtaBins].var_high,
                               values[iEffSample][0][iPTBins][iEtaBins].eff_low, values[iEffSample][0][iPTBins][iEtaBins].eff_high);
          std::cout << effName << " " << effSampleName[iEffSample] << " Mu5_Track2 " << points <<" ptbin = "  << iPTBins
                    << " mean = " << values[iEffSample][0][iPTBins][iEtaBins].var << " eff = " << values[iEffSample][0][iPTBins][iEtaBins].eff
                    << " low = " << values[iEffSample][0][iPTBins][iEtaBins].eff_low << " high = " << values[iEffSample][0][iPTBins][iEtaBins].eff_high
                    << std::endl;
          points++;
        }
      } // iPTBins

      graph->SetMarkerStyle(8);
      output->cd();
      graph->Write();

    } // iEtaBins
  } // iEffSample

    // ratio: DATA/MC
  TGraphAsymmErrors *ratio = new TGraphAsymmErrors();

  for(size_t iEtaBins = 0; iEtaBins < nEtaBins; iEtaBins++){

    std::stringstream name1;
    name1 << "RATIO";
    std::cout << name1.str() << std::endl;

    int points = 0;
    for(size_t iPTBins = 0; iPTBins < nPTBins; iPTBins++){

      // compute ratio
      double eff_ratio = 0,
        err_low_ratio = 0,
        err_high_ratio = 0,
        mean_var = 0;

      if(/*pTBins[iPTBins] < 8 &&*/ !isZero(values[1][0][iPTBins][iEtaBins].eff)){ // 7 GeV
        eff_ratio = values[0][0][iPTBins][iEtaBins].eff / values[1][0][iPTBins][iEtaBins].eff;
        err_low_ratio = eff_ratio*
          TMath::Sqrt(TMath::Power(values[0][0][iPTBins][iEtaBins].eff_low/values[0][0][iPTBins][iEtaBins].eff,2)
                      + TMath::Power(values[1][0][iPTBins][iEtaBins].eff_low/values[1][0][iPTBins][iEtaBins].eff,2));
        err_high_ratio = eff_ratio*
          TMath::Sqrt(TMath::Power(values[0][0][iPTBins][iEtaBins].eff_high/values[0][0][iPTBins][iEtaBins].eff,2)
                      + TMath::Power(values[1][0][iPTBins][iEtaBins].eff_high/values[1][0][iPTBins][iEtaBins].eff,2));
        mean_var = (values[0][0][iPTBins][iEtaBins].var + values[1][0][iPTBins][iEtaBins].var)/2;
      }

      double x_high = 0, x_low = 0;
      for (size_t i = 0; i < nPTBins; i++){
        if ((pTBins[i] <= mean_var) && (pTBins[i+1] >= mean_var)){
          //std::cout << "found bin: " << pTBins[i] << " " << pTBins[i+1] << std::endl;
          x_low = mean_var - pTBins[i];
          x_high = - mean_var + pTBins[i+1];
        }
      }

      //std::cout << x_low << " " << x_high << " " << mean_var << std::endl;
      //fill TGraphsAsymmErrors
      if (!TMath::IsNaN(eff_ratio) && eff_ratio!=0){
        ratio->SetPoint(points, mean_var, eff_ratio);
        ratio->SetPointError(points, TMath::Abs(x_low), TMath::Abs(x_high), err_low_ratio, err_high_ratio);
        points++;
        //std::cout << effName[iEff] << " etabin = "  << iPTBins << " mean = " << mean_var
        //          << " eff = " << eff_ratio << " low = " << err_low_ratio << " high = " << err_high_ratio
        //          << std::endl;
      }

      std::cout.precision(4);
      //std::cout << "Only Mu5_Track2 is shown. Be careful with ratios." << std::endl;
      std::cout << " & " << pTBins[iPTBins] << " $ < p_T < $ " << pTBins[iPTBins+1] << " & $"
                << std::fixed << values[0][0][iPTBins][iEtaBins].eff << "^{+"
                << values[0][0][iPTBins][iEtaBins].eff_high << "}_{-"
                << values[0][0][iPTBins][iEtaBins].eff_low << "}$ & $"
                << values[1][0][iPTBins][iEtaBins].eff << "^{+"
                << values[1][0][iPTBins][iEtaBins].eff_high << "}_{-"
                << values[1][0][iPTBins][iEtaBins].eff_low << "}$ & $"
                << eff_ratio<< "^{+"
                << err_high_ratio << "}_{-"
                << err_low_ratio << "}$\\" << std::endl;

    } //iPTBins

    ratio->SetMarkerStyle(8);
    ratio->SetName(name1.str().c_str());
    ratio->SetTitle("Ratio Data/MC");
    output->cd();
    ratio->Write();

  } //iEtaBins

  output->Close();
} //void

#ifndef __CINT__
int main (int argc, char* const argv[])
{
  if (argc < 3) {
    std::cerr << "Need at least three input arguments: input DATA file, input MC file, ID name." << std::endl;
    std::cerr << "Optional 4th argument: output file name" << std::endl;
    return 1;
  } else if (argc < 4) {
    createPtRootFiles(argv[1], argv[2], argv[3]);
  } else {
    createPtRootFiles(argv[1], argv[2], argv[3], argv[4]);
  }

  return 0;
}
#endif
