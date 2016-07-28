////////////////////////////////
//
//written by Ilse Kraetschmer
//last updated on 28th October 2013
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
void systCheck(){
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);

  const std::vector<std::string> effName = {"Loose2012", "Soft2012", "newSoft2012", "Tight2012"};
  int iEff = 0;

  //input files
  std::stringstream datafile, mcfile, datafile1, mcfile1;
  datafile << "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/MuonID_" << effName[iEff] << "_plateau_abseta" << "_2012.root";
  datafile1 << "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/MuonID_" << effName[iEff] << "_plateau_abseta" << "_seagulls_dist200_2012.root";

  // Name of samples: data and MC
  const std::vector<std::string> effSampleName = {"DATA", "MC", "RATIO"};
  const auto nEffSample = effSampleName.size();
  std::cout << "Number of samples: " << nEffSample << std::endl;

  //Declare bins according to efficiency
  //pt_abseta
  const std::vector<double> bins1 = {0.,0.9,1.2,2.1};
  const auto nBins1 = bins1.size();

  // structure to store values
  storage values[2][nEffSample][nBins1];

  // initialize storage
  for(size_t i = 0; i < 2; i++){ // loop through sample & refernce sample
    for(size_t iEffSample = 0; iEffSample < nEffSample; iEffSample++){
      for (size_t iBins1 = 0; iBins1 < nBins1; iBins1++){
        values[i][iEffSample][iBins1].null();
      } //iBins1
    } //iEffSample
  }//i

  for(size_t i = 0; i < 2; i++){ // loop through sample & refernce sample
    for(size_t iEffSample = 0; iEffSample < nEffSample; iEffSample++){

      // open input files
      TFile *file;
      if(i==0)
        file = open(datafile.str().c_str());
      else
        file = open(datafile1.str().c_str());
      if(!file) return;
      std:: cout << "Sucessfully opened data file" << std::endl;

      std::string plot;
      if(iEffSample==0) plot = "DATA";
      else if(iEffSample==1) plot = "MC";
      else plot = "RATIO";
      std::cout << plot.c_str() << std::endl;

      TGraphAsymmErrors *get_plot = (TGraphAsymmErrors*)gDirectory->Get(plot.c_str());
      //check if there are values in PLOT
      int N = get_plot->GetN();
      if (N == 0) continue;

      for(size_t iBins1 = 0; iBins1 < nBins1-1; iBins1++){

        //get values from plot
        double x = 0, y = 0;
        // double z = get_plot->GetPoint(iBins1, x, y);
        double err_high = get_plot->GetErrorYhigh(iBins1);
        double err_low = get_plot->GetErrorYlow(iBins1);
        double var_high = get_plot->GetErrorXhigh(iBins1);
        double var_low = get_plot->GetErrorXlow(iBins1);

        //store values
        for(size_t s = 0; s < nBins1; s++){
          if(x > bins1[s] && x < bins1[s+1]){
            values[i][iEffSample][s].setEff(y, err_low, err_high);
            values[i][iEffSample][s].setVar(x, var_low, var_high);
            std::cout << s << ": eff = " << values[i][iEffSample][s].eff
                      << " pt = " << values[i][iEffSample][s].var << std::endl;
            break;
          }
        }//s

      } //iBins1
    } // iEffSample
  }//i

  for(size_t iEffSample = 0; iEffSample < nEffSample; iEffSample++){
    // TGraphAsymmErrors *graph = new TGraphAsymmErrors();

    std::stringstream name;
    name << effSampleName[iEffSample];
    std::cout << name.str().c_str() << std::endl;

    for(size_t iBins1 = 0; iBins1 < nBins1-1; iBins1++){

      double syst = (values[1][iEffSample][iBins1].eff - values[0][iEffSample][iBins1].eff)/values[0][iEffSample][iBins1].eff;
      std::cout.precision(4);
      std::cout << bins1[iBins1] << " $ < p_T < $ " << bins1[iBins1+1] << " & "
                << std::fixed << syst << std::endl;

    } //iBins1
  } //iEffSample
} //void

#ifndef __CINT__
int main(int argc, char* const argv[])
{
  systCheck();

  return 0;
}

#endif
