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
void createPtRootFiles(){
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);


  // code set up so that MC file is split in different files for abseta bins
  // const std::vector<std::string> effName = {"Loose2012", "Soft2012", "newSoft2012", "Tight2012"};
  const std::vector<std::string> effName = {"Loose2016", "Loose2016_test"};
  const int iEff = 0;
  // give number of abseta bin
  const int abseta = 0;
  const std::vector<std::string> absEtaBin = { "0_0p9", "0p9_1p2", "1p2_2p1", "2p1_2p4" };

  // input files
  const std::string datafile = "/afs/hephy.at/work/t/tmadlener/CMSSW_8_0_12/src/data_rootfiles/TnP_MuonID__data_all__" + effName[iEff] + "_pt_abseta_abseta_" + absEtaBin[abseta] + "_all_pt.root";
  const std::string mcfile = "/afs/hephy.at/work/t/tmadlener/CMSSW_8_0_12/src/mc_rootfiles/TnP_MuonID__signal_mc__" + effName[iEff] + "_pt_abseta_abseta_" + absEtaBin[abseta] + "_all_pt.root";

  //output file
  std::stringstream outputfile;
  // outputfile << "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/multiplicity/MuonID_" << effName[iEff] << "_pt_abseta" << abseta << "_2012_multiplicity.root";
  outputfile << "/afs/hephy.at/work/t/tmadlener/CMSSW_8_0_12/src/outputfiles/createPTRootFiles_" << effName[iEff] << "_pt_abseta" << abseta << "_test.root";
  TFile *output = new TFile(outputfile.str().c_str(),"RECREATE");

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
  const std::vector<double> bins1 = {2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 20.0, 30.0, 40.0};
  //double bins1[] = {2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0};
  //double bins1[] = {2, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 9.0, 11.0, 14.0, 17.0, 20.0};
  //double bins1[] = {10, 20, 25, 30, 35, 40, 50, 60, 90, 140, 300, 500};
  const std::vector<double> bins2 = {0,2.4}; // does not matter since we always look at only one abseta bin
  //double bins2[] = {0.,0.9,1.2,2.1,2.4};

  const auto nBins1 = bins1.size();
  const auto nBins2 = bins2.size();

  // structure to store values
  storage values[nEffSample][nTrack][nBins1][nBins2];

  // initialize storage
  for(size_t iEffSample = 0; iEffSample < nEffSample; iEffSample++){
    for(size_t iTrack = 0; iTrack < nTrack; iTrack++){
      for (size_t iBins1 = 0; iBins1 < nBins1; iBins1++){
        for (size_t iBins2 = 0; iBins2 < nBins2; iBins2++){
          values[iEffSample][iTrack][iBins1][iBins2].null();
        } //iBins1
      } //iBins2
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
      directory << effName[iEff] << "_pt_abseta_" << trackName[iTrack] << "_Jpsi";
      //directory << "Soft_pt_abseta";
      TDirectory* dir_run=cd(dir_tpTree, directory.str().c_str());
      if (!dir_run) return;

      //Jump to fit directory
      TDirectory* dir_fit_eff = cd(dir_run,"fit_eff_plots");
      if (!dir_fit_eff) return;
      //std:: cout << "Found fit directory" << std::endl;

      std::string plot = "pt_PLOT_" + trackName[iTrack] + "_Jpsi_TK_pass_&_tag_" + trackName[iTrack] + "_Jpsi_MU_pass";

      for(size_t iBins2 = 0; iBins2 < nBins2-1; iBins2++){

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

        for(int iBins1 = 0; iBins1 < N; iBins1++){

          //get values from plot
          double x = 0, y = 0;
          if (get_plot->GetPoint(iBins1, x, y) != iBins1) {
            std::cout << "Error while getting point " << iBins1 << " from graph" << std::endl;
          }
          double err_high = get_plot->GetErrorYhigh(iBins1);
          double err_low = get_plot->GetErrorYlow(iBins1);
          double var_high = get_plot->GetErrorXhigh(iBins1);
          double var_low = get_plot->GetErrorXlow(iBins1);
          //if(err_high > 0.05) {err_high = err_low; std::cout << "changed high error" << std::endl;}
          //if(err_low > 0.05) {err_low = err_high; std::cout << "changed low error" << std::endl;}
          //if(iEffSample==1 && err_high > 0.05) {err_high = err_low; std::cout << "changed MC high error" << std::endl;}
          if(y + err_high > 1) err_high = 1 - y;

          //store values
          for(size_t s = 0; s < nBins1; s++){
            if(x > bins1[s] && x < bins1[s+1]){
              values[iEffSample][iTrack][s][iBins2].setEff(y, err_low, err_high);
              values[iEffSample][iTrack][s][iBins2].setVar(x, var_low, var_high);
              std::cout << s << ". " << trackName[iTrack] << ": eff = " << values[iEffSample][iTrack][s][iBins2].eff
                        << " pt = " << values[iEffSample][iTrack][s][iBins2].var << std::endl;
              break;
            }
          }//s

        } //iBins1
      } //iBins2
    } // iTrack


    for(size_t iBins2 = 0; iBins2 < nBins2-1; iBins2++){

      int points = 0;

      for(size_t iBins1 = 0; iBins1 < nBins1-1; iBins1++){

        //fill TGraphsAsymmErrors
        //fill with Mu5_Track2 for pt < 9 GeV - former 7 GeV
        //if(!isZero(values[iEffSample][0][iBins1][iBins2].eff)){
        if(bins1[iBins1] < 8 && !isZero(values[iEffSample][0][iBins1][iBins2].eff)){
          if((abseta == 0 && bins1[iBins1] >= 3.5) || (abseta == 1 && bins1[iBins1] >= 2.5) || (abseta == 2 && bins1[iBins1] >= 2.0)){
            // only tight
            //if((abseta == 0 && bins1[iBins1] >= 3.5) || (abseta == 1 && bins1[iBins1] >= 3.0) || (abseta == 2 && bins1[iBins1] >= 2.0)){
            graph->SetPoint(points, values[iEffSample][0][iBins1][iBins2].var, values[iEffSample][0][iBins1][iBins2].eff);
            graph->SetPointError(points,
                                 values[iEffSample][0][iBins1][iBins2].var_low, values[iEffSample][0][iBins1][iBins2].var_high,
                                 values[iEffSample][0][iBins1][iBins2].eff_low, values[iEffSample][0][iBins1][iBins2].eff_high);
            std::cout << effName[iEff] << " " << effSampleName[iEffSample] << " Mu5_Track2 " << points <<" ptbin = "  << iBins1
                      << " mean = " << values[iEffSample][0][iBins1][iBins2].var << " eff = " << values[iEffSample][0][iBins1][iBins2].eff
                      << " low = " << values[iEffSample][0][iBins1][iBins2].eff_low << " high = " << values[iEffSample][0][iBins1][iBins2].eff_high
                      << std::endl;
            points++;
            // }
          }
        }
        // fill with Mu7_Track7 for pt > 9 GeV - former 7 GeV
        else if(bins1[iBins1] >= 8 && !isZero(values[iEffSample][1][iBins1][iBins2].eff)){
          graph->SetPoint(points, values[iEffSample][1][iBins1][iBins2].var, values[iEffSample][1][iBins1][iBins2].eff);
          graph->SetPointError(points,
                               values[iEffSample][1][iBins1][iBins2].var_low, values[iEffSample][1][iBins1][iBins2].var_high,
                               values[iEffSample][1][iBins1][iBins2].eff_low, values[iEffSample][1][iBins1][iBins2].eff_high);
          std::cout << effName[iEff] << " " << effSampleName[iEffSample] << " Mu7_Track7 " << points << " ptbin = "  << iBins1
                    << " mean = " << values[iEffSample][1][iBins1][iBins2].var << " eff = " << values[iEffSample][1][iBins1][iBins2].eff
                    << " low = " << values[iEffSample][1][iBins1][iBins2].eff_low << " high = " << values[iEffSample][1][iBins1][iBins2].eff_high
                    << std::endl;
          points++;
        }

      } // iBins1

      graph->SetMarkerStyle(8);
      output->cd();
      graph->Write();

    } // iBins2
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
      double eff_ratio = 0,
        err_low_ratio = 0,
        err_high_ratio = 0,
        mean_var = 0;

      if(bins1[iBins1] < 8 && !isZero(values[1][0][iBins1][iBins2].eff)){ // 7 GeV
        if((abseta == 0 && bins1[iBins1] >= 3.5) || (abseta == 1 && bins1[iBins1] >= 2.5) || (abseta == 2 && bins1[iBins1] >= 2.0)){
          // only tight
          //if((abseta == 0 && bins1[iBins1] >= 3.5) || (abseta == 1 && bins1[iBins1] >= 3.0) || (abseta == 2 && bins1[iBins1] >= 2.0)){
          eff_ratio = values[0][0][iBins1][iBins2].eff / values[1][0][iBins1][iBins2].eff;
          err_low_ratio = eff_ratio*
            TMath::Sqrt(TMath::Power(values[0][0][iBins1][iBins2].eff_low/values[0][0][iBins1][iBins2].eff,2)
                        + TMath::Power(values[1][0][iBins1][iBins2].eff_low/values[1][0][iBins1][iBins2].eff,2));
          err_high_ratio = eff_ratio*
            TMath::Sqrt(TMath::Power(values[0][0][iBins1][iBins2].eff_high/values[0][0][iBins1][iBins2].eff,2)
                        + TMath::Power(values[1][0][iBins1][iBins2].eff_high/values[1][0][iBins1][iBins2].eff,2));
          mean_var = (values[0][0][iBins1][iBins2].var + values[1][0][iBins1][iBins2].var)/2;
          // }
        }
      }
      else{
        eff_ratio = values[0][1][iBins1][iBins2].eff / values[1][1][iBins1][iBins2].eff;
        err_low_ratio = eff_ratio*
          TMath::Sqrt(TMath::Power(values[0][1][iBins1][iBins2].eff_low/values[0][1][iBins1][iBins2].eff,2)
                      + TMath::Power(values[1][1][iBins1][iBins2].eff_low/values[1][1][iBins1][iBins2].eff,2));
        err_high_ratio = eff_ratio*
          TMath::Sqrt(TMath::Power(values[0][1][iBins1][iBins2].eff_high/values[0][1][iBins1][iBins2].eff,2)
                      + TMath::Power(values[1][1][iBins1][iBins2].eff_high/values[1][1][iBins1][iBins2].eff,2));
        mean_var = (values[0][1][iBins1][iBins2].var + values[1][1][iBins1][iBins2].var)/2;
      }

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
      if (!TMath::IsNaN(eff_ratio) && eff_ratio!=0){
        ratio->SetPoint(points, mean_var, eff_ratio);
        ratio->SetPointError(points, TMath::Abs(x_low), TMath::Abs(x_high), err_low_ratio, err_high_ratio);
        points++;
        //std::cout << effName[iEff] << " etabin = "  << iBins1 << " mean = " << mean_var
        //          << " eff = " << eff_ratio << " low = " << err_low_ratio << " high = " << err_high_ratio
        //          << std::endl;
      }

      std::cout.precision(4);
      //std::cout << "Only Mu5_Track2 is shown. Be careful with ratios." << std::endl;
      std::cout << " & " << bins1[iBins1] << " $ < p_T < $ " << bins1[iBins1+1] << " & $"
                << std::fixed << values[0][0][iBins1][iBins2].eff << "^{+"
                << values[0][0][iBins1][iBins2].eff_high << "}_{-"
                << values[0][0][iBins1][iBins2].eff_low << "}$ & $"
                << values[1][0][iBins1][iBins2].eff << "^{+"
                << values[1][0][iBins1][iBins2].eff_high << "}_{-"
                << values[1][0][iBins1][iBins2].eff_low << "}$ & $"
                << eff_ratio<< "^{+"
                << err_high_ratio << "}_{-"
                << err_low_ratio << "}$\\" << std::endl;

    } //iBins1

    ratio->SetMarkerStyle(8);
    ratio->SetName(name1.str().c_str());
    ratio->SetTitle("Ratio Data/MC");
    output->cd();
    ratio->Write();

  } //iBins2

  output->Close();
} //void

#ifndef __CINT__
int main (int argc, char* const argv[])
{
  createPtRootFiles();

  return 0;
}
#endif
