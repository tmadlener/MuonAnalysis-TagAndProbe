////////////////////////////////
//
//written by Ilse Kraetschmer
//last updated on 19th August 2013
//
////////////////////////////////
//change for combination of tracks

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

//Return opened root file (or 0 on fail)
TFile* open(std::string filename){
  TFile* file = new TFile(filename.c_str());
  if ((!file)||(file->IsZombie())) {
    std::cerr << "File "<<filename<<" is bad" << std::endl;
    return 0;
  }
  return file;
}

//Change to a subdirectory
//Return: directory (or 0 on fail)
TDirectory* cd(TDirectory* indir,std::string subdir)
{
  TDirectory* outdir = dynamic_cast<TDirectory*>(indir->Get(subdir.c_str()));
  if (!outdir) {
    std::cerr << "Directory "<< subdir << " not found" << std::endl;
    return 0;
  }
  return outdir;
}

//comparison to 0
bool isZero(double a) {return fabs(a)<1E-301; }

//Struct to store efficiency values
struct storage{
  double var;
  double var_low;
  double var_high;
  double eff;
  double eff_low;
  double eff_high;
  // Resets every member
  void null() {
    var = 0.0;
    var_low = 0.0;
    var_high = 0.0;
    eff = 0.0;
    eff_low = 0.0;
    eff_high = 0.0;
  }
  void print() {
    std::cout << "var = " << var << " + " << var_high << " - " << var_low << "\n"
              << "eff = " << eff << " + " << eff_high << " - " << eff_low << "\n"
              << std::endl;
  }
  void setEff(double e, double low, double high) {
    eff = e;
    eff_low = low;
    eff_high = high;
  }
  void setVar(double meanVar, double v_low, double v_high) {
    var = meanVar;
    var_low = v_low;
    var_high = v_high;
  }
};

//////////////////////////////////////////////
//MAIN FUNCTION
//////////////////////////////////////////////
void calcPtSF(){
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);

  // code set up so that it takes different mc files for different abseta bins as input
  const std::string effName[] = {"Loose2012", "Soft2012", "newSoft2012", "Tight2012"};
  int iEff = 3;
  // give number of abseta bin
  int abseta = 2;

  //input files
  std::stringstream datafile, mcfile, notfile;
  mcfile << "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/multiplicity/TnP_MuonID_signal_mc_" << effName[iEff] <<  "_pt_abseta" << abseta << "_multiplicity.root";
  notfile << "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/multiplicity/TnP_MuonID_signal_mc_" << effName[iEff] <<  "_pt_abseta" << abseta << "_multiplicity_noTrigger.root";

  //output file
  std::stringstream outputfile;
  outputfile << "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/multiplicity/SF_MuonID_" << effName[iEff] << "_pt_abseta" << abseta << "_multiplicity.root";
  TFile *output = new TFile(outputfile.str().c_str(),"RECREATE");

  // Name of samples: data and MC
  const std::string effSampleName[] = {"MC", "MC NO TRIGGER"};
  const int nEffSample = sizeof(effSampleName)/sizeof(effSampleName[0]);
  // Name of trigger: Mu5_Track2, Mu7_Track7
  const std::string trackName[] = {"Mu5_Track2", "Mu7_Track7"};
  const int nTrack = sizeof(trackName)/sizeof(trackName[0]);

  //Declare bins according to efficiency
  //pt_abseta
  //double bins1[] = {2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 20.0};
  double bins1[] = {2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0};
  //double bins1[] = {2, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 9.0, 11.0, 14.0, 17.0, 20.0};
  //double bins1[] = {2, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 9.0, 11.0, 14.0, 17.0, 20.0, 25.0, 30.0, 35.0, 40.0};
  //double bins1[] = {10, 20, 25, 30, 35, 40, 50, 60, 90, 140, 300, 500};
  double bins2[] = {0,2.4}; // does not matter since we always look at only one abseta bin
  //double bins2[] = {0.,0.9,1.2,2.1,2.4};

  const int nBins1 = sizeof(bins1)/sizeof(bins1[0]);
  const int nBins2 = sizeof(bins2)/sizeof(bins2[0]);

  // structure to store values
  storage values[nEffSample][nTrack][nBins1][nBins2];

  // initialize storage
  for(int iEffSample = 0; iEffSample < nEffSample; iEffSample++){
    for(int iTrack = 0; iTrack < nTrack; iTrack++){
      for (int iBins1 = 0; iBins1 < nBins1; iBins1++){
        for (int iBins2 = 0; iBins2 < nBins2; iBins2++){
          values[iEffSample][iTrack][iBins1][iBins2].null();
        } //iBins1
      } //iBins2
    }//iTrack
  } //iEffSample


  for(int iEffSample = 0; iEffSample < nEffSample; iEffSample++){

    // open input files
    TFile *file;
    if(iEffSample==0){
      file = open(mcfile.str().c_str());
      if(!file) return;
      std:: cout << "Sucessfully opened MC file"<< std::endl;
    }
    else{
      file = open(notfile.str().c_str());
      if(!file) return;
      std:: cout << "Sucessfully opened MC no trigger file"<< std::endl;
    }

    //Jump to tpTree
    TDirectory* dir_tpTree=cd(file,"tpTree");
    if (!dir_tpTree) return;

    //Jump to next directory
    for(int iTrack = 0; iTrack < nTrack; iTrack++){

      std::stringstream directory;
      directory << effName[iEff] << "_pt_abseta_" << trackName[iTrack];
      //directory << "Soft_pt_abseta";
      TDirectory* dir_run=cd(dir_tpTree, directory.str().c_str());
      if (!dir_run) return;

      //Jump to fit directory
      TDirectory* dir_fit_eff = cd(dir_run,"fit_eff_plots");
      if (!dir_fit_eff) return;
      //std:: cout << "Found fit directory" << std::endl;

      std::string plot;
      std::stringstream inter;
      if(iEffSample==1)
        inter << "pt_PLOT";
      //inter << "pt_PLOT_abseta_bin" << abseta;
      else
        inter << "pt_PLOT_" << trackName[iTrack] << "_Jpsi_TK_pass_&_tag_" << trackName[iTrack] << "_Jpsi_MU_pass";
      //inter << "pt_PLOT_abseta_bin" << abseta << "_&_" << trackName[iTrack] << "_Jpsi_TK_pass_&_tag_" << trackName[iTrack] << "_Jpsi_MU_pass";
      plot = inter.str();

      for(int iBins2 = 0; iBins2 < nBins2-1; iBins2++){

        std::cout << plot.c_str() << std::endl;
        TCanvas *c = dynamic_cast<TCanvas*>(dir_fit_eff->Get(plot.c_str()));
        if (!c) std::cout << "No " << plot.c_str() << " in " << effSampleName[iEffSample] << std::endl;
        TGraphAsymmErrors *get_plot = dynamic_cast<TGraphAsymmErrors*>(c->FindObject("hxy_fit_eff"));
        //check if there are values in PLOT
        int N = get_plot->GetN();
        if (N == 0) continue;

        for(int iBins1 = 0; iBins1 < nBins1-1; iBins1++){

          //get values from plot
          double x = 0, y = 0;
          double z = get_plot->GetPoint(iBins1, x, y);
          double err_high = get_plot->GetErrorYhigh(iBins1);
          double err_low = get_plot->GetErrorYlow(iBins1);
          double var_high = get_plot->GetErrorXhigh(iBins1);
          double var_low = get_plot->GetErrorXlow(iBins1);
          if(err_high > 0.05) {err_high = err_low; std::cout << "changed high error" << std::endl;}
          if(err_low > 0.05) err_low = err_high;

          //store values
          for(int s = 0; s < nBins1; s++){
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
  } //iEffSample

    // compute scale factors
  for(int iEffSample = 0; iEffSample < nEffSample-1; iEffSample++){

    //create TGraphAsymmErrors to store SF
    TGraphAsymmErrors *graph = new TGraphAsymmErrors();
    std::string name = "SF";
    graph->SetName(name.c_str());
    graph->SetTitle(name.c_str());

    for(int iBins2 = 0; iBins2 < nBins2-1; iBins2++){

      double eff_SF = 0,
        err_low_SF = 0,
        err_high_SF = 0;
      int points = 0;

      for(int iBins1 = 0; iBins1 < nBins1-1; iBins1++){

        //fill TGraphsAsymmErrors
        //fill with Mu5_Track2 for pt < 9 GeV - former 7 GeV
        if(bins1[iBins1] < 8 && !isZero(values[iEffSample][0][iBins1][iBins2].eff)){
          if((abseta == 0 && bins1[iBins1] >= 3.5) || (abseta == 1 && bins1[iBins1] >= 2.5) || (abseta == 2 && bins1[iBins1] >= 2.0)){
            //only tight
            //if((abseta == 0 && bins1[iBins1] >= 3.5) || (abseta == 1 && bins1[iBins1] >= 3.0) || (abseta == 2 && bins1[iBins1] >= 2.0)){
            eff_SF = values[1][0][iBins1][iBins2].eff / values[iEffSample][0][iBins1][iBins2].eff;
            err_low_SF = eff_SF*
              TMath::Sqrt(TMath::Power(values[1][0][iBins1][iBins2].eff_low/values[1][0][iBins1][iBins2].eff,2)
                          + TMath::Power(values[iEffSample][0][iBins1][iBins2].eff_low/values[iEffSample][0][iBins1][iBins2].eff,2));
            err_high_SF = eff_SF*
              TMath::Sqrt(TMath::Power(values[1][0][iBins1][iBins2].eff_high/values[1][0][iBins1][iBins2].eff,2)
                          + TMath::Power(values[iEffSample][0][iBins1][iBins2].eff_high/values[iEffSample][0][iBins1][iBins2].eff,2));
            std::cout << "err_high(MC) = " << values[0][0][iBins1][iBins2].eff_high << ", err_high(noTR) = " << values[1][0][iBins1][iBins2].eff_high
                      << ", err_SF = " << err_high_SF << std::endl;
            if(err_high_SF > 0.05){
              std::cout << "changed " << err_high_SF << " to " << err_low_SF << std::endl;
              err_high_SF = err_low_SF;
            }

            graph->SetPoint(points, values[iEffSample][0][iBins1][iBins2].var, eff_SF);
            graph->SetPointError(points,
                                 values[iEffSample][0][iBins1][iBins2].var_low, values[iEffSample][0][iBins1][iBins2].var_high,
                                 err_low_SF, err_high_SF);
            std::cout << effName[iEff] << " " << effSampleName[iEffSample] << " Mu5_Track2 " << points <<" ptbin = "  << iBins1
                      << " mean = " << values[iEffSample][0][iBins1][iBins2].var << "\n DATA: eff = " << values[iEffSample][0][iBins1][iBins2].eff
                      << " low = " << values[iEffSample][0][iBins1][iBins2].eff_low << " high = " << values[iEffSample][0][iBins1][iBins2].eff_high
                      << "\n NOTR: eff = " <<  values[1][0][iBins1][iBins2].eff
                      << " low = " << values[1][0][iBins1][iBins2].eff_low << " high = " << values[1][0][iBins1][iBins2].eff_high
                      << std::endl;
            points++;
            //}
          }
        }
        // fill with Mu7_Track7 for pt > 9 GeV - former 7 GeV
        else if(bins1[iBins1] >= 8 && !isZero(values[iEffSample][1][iBins1][iBins2].eff)){
          eff_SF = values[1][1][iBins1][iBins2].eff / values[iEffSample][1][iBins1][iBins2].eff;
          err_low_SF = eff_SF*
            TMath::Sqrt(TMath::Power(values[1][1][iBins1][iBins2].eff_low/values[1][1][iBins1][iBins2].eff,2)
                        + TMath::Power(values[iEffSample][1][iBins1][iBins2].eff_low/values[iEffSample][1][iBins1][iBins2].eff,2));
          err_high_SF = eff_SF*
            TMath::Sqrt(TMath::Power(values[1][1][iBins1][iBins2].eff_high/values[1][1][iBins1][iBins2].eff,2)
                        + TMath::Power(values[iEffSample][1][iBins1][iBins2].eff_high/values[iEffSample][1][iBins1][iBins2].eff,2));

          std::cout << "err_high(MC) = " << values[0][0][iBins1][iBins2].eff_high << ", err_high(noTR) = " << values[1][0][iBins1][iBins2].eff_high
                    << ", err_SF = " << err_high_SF << std::endl;
          if(err_high_SF > 0.05){
            err_high_SF = err_low_SF;
            std::cout << "changed " << err_high_SF << " to " << err_low_SF << std::endl;
          }

          graph->SetPoint(points, values[iEffSample][1][iBins1][iBins2].var, eff_SF);
          graph->SetPointError(points,
                               values[iEffSample][1][iBins1][iBins2].var_low, values[iEffSample][1][iBins1][iBins2].var_high,
                               err_low_SF, err_high_SF);
          std::cout << effName[iEff] << " " << effSampleName[iEffSample] << " Mu7_Track7 " << points <<" ptbin = "  << iBins1
                    << " mean = " << values[iEffSample][1][iBins1][iBins2].var << "\n DATA: eff = " << values[iEffSample][1][iBins1][iBins2].eff
                    << " low = " << values[iEffSample][1][iBins1][iBins2].eff_low << " high = " << values[iEffSample][1][iBins1][iBins2].eff_high
                    << "\n NOTR: eff = " <<  values[1][1][iBins1][iBins2].eff
                    << " low = " << values[1][1][iBins1][iBins2].eff_low << " high = " << values[1][1][iBins1][iBins2].eff_high
                    << std::endl;
          points++;
        }

        //cout.precision(4);
        //std::cout << " & " << bins1[iBins1] << " $ < p_T < $ " << bins1[iBins1+1] << " & $"
        //          << fixed << values[iEffSample][0][iBins1][iBins2].eff << "^{+"
        //          << values[iEffSample][0][iBins1][iBins2].eff_high << "}_{-"
        //          << values[iEffSample][0][iBins1][iBins2].eff_low << "}$ & $"
        //          << values[1][0][iBins1][iBins2].eff << "^{+"
        //          << values[1][0][iBins1][iBins2].eff_high << "}_{-"
        //          << values[1][0][iBins1][iBins2].eff_low << "}$ & $"
        //          << eff_SF << "^{+"
        //          << err_high_SF << "}_{-"
        //          << err_low_SF << "}$\\" << std::endl;


      } // iBins1

      graph->SetMarkerStyle(8);
      output->cd();
      graph->Write();

    } // iBins2
  } // iEffSample


} //void
