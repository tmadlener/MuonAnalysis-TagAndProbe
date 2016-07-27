////////////////////////////////
//
//written by Ilse Kraetschmer
//last updated on 19th August 2013
//
////////////////////////////////

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
void calcSF(){
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);

  const std::string effName[] = {"Loose2012", "Soft2012", "newSoft2012", "Tight2012"};
  int iEff = 3;
  const std::string scenario[] = {"_eta", "_vtx", "_plateau_abseta"};
  int iScen = 0;

  //input files
  std::stringstream datafile, mcfile, notfile;
  mcfile << "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/multiplicity/TnP_MuonID_signal_mc_" << effName[iEff] << scenario[iScen] << "_multiplicity.root";
  notfile << "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/multiplicity/TnP_MuonID_signal_mc_" << effName[iEff] << scenario[iScen] <<  "_noTrigger_multiplicity.root";

  //output file
  std::stringstream outputfile;
  outputfile << "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/multiplicity/SF_MuonID_" << effName[iEff] << scenario[iScen] << "_2012_multiplicity.root";
  TFile *output = new TFile(outputfile.str().c_str(),"RECREATE");

  // Name of samples: data and MC
  const std::string effSampleName[] = {"MC", "MC NO TRIGGER"};
  const int nEffSample = sizeof(effSampleName)/sizeof(effSampleName[0]);

  //Declare bins according to efficiency
  //eta
  double etabins[] = {-2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1};
  int nEtabins = sizeof(etabins)/sizeof(etabins[0]);
  //vtx
  double vtxbins[] ={0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 22.5, 24.5, 26.5, 28.5, 30.5}; //vtx
  int nVtxbins = sizeof(vtxbins)/sizeof(vtxbins[0]);
  //plateau_abseta
  double absetabins[] = {0., 0.9, 1.2, 2.1}; //abseta
  int nAbsetabins = sizeof(absetabins)/sizeof(absetabins[0]);

  std::vector <double> bins1;
  if(iScen == 0) bins1.assign(etabins, etabins + nEtabins);
  else if(iScen == 1) bins1.assign(vtxbins, vtxbins + nVtxbins);
  else if(iScen == 2) bins1.assign(absetabins, absetabins + nAbsetabins);
  else{
    std::cout << "This scenario is not possible!" << std::endl;
    return;}

  const int nBins1 = bins1.size();
  int nBins2 = 2;

  // structure to store values
  storage values[nEffSample][nBins1][nBins2];

  // initialize storage
  for(int iEffSample = 0; iEffSample < nEffSample; iEffSample++){
    for (int iBins1 = 0; iBins1 < nBins1; iBins1++){
      for (int iBins2 = 0; iBins2 < nBins2; iBins2++){
        values[iEffSample][iBins1][iBins2].null();
      } //iBins1
    } //iBins2
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
    std::stringstream directory;
    directory << effName[iEff] << scenario[iScen];
    TDirectory* dir_run=cd(dir_tpTree, directory.str().c_str());
    if (!dir_run) return;

    //Jump to fit directory
    TDirectory* dir_fit_eff = cd(dir_run,"fit_eff_plots");
    if (!dir_fit_eff) return;
    //std:: cout << "Found fit directory" << std::endl;

    //create TGraphAsymmErrors
    TGraphAsymmErrors *graph = new TGraphAsymmErrors();

    std::string plot;
    if(iScen == 0){
      if(iEffSample==0)
        plot = "eta_PLOT_Mu7_Track7_Jpsi_TK_pass_&_tag_Mu7_Track7_Jpsi_MU_pass";
      else
        plot = "eta_PLOT";
    }
    else if(iScen == 1){
      if(iEffSample==0)
        plot = "tag_nVertices_PLOT_Mu7_Track7_Jpsi_TK_pass_&_tag_Mu7_Track7_Jpsi_MU_pass";
      else
        plot = "tag_nVertices_PLOT";
    }
    else if(iScen == 2){
      if(iEffSample==0)
        plot = "abseta_PLOT_Mu7_Track7_Jpsi_TK_pass_&_tag_Mu7_Track7_Jpsi_MU_pass";
      else
        plot = "abseta_PLOT";
    }

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

        //if(iEffSample == 0 && err_high > 0.02) {err_high = err_low; std::cout << "changed high error" << std::endl;}
        //if(iEffSample == 1 && err_high > 0.01) {err_high = err_low; std::cout << "changed high error" << std::endl;}

        //store values
        for(int s = 0; s < nBins1; s++){
          if(x > bins1[s] && x < bins1[s+1]){
            values[iEffSample][s][iBins2].setEff(y, err_low, err_high);
            values[iEffSample][s][iBins2].setVar(x, var_low, var_high);
            std::cout << s << ": eff = " << values[iEffSample][s][iBins2].eff
                      << " pt = " << values[iEffSample][s][iBins2].var << std::endl;
            break;
          }
        }//s

      } //iBins1
    } //iBins2
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
        if(!isZero(values[iEffSample][iBins1][iBins2].eff)){
          eff_SF = values[1][iBins1][iBins2].eff / values[iEffSample][iBins1][iBins2].eff;
          err_low_SF = eff_SF*
            TMath::Sqrt(TMath::Power(values[1][iBins1][iBins2].eff_low/values[1][iBins1][iBins2].eff,2)
                        + TMath::Power(values[iEffSample][iBins1][iBins2].eff_low/values[iEffSample][iBins1][iBins2].eff,2));
          err_high_SF = eff_SF*
            TMath::Sqrt(TMath::Power(values[1][iBins1][iBins2].eff_high/values[1][iBins1][iBins2].eff,2)
                        + TMath::Power(values[iEffSample][iBins1][iBins2].eff_high/values[iEffSample][iBins1][iBins2].eff,2));

          graph->SetPoint(points, values[iEffSample][iBins1][iBins2].var, eff_SF);
          graph->SetPointError(points,
                               values[iEffSample][iBins1][iBins2].var_low, values[iEffSample][iBins1][iBins2].var_high,
                               err_low_SF, err_high_SF);
          points++;
        }

        cout.precision(4);
        std::string variable;
        if(iScen==0)
          variable = " $ < eta < $ ";
        else if(iScen==1)
          variable = " $ < nb of vertices < $ ";
        else
          variable = " $ < |eta| < $ ";
        std::cout << " & " << bins1[iBins1] << variable.c_str() << bins1[iBins1+1] << " & $"
                  << fixed //<< values[iEffSample][iBins1][iBins2].eff << "^{+"
          //<< values[iEffSample][iBins1][iBins2].eff_high << "}_{-"
          //    << values[iEffSample][iBins1][iBins2].eff_low << "}$ & $"
                  << eff_SF << "^{+"
                  << err_high_SF << "}_{-"
                  << err_low_SF << "}$\\" << std::endl;


      } // iBins1

      graph->SetMarkerStyle(8);
      output->cd();
      graph->Write();

    } // iBins2
  } // iEffSample


} //void
