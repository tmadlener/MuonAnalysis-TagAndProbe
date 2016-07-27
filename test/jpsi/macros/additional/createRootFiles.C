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
#include <iomanip>
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
void createRootFiles(){
    gROOT->SetStyle("Plain");
    gStyle->SetTitleBorderSize(0);

    const std::string effName[] = {"Loose2012", "Soft2012", "newSoft2012", "Tight2012"};
    int iEff = 0;
    const std::string scenario[] = {"_eta", "_vtx", "_plateau_abseta"};
    int iScen = 1;

    //input files
    std::stringstream datafile, mcfile;
    datafile << "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/multiplicity/TnP_MuonID_data_all_" << effName[iEff] << scenario[iScen] << "_multiplicity.root";
    mcfile << "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/multiplicity/TnP_MuonID_signal_mc_" << effName[iEff] << scenario[iScen] << "_multiplicity.root";

    //output file
    std::stringstream outputfile;
    outputfile << "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/multiplicity/MuonID_" << effName[iEff] << scenario[iScen] << "_multiplicity.root";
    TFile *output = new TFile(outputfile.str().c_str(),"RECREATE");

    // Name of samples: data and MC
    const std::string effSampleName[] = {"DATA", "MC"};
    //const std::string effSampleName[] = {"MC NO TRIGGER", "MC"};
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
            file = open(datafile.str().c_str());
            std:: cout << "Sucessfully opened data file" << std::endl;
            if(!file) return;
        }
        else{
            file = open(mcfile.str().c_str());
            if(!file) return;
            std:: cout << "Sucessfully opened MC file"<< std::endl;
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
        if(iScen == 0)
            plot = "eta_PLOT_Mu7_Track7_Jpsi_TK_pass_&_tag_Mu7_Track7_Jpsi_MU_pass";
        else if(iScen == 1)
            plot = "tag_nVertices_PLOT_Mu7_Track7_Jpsi_TK_pass_&_tag_Mu7_Track7_Jpsi_MU_pass";
        else if(iScen == 2){
            //if(iEffSample==1)
                //plot = "abseta_PLOT_Mu5_Track2_Jpsi_TK_pass_&_tag_Mu5_Track2_Jpsi_MU_pass";
                plot = "abseta_PLOT_Mu7_Track7_Jpsi_TK_pass_&_tag_Mu7_Track7_Jpsi_MU_pass";
                //else
                //plot = "abseta_PLOT";
        }

        for(int iBins2 = 0; iBins2 < nBins2-1; iBins2++){

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

            int points = 0;
            for(int iBins1 = 0; iBins1 < nBins1-1; iBins1++){

                //get values from plot
                double x = 0, y = 0;
                double z = get_plot->GetPoint(iBins1, x, y);
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

    for(int iBins2 = 0; iBins2 < nBins2-1; iBins2++){

        std::stringstream name1;
        name1 << "RATIO";
        std::cout << name1.str().c_str() << std::endl;

        int points = 0;
        for(int iBins1 = 0; iBins1 < nBins1-1; iBins1++){

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
            for (int i = 0; i < nBins1-1; i++){
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

            cout.precision(4);
            std::cout << " & " << bins1[iBins1] << " $ < |eta| < $ " << bins1[iBins1+1] << " & $"
                      << fixed << values[0][iBins1][iBins2].eff << "^{+"
                      << values[0][iBins1][iBins2].eff_high << "}_{-"
                      << values[0][iBins1][iBins2].eff_low << "}$ & $"
                      << values[1][iBins1][iBins2].eff << "^{+"
                      << values[1][iBins1][iBins2].eff_high << "}_{-"
                      << values[1][iBins1][iBins2].eff_low << "}$ & $"
//                      << fixed << values[1][iBins1][iBins2].eff << "^{+"
//                      << values[1][iBins1][iBins2].eff_high << "}_{-"
//                      << values[1][iBins1][iBins2].eff_low << "}$ & $"
//                      << values[0][iBins1][iBins2].eff << "^{+"
//                      << values[0][iBins1][iBins2].eff_high << "}_{-"
//                      << values[0][iBins1][iBins2].eff_low << "}$ & $"
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
