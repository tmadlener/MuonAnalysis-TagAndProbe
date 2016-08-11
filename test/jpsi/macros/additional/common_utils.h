// header providing functionality that is used throughout the other macros
//
// written by Ilse Kraetschmer in the first place and moved to a common header by Thomas Madlener (27th July 2016)
// only slightly adapted (replaced 0 with nullptr, doc comments, passing const refs to strings instead of strings)

#ifndef COMMON_UTILS_TNP_H__
#define COMMON_UTILS_TNP_H__

#include "TFile.h"
#include "TDirectory.h"

// stl
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

/** Return opened root file or nullptr on fail. */
TFile* open(const std::string& filename)
{
  TFile* file = TFile::Open(filename.c_str());
  if ((!file)||(file->IsZombie())) {
    std::cerr << "File "<<filename<<" is bad" << std::endl;
    return nullptr;
  }
  return file;
}

/** change to a subdirectory and return that directory or nullptr on fail. */
TDirectory* cd(TDirectory* indir, const std::string& subdir)
{
  TDirectory* outdir = dynamic_cast<TDirectory*>(indir->Get(subdir.c_str()));
  if (!outdir) {
    std::cerr << "Directory "<< subdir << " not found" << std::endl;
    return nullptr;
  }
  return outdir;
}

/** comparison to 0. */
inline bool isZero(double a) { return fabs(a) < 1E-301; }

/** Struct to store efficiency values. */
struct storage {
  double var;
  double var_low;
  double var_high;
  double eff;
  double eff_low;
  double eff_high;
  /** Resets every member. */
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


/** struct to store the different settings for the different plots in one instance. */
struct PlotSettings {
  PlotSettings(const std::string& axis, const std::string& val, double x1, double x2, double y1, double y2,
               double t1, double t2, double t3, double t4, double lx, double hy1, double hy2 ) :
    xtitle(axis), values(val), x1(x1), x2(x2), y1(y1), y2(y2), t1(t1), t2(t2), t3(t3), t4(t4), lx(lx), hy1(hy1), hy2(hy2) {}

  PlotSettings() = default;

  std::string xtitle{}; /**< x-axis name. */
  std::string values{}; /**< string containing fixed parameter information. printed on pad. */
  double x1{}; /**< x-range low value. */
  double x2{}; /**< x-range high value. */
  double y1{}; /**< y-range low value ratio plot. */
  double y2{}; /**< y-range high value ratio plot. */
  double t1{}; /**< left coordinate of the text patch containing the ID. */
  double t2{}; /**< lower coordinate of the text patch containing the ID. */
  double t3{}; /**< right coordinate of the text patch containing the ID. */
  double t4{}; /**< upper coordinate of the text patch containing the ID. */
  double lx{}; /**< right boundary of legend containing the values string. */
  double hy1{}; /**< y-range low value of efficiency plot.*/
  double hy2{}; /**< y-range high value of efficiency plot. */
};

PlotSettings getDefaultSettings(int scenario, const std::string& absetaBin = "")
{
  // define some default values for better readability
  const double t1 = 0.43;
  const double t2 = 0.8;
  const double t3 = 0.63;
  const double t4 = 0.9;
  const double lx = 0.65;

  const std::vector<std::string> etaRanges = {"|#eta| < 0.9", "0.9 < |#eta| < 1.2",
                                                "1.2 < |#eta| < 2.1", "2.1 < |#eta| <2.4"};
  std::string etaRange;

  switch(scenario) {
  case 1:
    return PlotSettings("#eta", "p_{T} > 8 GeV/c", -2.5, 2.5, 0.9, 1.1, t1, t2, t3, t4, lx, 0.75, 1.1);
  case 2:
    return PlotSettings("|#eta|", "p_{T} > 8 GeV/c", 0, 2.4, 0.3, 1.3, t1, 0.6, t3, 0.7, lx, 0., 1.05);
  case 3:
    etaRange = etaRanges[std::stoi(absetaBin)];
    return PlotSettings("p_{T} [GeV/c]", etaRange, 1., 41., 0.6, 1.2, lx, 0.4, 0.9, 0.5, lx, 0., 1.1);
  case 4:
    return PlotSettings("number of vertices", "|#eta| < 2.4, p_{T} > 8 GeV/c",
                        0., 31., 0.9, 1.1, t1, t2, t3, t4, 0.56, 0.75, 1.1);
  default:
    return PlotSettings();
  }
}

#endif
