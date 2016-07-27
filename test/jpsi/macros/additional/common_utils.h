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

#endif
