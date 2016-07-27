#!/usr/bin/env python
import sys
import os
import ROOT
from array import *
import math
import pickle

def getBin(binning, x):
    for i in range(len(binning)):
        if x>binning[i] and x<binning[i+1] :
            return [binning[i],binning[i+1]]

def getFile( file, dict, rootoutput, cat = "Loose", etaPtVtx = "eta", subCatSuffix = ""):
    nan = float('nan')
    f = ROOT.TFile.Open(file)
    print "opened file"
    h_data = f.Get("DATA")
   # print "found DATA"
    h_mc = f.Get("MC")
   # print "found MC"
    ratio = f.Get("RATIO")
   # print "found RATIO"
    SF = f.Get("TriggerBias_MC")

    if not dict.has_key(cat): dict[cat] = {}
    if not dict[cat].has_key(etaPtVtx+subCatSuffix): dict[cat][etaPtVtx+subCatSuffix] = {}

    etabins = [-2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1]
    ptbins = [2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0]
    ptbins1 = [2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0,20.0]
    #ptbins = [2, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 9.0, 11.0, 14.0, 17.0, 20.0]
    vtxbins = [0.5,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5,24.5,26.5,28.5,30.5]
    absetabins = [0, 0.9, 1.2, 2.1]

    if etaPtVtx == "eta": binning = etabins
    if etaPtVtx == "pt_abseta":  binning = ptbins1
    else if etaPtVtx == "pt_abseta0":  binning = ptbins
    if etaPtVtx == "vtx": binning = vtxbins
    if etaPtVtx == "abseta": binning = absetabins


    remove_data = []
    for i in range(h_data.GetN()):
        x = ROOT.Double(0)
        y = ROOT.Double(0)
        h_data.GetPoint(i, x, y)
        y_errlo = h_data.GetErrorYlow(i)
        y_errhi = h_data.GetErrorYhigh(i)

        print "x = ", x, " y = ", y

        [bl,bh] = getBin(binning, x.real)
        if not dict[cat][etaPtVtx+subCatSuffix].has_key(str(bl)+"_"+str(bh)): dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)] = {}
        if not dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)].has_key("data"): dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data"] = {}
        if not ((y == 0.0) and (y_errlo == 0) and (y_errhi == 0)):
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data"][etaPtVtx] = x.real
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data"]["efficiency"] = y.real
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data"]["err_low"] = y_errlo
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data"]["err_hi"] = y_errhi
        else:
            remove_data.append(i)
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data"][etaPtVtx] = nan
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data"]["efficiency"] = nan
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data"]["err_low"] = nan
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data"]["err_hi"] = nan

    remove_mc = []
    for i in range(h_mc.GetN()):
        x = ROOT.Double(0)
        y = ROOT.Double(0)
        h_mc.GetPoint(i, x, y)
        y_errlo = h_mc.GetErrorYlow(i)
        y_errhi = h_mc.GetErrorYhigh(i)

        print "x = ", x, " y = ", y

        [bl,bh] = getBin(binning, x)
        if not dict[cat][etaPtVtx+subCatSuffix].has_key(str(bl)+"_"+str(bh)): dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)] = {}
        if not dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)].has_key("mc"): dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["mc"] = {}
        if not ((y == 0.0001) and (y_errlo == 0.0001) and (y_errhi == 0.0001)):
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["mc"][etaPtVtx] = x.real
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["mc"]["efficiency"] = y.real
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["mc"]["err_low"] = y_errlo
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["mc"]["err_hi"] = y_errhi
        else:
            remove_mc.append(i)
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["mc"][etaPtVtx] = nan
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["mc"]["efficiency"] = nan
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["mc"]["err_low"] = nan
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["mc"]["err_hi"] = nan

    remove_ratio = []
    for i in range(ratio.GetN()):
        x = ROOT.Double(0)
        y = ROOT.Double(0)
        ratio.GetPoint(i, x, y)
        y_errlo = ratio.GetErrorYlow(i)
        y_errhi = ratio.GetErrorYhigh(i)

        [bl,bh] = getBin(binning, x)
        if not dict[cat][etaPtVtx+subCatSuffix].has_key(str(bl)+"_"+str(bh)): dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)] = {}
        if not dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)].has_key("data/mc"): dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data/mc"] = {}
        if not ((y == 0.0) and (y_errlo == 0) and (y_errhi == 0)):
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data/mc"][etaPtVtx] = x.real
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data/mc"]["efficiency_ratio"] = y.real
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data/mc"]["err_low"] = y_errlo
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data/mc"]["err_hi"] = y_errhi
        else:
            remove_ratio.append(i)
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data/mc"][etaPtVtx] = nan
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data/mc"]["efficiency_ratio"] = nan
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data/mc"]["err_low"] = nan
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["data/mc"]["err_hi"] = nan

    remove_SF = []
    for i in range(SF.GetN()):
        x = ROOT.Double(0)
        y = ROOT.Double(0)
        SF.GetPoint(i, x, y)
        y_errlo = SF.GetErrorYlow(i)
        y_errhi = SF.GetErrorYhigh(i)

        print "x = ", x, " y = ", y

        [bl,bh] = getBin(binning, x.real)
        if not dict[cat][etaPtVtx+subCatSuffix].has_key(str(bl)+"_"+str(bh)): dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)] = {}
        if not dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)].has_key("TriggerBias_MC"): dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["TriggerBias_MC"] = {}
        if not ((y == 0.0) and (y_errlo == 0) and (y_errhi == 0)):
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["TriggerBias_MC"][etaPtVtx] = x.real
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["TriggerBias_MC"]["efficiency"] = y.real
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["TriggerBias_MC"]["err_low"] = y_errlo
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["TriggerBias_MC"]["err_hi"] = y_errhi
        else:
            remove_data.append(i)
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["TriggerBias_MC"][etaPtVtx] = nan
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["TriggerBias_MC"]["efficiency"] = nan
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["TriggerBias_MC"]["err_low"] = nan
            dict[cat][etaPtVtx+subCatSuffix][str(bl)+"_"+str(bh)]["TriggerBias_MC"]["err_hi"] = nan

    rm = 0
    for i in remove_data:
        h_data.RemovePoint(i-rm)
        rm = rm+1
    rm = 0
    for i in remove_mc:
        h_mc.RemovePoint(i-rm)
        rm = rm+1
        print "mc "+cat+etaPtVtx+subCatSuffix + "  " +str(i)
    rm = 0
    for i in remove_ratio:
        ratio.RemovePoint(i-rm)
        rm = rm+1
    rm = 0
    for i in remove_SF:
        SF.RemovePoint(i-rm)
        rm = rm+1

    rootoutput.cd()
    h_data.SetName("DATA_"+cat+"_"+etaPtVtx+"_"+subCatSuffix)
    h_data.SetTitle("DATA_"+cat+"_"+etaPtVtx+"_"+subCatSuffix)
    h_mc.SetName("MC_"+cat+"_"+etaPtVtx+"_"+subCatSuffix)
    h_mc.SetTitle("MC_"+cat+"_"+etaPtVtx+"_"+subCatSuffix)
    ratio.SetName("DATA/MC_"+cat+"_"+etaPtVtx+"_"+subCatSuffix)
    ratio.SetTitle("DATA/MC_"+cat+"_"+etaPtVtx+"_"+subCatSuffix)
    SF.SetName("TriggerBias_MC_"+cat+"_"+etaPtVtx+"_"+subCatSuffix)
    SF.SetTitle("TriggerBias_MC_"+cat+"_"+etaPtVtx+"_"+subCatSuffix)
    h_data.Write()
    h_mc.Write()
    ratio.Write()
    SF.Write()


def runAll():
    dict = {}

    rootoutput = ROOT.TFile.Open("MuonEfficiencies_Jpsi_26June2014_run2012ABCD_53X.root","Recreate")

    getFile( 'MuonID_Loose_eta_April2014.root' , dict, rootoutput, "Loose", "eta", "pt8-20"  )
    getFile( 'MuonID_Loose_pt_abseta0_April2014.root' , dict, rootoutput, "Loose", "pt", "abseta<0.9" )
    getFile( 'MuonID_Loose_pt_abseta1_June2014.root' , dict, rootoutput, "Loose", "pt", "abseta0.9-1.2" )
    getFile( 'MuonID_Loose_pt_abseta2_June2014.root' , dict, rootoutput, "Loose", "pt", "abseta1.2-2.1" )
    getFile( 'MuonID_Loose_plateau_abseta_April2014.root' , dict, rootoutput, "Loose", "abseta", "pt8-20" )
    getFile( 'MuonID_Loose_vtx_April2014.root' , dict, rootoutput, "Loose", "vtx", "pt8-20_abseta<2.1" )
    #getFile( 'MuonID_Soft_eta_April2014.root' , dict, rootoutput, "Soft", "eta", "pt8-20"  )
    #getFile( 'MuonID_Soft_pt_abseta0_April2014.root' , dict, rootoutput, "Soft", "pt", "abseta<0.9" )
    #getFile( 'MuonID_Soft_pt_abseta1_April2014.root' , dict, rootoutput, "Soft", "pt", "abseta0.9-1.2" )
    #getFile( 'MuonID_Soft_pt_abseta2_April2014.root' , dict, rootoutput, "Soft", "pt", "abseta1.2-2.1" )
    #getFile( 'MuonID_Soft_plateau_abseta_April2014.root' , dict, rootoutput, "Soft", "abseta", "pt8-20" )
    getFile( 'MuonID_newSoft_eta_April2014.root' , dict, rootoutput, "Soft", "eta", "pt8-20"  )
    getFile( 'MuonID_newSoft_pt_abseta0_April2014.root' , dict, rootoutput, "Soft", "pt", "abseta<0.9" )
    getFile( 'MuonID_newSoft_pt_abseta1_June2014.root' , dict, rootoutput, "Soft", "pt", "abseta0.9-1.2" )
    getFile( 'MuonID_newSoft_pt_abseta2_June2014.root' , dict, rootoutput, "Soft", "pt", "abseta1.2-2.1" )
    getFile( 'MuonID_newSoft_plateau_abseta_April2014.root' , dict, rootoutput, "Soft", "abseta", "pt8-20" )
    getFile( 'MuonID_newSoft_vtx_April2014.root' , dict, rootoutput, "Soft", "vtx", "pt8-20_abseta<2.1" )
    getFile( 'MuonID_Tight_eta_April2014.root' , dict, rootoutput, "Tight", "eta", "pt8-20"  )
    getFile( 'MuonID_Tight_pt_abseta0_April2014.root' , dict, rootoutput, "Tight", "pt", "abseta<0.9" )
    getFile( 'MuonID_Tight_pt_abseta1_June2014.root' , dict, rootoutput, "Tight", "pt", "abseta0.9-1.2" )
    getFile( 'MuonID_Tight_pt_abseta2_June2014.root' , dict, rootoutput, "Tight", "pt", "abseta1.2-2.1" )
    getFile( 'MuonID_Tight_plateau_abseta_April2014.root' , dict, rootoutput, "Tight", "abseta", "pt8-20" )
    getFile( 'MuonID_Tight_vtx_April2014.root' , dict, rootoutput, "Tight", "vtx", "pt8-20_abseta<2.1" )

    f = open('MuonEfficiencies_Jpsi_26June2014_run2012ABCD_53X.pkl', 'w')
    pickle.dump(dict, f)
    f.close()

    rootoutput.Close()

    return dict

runAll()
