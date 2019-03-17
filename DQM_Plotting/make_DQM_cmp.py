import ROOT
from ROOT import *

from optparse import OptionParser
import sys

def make_ratio_plot(h1, h2, label1, label2, outfile, rebin=0):
    c1 = TCanvas("c1", "")
    p1 = TPad("pad1", "", 0., 0.3, 0.98,1.)
    p1.SetBottomMargin(0)
    p1.Draw()
    p1.cd()

    x_start = 0.3
    y_start = 0.77
    tag5 = ROOT.TLatex(x_start,y_start + 0.09,"Default Mean: %.3e +/- %.3e "%(h1.GetMean(),h1.GetMeanError()))
    tag5.SetNDC(); tag5.SetTextFont(42); tag5.SetTextSize(0.02); tag5.SetTextColor(1);
    tag6 = ROOT.TLatex(x_start,y_start + 0.06,"Default RMS: %.2e +/- %.2e "%(h1.GetRMS(),h1.GetRMSError()))
    tag6.SetNDC(); tag6.SetTextFont(42); tag6.SetTextSize(0.02); tag6.SetTextColor(1);

    tag5_CR = ROOT.TLatex(x_start,y_start + 0.03,"CR Mean: %.3e +/- %.3e "%(h2.GetMean(),h2.GetMeanError()))
    tag5_CR.SetNDC(); tag5_CR.SetTextFont(42); tag5_CR.SetTextSize(0.02); tag5_CR.SetTextColor(1);
    tag6_CR = ROOT.TLatex(x_start,y_start,"CR RMS: %.2e +/- %.2e "%(h2.GetRMS(),h2.GetRMSError()))
    tag6_CR.SetNDC(); tag6_CR.SetTextFont(42); tag6_CR.SetTextSize(0.02); tag6_CR.SetTextColor(1);

    if(rebin !=0):
        h1.Rebin(rebin)
        h2.Rebin(rebin)
    h1.SetLineColor(kBlue)
    h2.SetLineColor(kGreen)

    h1.SetMarkerSize(0)
    h2.SetMarkerSize(0)
    #rp = TRatioPlot(h1,h2)
    #rp.Draw("hist")
    h1.Draw("hist e")
    h2.Draw("hist e same")



    tag5.Draw()
    tag6.Draw()
    tag5_CR.Draw()
    tag6_CR.Draw()

    leg_size = 0.2
    leg_x = 0.53
    leg_y = 0.7
    leg = TLegend(leg_x, leg_y, leg_x + leg_size, leg_y + leg_size)
    leg.AddEntry(h1, label1)
    leg.AddEntry(h2, label2)
    leg.Draw()

    c1.cd()
    p2 = TPad("pad2", "", 0., 0., 0.98, 0.3)
    p2.SetTopMargin(0.2)
    p2.SetGridy()
    p2.Draw()
    p2.cd()
    h_ratio = h1.Clone("h_ratio")
    h_ratio.Sumw2()
    h_ratio.SetStats(0)
    h_ratio.Divide(h2)
    h_ratio.Draw("ep")


    c1.Print(outfile)



if __name__ == "__main__":

    inputFile1 = "step4_CRON.root"
    inputFile2 = "step4_CROFF.root"
    f_CR = TFile.Open(inputFile1)
    f_def = TFile.Open(inputFile2)

    if(not f_def):
        print("Cant open file %s. Exitting \n" % inputFile1)
        exit(1)

    if(not f_CR):
        print("Cant open file %s. Exitting \n" % inputFile2)
        exit(1)

    outDir = 'plots/'

    label_def = "CMSSW_10_5_pre2"
    label_CR = "CMSSW_10_5_0_pre2 + CR (bugfix)"

    #import tdrstyle
    #tdrstyle.setTDRStyle()
    gStyle.SetOptStat(0)
    #ROOT.gStyle.SetPadTopMargin(0.10)
    #ROOT.gStyle.SetPadLeftMargin(0.16)
    #ROOT.gStyle.SetPadRightMargin(0.10)
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetPaintTextFormat("1.1f")
    ROOT.gStyle.SetOptFit(0000)
    ROOT.gROOT.SetBatch()

    base_dir = "DQMData/Run 321833/"
    #my_hist = "Tracking/Run summary/TrackParameters/highPurityTracks/dzPV0p1/GeneralProperties/"
    my_hist = "PixelPhase1/Run summary/Tracks/PXBarrel/"
    my_dir = base_dir + my_hist
    name = "residual_y_PXLayer_4"
    rebin = 5
    
    f_def.cd()
    gDirectory.cd(my_dir)
    gDirectory.ls()
    h_def = gDirectory.Get(name).Clone(name + "def")
    h_def.SetDirectory(0)

    f_CR.cd()
    gDirectory.cd(my_dir)
    gDirectory.pwd()
    gDirectory.ls()
    h_CR = gDirectory.Get(name).Clone(name + "CR")
    h_CR.SetDirectory(0)


    if(not(h_def) or not(h_CR)):
        print("couldn't get hists! \n")
        exit(0)

    h_def.Print()
    h_CR.Print()

    make_ratio_plot(h_CR, h_def, label_CR, label_def, outDir + name + ".png", rebin)


