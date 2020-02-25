import ROOT
from ROOT import TFile, TTree, TChain, gPad, gDirectory, TVirtualFitter, kBlack
import math
import sys
from math import sqrt
import time
import array

# student-t                                                                                                                                  
isMiniBiasScan = False;
rej = 6.;
justPlot = True;
isDrawConstituents = False;
numberOfFits = 5;

def studentT(x, p):
    nn = 0;
    nn+= 1;
    dx = 0.1; 
    b1  = 0;
    if nn == 1: b1 = x[0];
    if nn == 2: dx = x[0] - b1;
    xm = p[0]
    t = (x[0] - xm ) /p[1];
    tt = t*t

    xn = 0.5 * (p[2] + 1.0);
    pk = 0.
    if( p[2] > 0.0 and p[1] > 0.0):
        #print("%.2f %.2f \n", p[2], p[1])
        aa = dx / p[1] / sqrt(p[2]*math.pi) * math.gamma(xn) / math.gamma(0.5*p[2]);
        pk = p[3] * aa * math.exp( -xn * math.log( 1.0 + tt/p[2] ) );
    return pk + p[4];


def make1D(iTmp,iLegend,iColor,iStyle,iName,iName2,iOdir, doFit = True, draw_opt = "ep same"):
    #print iTmp.GetName()
    print "Integral of tree %s is %.0f " % (iTmp.GetName(), iTmp.Integral())
    pLeg = ROOT.TLegend(0.65,0.83,0.85,0.88)
    pLeg.SetFillStyle(0)
    pLeg.SetBorderSize(0)
    pLeg.SetTextSize(0.025)
    pLeg.SetTextFont(42)

    iTmp.GetXaxis().SetTitleSize(0.04)
    lmean = lmeanErr = lwidth = lwidthErr = -1

    if(doFit):
        pFunc = ROOT.TF1("student_t", studentT, iTmp.GetBinCenter(1), iTmp.GetBinCenter(iTmp.GetNbinsX()), 5 );
        pFunc.SetParameter( 0, 0.);  pFunc.SetParName( 0, "mean"); # peak position
        print("RMS IS %.2f \n" % iTmp.GetRMS());
        #pFunc.SetParameter( 1, iTmp.GetRMS()); pFunc.SetParName( 1, "sigma"); # width
        pFunc.SetParameter( 1, 25.); pFunc.SetParName( 1, "sigma"); # width
        pFunc.SetParLimits(1, 0.001, 200.)
        pFunc.SetParameter( 2, 2.2); pFunc.SetParName( 2, "nu");# nu
        pFunc.SetParLimits(2, 0.001, 200.)
        pFunc.SetParameter( 3, iTmp.Integral()); pFunc.SetParName( 3, "area"); # N
        pFunc.SetParameter( 4, 0.); pFunc.SetParName( 4, "bkg");
        iTmp.Fit(pFunc, "RN" ,"ep");
        pFunc.SetNpx(500);
        pFunc.SetLineColor(ROOT.kGreen);
        pChi2 = pFunc.GetChisquare()/(iTmp.GetNbinsX()-5);
        pFunc.SetLineWidth(3);
        lmean = pFunc.GetParameter(0)
        lmeanErr = pFunc.GetParError(0)
        lwidth = pFunc.GetParameter(1)
        lwidthErr = pFunc.GetParError(1)
    iTmp.SetMarkerStyle(21);
    iTmp.SetMarkerSize(0.8);
    iTmp.SetMarkerColor(kBlack);
    iTmp.SetLineColor(kBlack);
    iTmp.SetTitle("");

    pLeg.AddEntry(iTmp,iLegend,"p");

    c =ROOT.TCanvas("cfit1d%s"%(iTmp.GetName()),"",800,800)
    ROOT.TGaxis.SetMaxDigits(3);
    iTmp.Draw(draw_opt);
    if(doFit): pFunc.Draw("ep same");
    pLeg.Draw()
    if(doFit):
        tag1 = ROOT.TLatex(0.70,0.80,"#mu_{r}: %.2f +/- %.2f "%(pFunc.GetParameter(0),pFunc.GetParError(0)))
        tag1.SetNDC(); tag1.SetTextFont(42); tag1.SetTextSize(0.02); tag1.SetTextColor(4);
        tag2 = ROOT.TLatex(0.70,0.77,"#sigma_{r}: %.2f +/- %.2f "%(pFunc.GetParameter(1),pFunc.GetParError(1)))
        tag2.SetNDC(); tag2.SetTextFont(42); tag2.SetTextSize(0.02); tag2.SetTextColor(4);
        tag1.Draw()
        tag2.Draw()
    tag3 = ROOT.TLatex(0.65,0.92,iName)
    tag3.SetNDC(); tag3.SetTextFont(42)
    tag4 = ROOT.TLatex(0.20,0.82,iName2)
    tag4.SetNDC(); tag4.SetTextFont(42); tag4.SetTextSize(0.025);
    tag5 = ROOT.TLatex(0.70,0.73,"Mean: %.2f +/- %.2f "%(iTmp.GetMean(),iTmp.GetMeanError()))
    tag5.SetNDC(); tag5.SetTextFont(42); tag5.SetTextSize(0.02); tag5.SetTextColor(4);
    tag6 = ROOT.TLatex(0.70,0.70,"RMS: %.2f +/- %.2f "%(iTmp.GetRMS(),iTmp.GetRMSError()))
    tag6.SetNDC(); tag6.SetTextFont(42); tag6.SetTextSize(0.02); tag6.SetTextColor(4);
    tag3.Draw()
    tag4.Draw()
    tag5.Draw()
    tag6.Draw()
    c.SaveAs("%s/Residuals_%s_%s.png"%(iOdir,iLegend, iTmp.GetName()))

    return lmean,lmeanErr,lwidth,lwidthErr
