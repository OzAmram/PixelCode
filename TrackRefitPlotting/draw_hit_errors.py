
from ROOT import *

from optparse import OptionParser
from plotUtils import *
import sys

def outputResid(h, outputName):
    c = TCanvas("c", "", 0,0, 800, 800)
    c.cd()
    h.Draw("hist")
    c.Print(outputName)
    return True

def getHist(tree, var, name, cut, nBins = 100, binLow = -300., binHigh = 300.):
    h = TH1F(name, name, nBins, binLow, binHigh)
    tree.Draw("%s>>%s" %(var, name),cut)
    h = gDirectory.Get(name)
    h.SetDirectory(0)
    return h

def getOnTrackHist(h, pTree, etamin = 0., etamax = 3.0):
    max_resid = 2000.
    for e in xrange(pTree.GetEntries()):
        pTree.GetEntry(e)
        #13x21 array of floats
        resid1Dy = pTree.resid1Dy
        resid1Dx = pTree.resid1Dx

        cut = (not pTree.isOnEdge2D) and pTree.assocTrack and abs(pTree.trkEta) > etamin and abs(pTree.trkEta) < etamax
        if((abs(resid1Dy) < max_resid) and (abs(resid1Dx) < max_resid) and cut):
                h.Fill(abs(pTree.resid1Dy), pTree.onTrack)
    return h

def getOnTrackHist_v2(h, pTree, etamin, etamax):
    max_resid = 500.
    for e in xrange(pTree.GetEntries()):
        pTree.GetEntry(e)
        #13x21 array of floats
        resid1Dy = pTree.resid1Dy
        resid1Dx = pTree.resid1Dx

        cut = (not pTree.isOnEdge2D) and pTree.assocTrack and abs(pTree.trkEta) > etamin and abs(pTree.trkEta) < etamax
        if((abs(resid1Dy) < max_resid) and (abs(resid1Dx) < max_resid) and cut):
                h.Fill(sqrt(pTree.resid1Dy**2 + pTree.resid1Dx**2), pTree.onTrack)
    return h


if __name__ == "__main__":
    if(len(sys.argv) < 3):
        print("Requires input filename and output label \n")
        exit(1)

    inputFile = sys.argv[1]
    label = sys.argv[2]
    layer = 1
    if(len(sys.argv) > 3):
        layer = int(sys.argv[3])
    fin = TFile.Open(inputFile)

    outDir = 'plots/feb28'

    import tdrstyle
    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.10)
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetPaintTextFormat("1.1f")
    ROOT.gStyle.SetOptFit(0000)
    ROOT.gROOT.SetBatch()

    if(not fin):
        print("Cant open file %s. Exitting \n" % inputFile)
        exit(1)

    pTree = fin.Get("Layer1_Residuals/tree")

    #cut = "1"
    cut = "!isOnEdge2D && hasBadPixels2D && numLayers > 6"
    #cut = "(hit_type == 1)"
    #cut = "(hasBadPixels2D)"
    #cut = "(hit_type ==1  && !hasBadPixels2D && !isOnEdge2D)"
    #cut = "(isOnEdge2D)"

    cutleft = cut + "&& (clustSizeY % 2 == 1) && (clustymin %2  == 0)"
    cutright = cut + "&& (clustSizeY % 2 == 1) && (clustymin %2  == 1)"

    cuteven = cut + "&& (clustSizeY % 2 == 0)"

    nbins = 50

    h1_err1Dx = getHist(pTree, "err1Dx", "err1Dx", cut,40, 0., 100.);
    h1_err2Dx = getHist(pTree, "err2Dx", "err1Dx", cut, 40, 0., 100.);

    h1_err1Dy = getHist(pTree, "err1Dy", "err1Dy", cut,40, 0., 100.);
    h1_err2Dy = getHist(pTree, "err2Dy", "err1Dy", cut, 40, 0., 100.);

    h1_err1Dx.GetXaxis().SetTitle("Hit Error X [#mum]");
    h1_err2Dx.GetXaxis().SetTitle("Hit Error X [#mum]");

    h1_err1Dy.GetXaxis().SetTitle("Hit Error Y [#mum]");
    h1_err2Dy.GetXaxis().SetTitle("Hit Error Y [#mum]");


    outDir = 'plots/feb28'

    lstyle = 0
    lColor = kBlack
    lTag2 = 'Track p_{T}>10 GeV'
    etamin = 0.0
    etamax = 3.0

    output1D(h1_err1Dx,label,kBlue,lstyle,"1D",lTag2,outDir)
    output1D(h1_err2Dx,label,kBlue,lstyle,"CR",lTag2,outDir)

    output1D(h1_err1Dy,label,kBlue,lstyle,"1D",lTag2,outDir)
    output1D(h1_err2Dy,label,kBlue,lstyle,"CR",lTag2,outDir)



