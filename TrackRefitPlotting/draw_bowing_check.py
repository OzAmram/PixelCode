

import ROOT
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

def print_edges(tree, num=50):
    count = 0
    for e in xrange(tree.GetEntries()):
        tree.GetEntry(e)

        if(tree.isOnEdge2D):
            print("X low high = %i %i. Y low high = %i %i \n"% (tree.clustxmin,tree.clustxmax, tree.clustymin, tree.clustymax))
            count+=1
            if(count >num): return None


if __name__ == "__main__":
    if(len(sys.argv) < 3):
        print("Requires input filename and output label \n")
        exit(1)

    inputFile = sys.argv[1]
    label = sys.argv[2]
    if(len(sys.argv) > 3):
        layer = int(sys.argv[3])
    fin = TFile.Open(inputFile)

    outDir = 'plots/'

    import tdrstyle
    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.10)
    ROOT.gStyle.SetPadBottomMargin(0.16)
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetPaintTextFormat("1.1f")
    ROOT.gStyle.SetOptFit(0000)
    ROOT.gROOT.SetBatch()

    if(not fin):
        print("Cant open file %s. Exitting \n" % inputFile)
        exit(1)

    pTree = fin.Get("Layer1_Residuals/tree")

    max_resid = 500
    outDir = 'plots/Dec4/'
    #print_edges(pTree)
    #exit(1)
    y_min = 35000
    h_residy_trk = TH2F("h_residy_trk","", 50, -max_resid, max_resid, 50, -y_min, y_min)
    h_residy_clust = TH2F("h_residy_clust","", 50, -max_resid, max_resid, 50, -y_min, y_min)


    for e in xrange(pTree.GetEntries()):
        pTree.GetEntry(e)
        #13x21 array of floats
        resid1Dy = pTree.resid1Dy
        resid1Dx = pTree.resid1Dx

        if((abs(resid1Dy) < max_resid) and (abs(resid1Dx) < max_resid) and pTree.hasBadPixels2D):
            if(pTree.cotbeta >0):
                h_residy_trk.Fill(resid1Dy, pTree.trackIntersecty)
                h_residy_clust.Fill(resid1Dy, pTree.clust_start_y)
            else:
                h_residy_trk.Fill(-resid1Dy, pTree.trackIntersecty)
                h_residy_clust.Fill(-resid1Dy, pTree.clust_start_y)
    prof_y_trk = h_residy_trk.ProfileY("pfx_trk")
    c = TCanvas("c", "", 0,0, 1200, 800)
    c.cd()
    prof_y_trk.Draw("")
    prof_y_trk.GetYaxis().SetTitle("#Delta Y * sign(cotbeta) (#mu m)")
    prof_y_trk.GetXaxis().SetTitle("Track Local Y (#mu m)")
    c.Print(outDir + label + 'residy_cotb_trk.png')

    prof_y_clust = h_residy_clust.ProfileY("pfx_clust")
    c2 = TCanvas("c2", "", 0,0, 1200, 800)
    c2.cd()
    prof_y_clust.Draw("")
    prof_y_clust.GetYaxis().SetTitle("#Delta Y * sign(cotbeta) (#mu m)")
    prof_y_clust.GetXaxis().SetTitle("Cluster Local Y (#mu m)")
    c2.Print(outDir + label + 'residy_cotb_clust.png')

