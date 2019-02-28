
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

def getHist(tree, var, name, cut, nBins = 30, binLow = -300., binHigh = 300.):
    h = TH1F(name, name, nBins, binLow, binHigh)
    tree.Draw("%s>>%s" %(var, name),cut)
    h = gDirectory.Get(name)
    h.SetDirectory(0)
    return h


if __name__ == "__main__":
    if(len(sys.argv) != 3):
        print("Requires input filename and output label \n")
        exit(1)

    inputFile = sys.argv[1]
    label = sys.argv[2]+"_edge"
    fin = TFile.Open(inputFile)


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

    pTree = fin.Get("BPixResolution_Template/tree")

    cut12 = "((pxn1*pxn2*pxn3) > 0 ) && (trkPt > 12)"
    cut34 = "((pxn2*pxn3*pxn4) > 0 ) && (trkPt > 12)"
    
    cut1 = cut12 + " && layer1HasBadPixels"
    cut2 = cut12 + " && layer2HasBadPixels"
    cut3 = cut34 + " && layer3HasBadPixels"
    cut4 = cut34 + " && layer4HasBadPixels"

    cut1_edge = cut12 + " && ((layer1ymin ==0) || (layer1ymax == 415))"
    cut2_edge = cut12 + " && ((layer2ymin ==0) || (layer2ymax == 415))"
    cut3_edge = cut34 + " && ((layer3ymin ==0) || (layer3ymax == 415))"
    cut4_edge = cut34 + " && ((layer4ymin ==0) || (layer4ymax == 415))"


    cut1leftedge = cut12 + "&&  (layer1ymin  == 0)"
    cut1rightedge = cut12 + "&&  (layer1ymax  == 415)"
    cut2leftedge = cut12 + "&&  (layer2ymin  == 0)"
    cut2rightedge = cut12 + "&&  (layer2ymax  == 415)"
    cut3leftedge = cut34 + "&&  (layer3ymin  == 0)"
    cut3rightedge = cut34 + "&&  (layer3ymax  == 415)"
    cut4leftedge = cut34 + "&&  (layer4ymin  == 0)"
    cut4rightedge = cut34 + "&&  (layer4ymax  == 415)"




    h1_resz_edge = getHist(pTree, "layer1dz", "layer1z_edge", cut1_edge)

    h2_resz_edge = getHist(pTree, "layer2dz", "layer2z_edge", cut2_edge)

    h3_resz_edge = getHist(pTree, "layer3dz", "layer3z_edge", cut3_edge)

    h4_resz_edge = getHist(pTree, "layer4dz", "layer4z_edge", cut4_edge)

    h1_resz_edge.GetXaxis().SetTitle("#Delta Z [#mum]")
    h2_resz_edge.GetXaxis().SetTitle("#Delta Z [#mum]")
    h3_resz_edge.GetXaxis().SetTitle("#Delta Z [#mum]")
    h4_resz_edge.GetXaxis().SetTitle("#Delta Z [#mum]")

    h1_leftz_edge = getHist(pTree, "layer1dz", "layer1zleftedge", cut1leftedge)
    h1_rightz_edge = getHist(pTree, "layer1dz", "layer1zrightedge", cut1rightedge)

    h2_leftz_edge = getHist(pTree, "layer2dz", "layer2zleftedge", cut2leftedge)
    h2_rightz_edge = getHist(pTree, "layer2dz", "layer2zrightedge", cut2rightedge)

    h3_leftz_edge = getHist(pTree, "layer3dz", "layer3zleftedge", cut3leftedge)
    h3_rightz_edge = getHist(pTree, "layer3dz", "layer3zrightedge", cut3rightedge)

    h4_leftz_edge = getHist(pTree, "layer4dz", "layer4zleftedge", cut4leftedge)
    h4_rightz_edge = getHist(pTree, "layer4dz", "layer4zrightedge", cut4rightedge)

    h1_leftz_edge.GetXaxis().SetTitle("#Delta Z [#mum]")
    h2_leftz_edge.GetXaxis().SetTitle("#Delta Z [#mum]")
    h3_leftz_edge.GetXaxis().SetTitle("#Delta Z [#mum]")
    h4_leftz_edge.GetXaxis().SetTitle("#Delta Z [#mum]")

    h1_rightz_edge.GetXaxis().SetTitle("#Delta Z [#mum]")
    h2_rightz_edge.GetXaxis().SetTitle("#Delta Z [#mum]")
    h3_rightz_edge.GetXaxis().SetTitle("#Delta Z [#mum]")
    h4_rightz_edge.GetXaxis().SetTitle("#Delta Z [#mum]")



    lstyle = 0
    lColor = kBlack
    lTag2 = 'Track p_{T}>12 GeV'
    outDir = 'plots/Jan30Edge'



    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_resz_edge,label,lColor,lstyle,"Layer 1 EdgeY",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_resz_edge,label,lColor,lstyle,"Layer 2 EdgeY",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h3_resz_edge,label,lColor,lstyle,"Layer 3 EdgeY",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h4_resz_edge,label,lColor,lstyle,"Layer 4 EdgeY",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_leftz_edge,label,lColor,lstyle,"Layer 1 left EdgeY",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_rightz_edge,label,lColor,lstyle,"Layer 1 right EdgeY",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_leftz_edge,label,lColor,lstyle,"Layer  2 left EdgeY",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_rightz_edge,label,lColor,lstyle,"Layer 2 right EdgeY",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h3_leftz_edge,label,lColor,lstyle,"Layer 3 left EdgeY",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h3_rightz_edge,label,lColor,lstyle,"Layer 3 right EdgeY",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h4_leftz_edge,label,lColor,lstyle,"Layer 4 left EdgeY",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h4_rightz_edge,label,lColor,lstyle,"Layer 4 right EdgeY",lTag2,outDir,0,100)





