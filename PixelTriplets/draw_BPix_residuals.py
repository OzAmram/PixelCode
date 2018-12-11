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

def getHist(tree, var, name, cut, nBins = 50, binLow = -300., binHigh = 300.):
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
    label = sys.argv[2]
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

    cut1left = cut12 + "&& (layer1SizeY % 2 == 1) && (layer1ymin %2  == 0)"
    cut1right = cut12 + "&& (layer1SizeY % 2 == 1) && (layer1ymin %2  == 1)"

    cut1even = cut12 + "&& (layer1SizeY % 2 == 0)"

    cut2left = cut12 + "&& (layer2SizeY % 2 == 1) && (layer2ymin %2  == 0)"
    cut2right = cut12 + "&& (layer2SizeY % 2 == 1) && (layer2ymin %2  == 1)"

    h1_resx = getHist(pTree, "layer1dx", "layer1x", cut12)
    h1_resz = getHist(pTree, "layer1dz", "layer1z", cut12)

    h2_resx = getHist(pTree, "layer2dx", "layer2x", cut12)
    h2_resz = getHist(pTree, "layer2dz", "layer2z", cut12)

    h3_resx = getHist(pTree, "layer3dx", "layer3x", cut34)
    h3_resz = getHist(pTree, "layer3dz", "layer3z", cut34)

    h4_resx = getHist(pTree, "layer4dx", "layer4x", cut34)
    h4_resz = getHist(pTree, "layer4dz", "layer4z", cut34)

    h1_resx.GetXaxis().SetTitle("#Delta X [#mum]")
    h1_resz.GetXaxis().SetTitle("#Delta Z [#mum]")

    h2_resx.GetXaxis().SetTitle("#Delta X [#mum]")
    h2_resz.GetXaxis().SetTitle("#Delta Z [#mum]")

    h3_resx.GetXaxis().SetTitle("#Delta X [#mum]")
    h3_resz.GetXaxis().SetTitle("#Delta Z [#mum]")

    h4_resx.GetXaxis().SetTitle("#Delta X [#mum]")
    h4_resz.GetXaxis().SetTitle("#Delta Z [#mum]")

    h1_leftz = getHist(pTree, "layer1dz", "layer1zleft", cut1left)
    h1_rightz = getHist(pTree, "layer1dz", "layer1zright", cut1right)
    h1_evenz = getHist(pTree, "layer1dz", "layer1zeven", cut1even)

    h2_leftz = getHist(pTree, "layer2dz", "layer2zleft", cut2left)
    h2_rightz = getHist(pTree, "layer2dz", "layer2zright", cut2right)

    h1_leftx = getHist(pTree, "layer1dx", "layer1xleft", cut1left)
    h1_rightx = getHist(pTree, "layer1dx", "layer1xright", cut1right)



    lstyle = 0
    lColor = kBlack
    lTag2 = 'Track p_{T}>12 GeV'
    outDir = 'plots/BPix'

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_resx,label,lColor,lstyle,"Layer 1",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_resz,label,lColor,lstyle,"Layer 1",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_resx,label,lColor,lstyle,"Layer 2",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_resz,label,lColor,lstyle,"Layer 2",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h3_resx,label,lColor,lstyle,"Layer 3",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h3_resz,label,lColor,lstyle,"Layer 3",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h4_resx,label,lColor,lstyle,"Layer 4",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h4_resz,label,lColor,lstyle,"Layer 4",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_leftz,label,lColor,lstyle,"Layer 1 left",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_rightz,label,lColor,lstyle,"Layer 1 right",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_evenz,label,lColor,lstyle,"Layer 1 even length",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_leftz,label,lColor,lstyle,"Layer 2 left",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_rightz,label,lColor,lstyle,"Layer 2 right",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_leftx,label,lColor,lstyle,"Layer 1 left",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_rightx,label,lColor,lstyle,"Layer 1 right",lTag2,outDir,0,100)



