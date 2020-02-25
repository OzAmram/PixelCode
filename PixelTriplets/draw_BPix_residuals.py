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
    outDir = 'plots/mar20_lowpt'

    lTag2 = 'Track p_{T}>0.75 GeV'
    cut12 = "((pxn1*pxn2*pxn3) > 0 ) && (trkPt > 0.75)"
    cut34 = "((pxn2*pxn3*pxn4) > 0 ) && (trkPt > 0.75)"
    
    cut1 = cut12 #+ " && layer1HasBadPixels"
    cut2 = cut12 #+ " && layer2HasBadPixels"
    cut3 = cut34 #+ " && layer3HasBadPixels"
    cut4 = cut34 #+ " && layer4HasBadPixels"


    cut1left = cut1 + "&& (layer1SizeY % 2 == 1) && (layer1ymin %2 ==0)"
    cut2left = cut2 + "&& (layer2SizeY % 2 == 1) && (layer2ymin %2 ==0)"
    cut3left = cut3 + "&& (layer3SizeY % 2 == 1) && (layer3ymin %2 ==0)"
    cut4left = cut4 + "&& (layer4SizeY % 2 == 1) && (layer4ymin %2 ==0)"

    cut1right = cut1 + "&& (layer1SizeY % 2 == 1) && (layer1ymin %2 ==1)"
    cut2right = cut2 + "&& (layer2SizeY % 2 == 1) && (layer2ymin %2 ==1)"
    cut3right = cut3 + "&& (layer3SizeY % 2 == 1) && (layer3ymin %2 ==1)"
    cut4right = cut4 + "&& (layer4SizeY % 2 == 1) && (layer4ymin %2 ==1)"

    cut1even = cut1 + "&& (layer1SizeY % 2 == 0)"
    cut2even = cut2 + "&& (layer2SizeY % 2 == 0)"
    cut3even = cut3 + "&& (layer3SizeY % 2 == 0)"
    cut4even = cut4 + "&& (layer4SizeY % 2 == 0)"


    h1_resx = getHist(pTree, "layer1dx", "layer1x", cut1)
    h1_resz = getHist(pTree, "layer1dz", "layer1z", cut1)

    h2_resx = getHist(pTree, "layer2dx", "layer2x", cut2)
    h2_resz = getHist(pTree, "layer2dz", "layer2z", cut2)

    h3_resx = getHist(pTree, "layer3dx", "layer3x", cut3)
    h3_resz = getHist(pTree, "layer3dz", "layer3z", cut3)

    h4_resx = getHist(pTree, "layer4dx", "layer4x", cut4)
    h4_resz = getHist(pTree, "layer4dz", "layer4z", cut4)

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
    h2_evenz = getHist(pTree, "layer2dz", "layer2zeven", cut2even)

    h3_leftz = getHist(pTree, "layer3dz", "layer3zleft", cut3left)
    h3_rightz = getHist(pTree, "layer3dz", "layer3zright", cut3right)
    h3_evenz = getHist(pTree, "layer3dz", "layer3zeven", cut3even)

    h4_leftz = getHist(pTree, "layer4dz", "layer4zleft", cut4left)
    h4_rightz = getHist(pTree, "layer4dz", "layer4zright", cut4right)
    h4_evenz = getHist(pTree, "layer4dz", "layer4zeven", cut4even)



    lstyle = 0
    lColor = kBlack

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_resx,label,lColor,lstyle,"Layer 1",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_resz,label,lColor,lstyle,"Layer 1",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_resx,label,lColor,lstyle,"Layer 2",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_resz,label,lColor,lstyle,"Layer 2",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h3_resx,label,lColor,lstyle,"Layer 3",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h3_resz,label,lColor,lstyle,"Layer 3",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h4_resx,label,lColor,lstyle,"Layer 4",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h4_resz,label,lColor,lstyle,"Layer 4",lTag2,outDir,0,100)



    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_leftz,label,lColor,lstyle,"Layer 1 left ",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_rightz,label,lColor,lstyle,"Layer 1 right ",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_leftz,label,lColor,lstyle,"Layer  2 left ",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_rightz,label,lColor,lstyle,"Layer 2 right ",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h3_leftz,label,lColor,lstyle,"Layer 3 left",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h3_rightz,label,lColor,lstyle,"Layer 3 right ",lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h4_leftz,label,lColor,lstyle,"Layer 4 left ",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h4_rightz,label,lColor,lstyle,"Layer 4 right ",lTag2,outDir,0,100)




    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_evenz,label,lColor,lstyle,"Layer 1 even ",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_evenz,label,lColor,lstyle,"Layer 2 even ",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h3_evenz,label,lColor,lstyle,"Layer 3 even ",lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h4_evenz,label,lColor,lstyle,"Layer 4 even ",lTag2,outDir,0,100)


