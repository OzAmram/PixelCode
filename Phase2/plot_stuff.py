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

def getHist(tree, var, name, cut, nBins = 30, binLow = -200., binHigh = 200.):
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

    pTree = fin.Get("ReadLocalMeasurement/PixelNtupleOnTrack")
    pTree.Print()
    outDir = 'plots_feb25'
    fpix_cut = "(subid == 2)"
    bpix_cut = "(subid == 1)"

    
    fpix_resx = getHist(pTree, "1e4*(xloc-hx)", "fpix_dx", fpix_cut)
    fpix_resy = getHist(pTree, "1e4*(yloc-hy)", "fpix_dy", fpix_cut)

    bpix_resx = getHist(pTree, "1e4*(xloc-hx)", "bpix_dx", bpix_cut)
    bpix_resy = getHist(pTree, "1e4*(y-hy)", "bpix_dy", bpix_cut)

    fpix_pullx = getHist(pTree, "(xloc-hx)/TMath::Sqrt(xxloc)", "fpix_pullx", fpix_cut, binLow=-5, binHigh=5)
    fpix_pully = getHist(pTree, "(yloc-hy)/TMath::Sqrt(yyloc)", "fpix_pully", fpix_cut, binLow=-5, binHigh=5)

    bpix_pullx = getHist(pTree, "(xloc-hx)/TMath::Sqrt(xxloc)", "bpix_pullx", bpix_cut, binLow=-5, binHigh=5)
    bpix_pully = getHist(pTree, "(yloc-hy)/TMath::Sqrt(yyloc)", "bpix_pully", bpix_cut, binLow=-5, binHigh=5)

    fpix_errx = getHist(pTree, "1e4 * TMath::Sqrt(xxloc)", "fpix_errx", fpix_cut, binLow=0., binHigh = 100.)
    fpix_erry = getHist(pTree, "1e4 * TMath::Sqrt(yyloc)", "fpix_erry", fpix_cut, binLow=0., binHigh = 100.)

    bpix_errx = getHist(pTree, "1e4 * TMath::Sqrt(xxloc)", "bpix_errx", fpix_cut, binLow=0., binHigh = 100.)
    bpix_erry = getHist(pTree, "1e4 * TMath::Sqrt(yyloc)", "bpix_erry", fpix_cut, binLow=0., binHigh = 100.)

    fpix_resx.GetXaxis().SetTitle("#Delta X [#mum]")
    fpix_resy.GetXaxis().SetTitle("#Delta Y [#mum]")
    bpix_resx.GetXaxis().SetTitle("#Delta X [#mum]")
    bpix_resy.GetXaxis().SetTitle("#Delta Y [#mum]")

    fpix_pullx.GetXaxis().SetTitle("Pull X ")
    fpix_pully.GetXaxis().SetTitle("Pull Y ")
    bpix_pullx.GetXaxis().SetTitle("Pull X ")
    bpix_pully.GetXaxis().SetTitle("Pull Y ")

    fpix_errx.GetXaxis().SetTitle("Error X [#mum]")
    fpix_erry.GetXaxis().SetTitle("Error Y [#mum]")
    bpix_errx.GetXaxis().SetTitle("Error X [#mum]")
    bpix_erry.GetXaxis().SetTitle("Error Y [#mum]")

    lstyle = 0
    lColor = kBlack
    lTag2 = ''

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(fpix_resx,label,lColor,lstyle,"FPix",lTag2,outDir)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(fpix_resy,label,lColor,lstyle,"FPix",lTag2,outDir)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(bpix_resx,label,lColor,lstyle,"BPix",lTag2,outDir)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(bpix_resy,label,lColor,lstyle,"BPix",lTag2,outDir)


    lmean,lmeanerr,lsigma,lsigmaerr = make1D(fpix_pullx,label,lColor,lstyle,"FPix",lTag2,outDir)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(fpix_pully,label,lColor,lstyle,"FPix",lTag2,outDir)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(bpix_pullx,label,lColor,lstyle,"BPix",lTag2,outDir)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(bpix_pully,label,lColor,lstyle,"BPix",lTag2,outDir)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(fpix_errx,label,lColor,lstyle,"FPix",lTag2,outDir, doFit=False, draw_opt = "hist")
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(fpix_erry,label,lColor,lstyle,"FPix",lTag2,outDir, doFit=False, draw_opt = "hist")
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(bpix_errx,label,lColor,lstyle,"BPix",lTag2,outDir, doFit=False, draw_opt = "hist")
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(bpix_erry,label,lColor,lstyle,"BPix",lTag2,outDir, doFit=False, draw_opt = "hist")
