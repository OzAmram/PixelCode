import ROOT
from ROOT import *

from optparse import OptionParser
from plotUtils import *
import sys
import os

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
    #pTree.Print()
    i = 0
    
    outDir = 'T25_100events'
    fpix_cut = "(subid == 2)" #&& trkPt > 5."
    bpix_cut = "(subid == 1) && layer == 1" #&& trkPt > 5."
    #bpix_cut = "(subid == 1) && (layer == 2 || layer == 3 || layer == 4)" #&& trkPt > 5."

    os.system("mkdir %s" % outDir)

    
    fpix_resx = getHist(pTree, "1e4*(x-hx)", "fpix_dx", fpix_cut, binLow = -200, binHigh = 200)
    fpix_resy = getHist(pTree, "1e4*(y-hy)", "fpix_dy", fpix_cut, binLow = -200, binHigh = 200)

    bpix_resx = getHist(pTree, "1e4*(x-hx)", "bpix_dx", bpix_cut, binLow = -200, binHigh = 200)
    bpix_resy = getHist(pTree, "1e4*(y-hy)", "bpix_dy", bpix_cut, binLow = -200, binHigh = 200)

    fpix_pullx = getHist(pTree, "(x-hx)/TMath::Sqrt(xx)", "fpix_pullx", fpix_cut, binLow=-5, binHigh=5)
    fpix_pully = getHist(pTree, "(y-hy)/TMath::Sqrt(yy)", "fpix_pully", fpix_cut, binLow=-5, binHigh=5)

    bpix_pullx = getHist(pTree, "(x-hx)/TMath::Sqrt(xx)", "bpix_pullx", bpix_cut, binLow=-5, binHigh=5)
    bpix_pully = getHist(pTree, "(y-hy)/TMath::Sqrt(yy)", "bpix_pully", bpix_cut, binLow=-5, binHigh=5)

    fpix_errx = getHist(pTree, "1e4 * TMath::Sqrt(xx)", "fpix_errx", fpix_cut, binLow=0., binHigh = 100.)
    fpix_erry = getHist(pTree, "1e4 * TMath::Sqrt(yy)", "fpix_erry", fpix_cut, binLow=0., binHigh = 100.)

    bpix_errx = getHist(pTree, "1e4 * TMath::Sqrt(xx)", "bpix_errx", fpix_cut, binLow=0., binHigh = 100.)
    bpix_erry = getHist(pTree, "1e4 * TMath::Sqrt(yy)", "bpix_erry", fpix_cut, binLow=0., binHigh = 100.)

    #fpix_tcota = getHist(pTree, "trkcota", "fpix_track_cota", fpix_cut, binLow = -3., binHigh = 3.)
    #fpix_tcotb = getHist(pTree, "trkcotb", "fpix_track_cotb", fpix_cut, binLow = -3., binHigh = 3.)

    #bpix_tcota = getHist(pTree, "trkcota", "bpix_track_cota", bpix_cut, binLow = -3., binHigh = 3.)
    #bpix_tcotb = getHist(pTree, "trkcotb", "bpix_track_cotab", bpix_cut, binLow = -3., binHigh = 3.)

    fpix_dcota = getHist(pTree, "cotAlphaFromDet", "fpix_det_cota", fpix_cut, binLow = -3., binHigh = 3.)
    fpix_dcotb = getHist(pTree, "cotBetaFromDet", "fpix_det_cotb", fpix_cut, binLow = -3., binHigh = 3.)

    bpix_dcota = getHist(pTree, "cotAlphaFromDet", "bpix_det_cota", bpix_cut, binLow = -3., binHigh = 3.)
    bpix_dcotb = getHist(pTree, "cotBetaFromDet", "bpix_det_cotb", bpix_cut, binLow = -3., binHigh = 3.)


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

    #fpix_tcota.GetXaxis().SetTitle("Track Cot(alpha)")
    #fpix_tcotb.GetXaxis().SetTitle("Track Cot(beta)")
    fpix_dcota.GetXaxis().SetTitle("Det Cot(alpha)")
    fpix_dcotb.GetXaxis().SetTitle("Det Cot(beta)")

    #bpix_tcota.GetXaxis().SetTitle("Track Cot(alpha)")
    #bpix_tcotb.GetXaxis().SetTitle("Track Cot(beta)")

    bpix_dcota.GetXaxis().SetTitle("Det Cot(alpha)")
    bpix_dcotb.GetXaxis().SetTitle("Det Cot(beta)")

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

    #lmean,lmeanerr,lsigma,lsigmaerr = make1D(fpix_tcota,label,lColor,lstyle,"FPix",lTag2,outDir, doFit=False, draw_opt = "hist")
    #lmean,lmeanerr,lsigma,lsigmaerr = make1D(fpix_tcotb,label,lColor,lstyle,"FPix",lTag2,outDir, doFit=False, draw_opt = "hist")
    #lmean,lmeanerr,lsigma,lsigmaerr = make1D(bpix_tcota,label,lColor,lstyle,"BPix",lTag2,outDir, doFit=False, draw_opt = "hist")
    #lmean,lmeanerr,lsigma,lsigmaerr = make1D(bpix_tcotb,label,lColor,lstyle,"BPix",lTag2,outDir, doFit=False, draw_opt = "hist")

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(fpix_dcota,label,lColor,lstyle,"FPix",lTag2,outDir, doFit=False, draw_opt = "hist")
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(fpix_dcotb,label,lColor,lstyle,"FPix",lTag2,outDir, doFit=False, draw_opt = "hist")
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(bpix_dcota,label,lColor,lstyle,"BPix",lTag2,outDir, doFit=False, draw_opt = "hist")
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(bpix_dcotb,label,lColor,lstyle,"BPix",lTag2,outDir, doFit=False, draw_opt = "hist")


