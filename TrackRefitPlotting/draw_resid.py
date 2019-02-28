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

    outDir = 'plots/jan21'

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
    h1_resx = getHist(pTree, "resid1Dx", "1Dx", cut, nbins, -300., 300.)
    h1_resy = getHist(pTree, "resid1Dy", "1Dy", cut, nbins, -300.,300.)

    h2_resx = getHist(pTree, "resid2Dx", "CRx", cut, nbins, -300., 300.)
    h2_resy = getHist(pTree, "resid2Dy", "CRy", cut, nbins, -300., 300.)


    h1_lefty = getHist(pTree, "resid1Dy", "1Dyleft", cutleft, nbins, -300., 300. )
    h2_lefty = getHist(pTree, "resid2Dy", "2Dyleft", cutleft,  nbins, -300., 300.)

    h1_righty = getHist(pTree, "resid1Dy", "1Dyright", cutright,  nbins, -300., 300.)
    h2_righty = getHist(pTree, "resid2Dy", "2Dyright", cutright,  nbins, -300., 300.)

    h1_eveny = getHist(pTree, "resid1Dy", "1Dyeven", cuteven,  nbins, -300., 300.)
    h2_eveny = getHist(pTree, "resid2Dy", "2Dyeven", cuteven,  nbins, -300., 300.)

    h1_probQ = getHist(pTree, "probQ1D", "probQ1D", cut,40, 0., 1.1);
    h2_probQ = getHist(pTree, "probQ2D", "probQCR", cut, 40, 0., 1.1);

    h1_probXY = getHist(pTree, "probXY1D", "probXY1D", cut,40, 0., 1.1);
    h2_probXY = getHist(pTree, "probXY2D", "probXYCR", cut, 40, 0., 1.1);

    h1_probQ.GetXaxis().SetTitle("probQ");
    h2_probQ.GetXaxis().SetTitle("probQ");

    h1_probXY.GetXaxis().SetTitle("probXY");
    h2_probXY.GetXaxis().SetTitle("probXY");

    h1_resx.GetXaxis().SetTitle("#Delta X [#mum]")
    h1_resy.GetXaxis().SetTitle("#Delta Y [#mum]")

    h2_resx.GetXaxis().SetTitle("#Delta X [#mum]")
    h2_resy.GetXaxis().SetTitle("#Delta Y [#mum]")

    h2_righty.GetXaxis().SetTitle("#Delta Y [#mum]")
    h2_lefty.GetXaxis().SetTitle("#Delta Y [#mum]")
    h1_lefty.GetXaxis().SetTitle("#Delta Y [#mum]")
    h1_righty.GetXaxis().SetTitle("#Delta Y [#mum]")


    lstyle = 0
    lColor = kBlack
    lTag2 = 'Track p_{T}>10 GeV'
    etamin = 0.0
    etamax = 3.0

    h_residy_ontrack = TH2F("residy_ontrack","", 50, 0, 500, 10000,0,2)
    getOnTrackHist(h_residy_ontrack, pTree, etamin, etamax)
    prof_y_trk = h_residy_ontrack.ProfileX("pfx_trk")
    c = TCanvas("c", "", 0,0, 1200, 800)
    c.cd()
    prof_y_trk.Draw("PE")
    prof_y_trk.GetYaxis().SetTitle("Fraction onTrack")
    prof_y_trk.GetXaxis().SetTitle("#DeltaY (1D) (#mu m)")
    tag = ROOT.TLatex(0.80,0.82, "Layer %i" %layer)
    tag.SetNDC(); tag.SetTextFont(42); tag.SetTextSize(0.025);
    tag.Draw()
    c.SaveAs("%s/Residuals_%s_%s.png"%(outDir,label, "residy_ontrack"))

    h_resid_ontrack = TH2F("resid_ontrack","", 50, 0, 500, 10000,0,2)
    getOnTrackHist_v2(h_resid_ontrack, pTree, etamin, etamax)
    prof_trk = h_resid_ontrack.ProfileX("pfx_trk")
    c = TCanvas("c", "", 0,0, 1200, 800)
    c.cd()
    prof_trk.Draw("PE")
    prof_trk.GetYaxis().SetTitle("Fraction onTrack")
    prof_trk.GetXaxis().SetTitle("#Delta (1D) (#mu m)")
    tag = ROOT.TLatex(0.80,0.82, "Layer %i" %layer)
    tag.SetNDC(); tag.SetTextFont(42); tag.SetTextSize(0.025);
    tag.Draw()
    c.SaveAs("%s/Residuals_%s_%s.png"%(outDir,label, "resid_ontrack"))

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_resx,label,lColor,lstyle,"1D: Layer %i"%layer,lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_resy,label,lColor,lstyle,"1D: Layer %i"%layer,lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_resx,label,lColor,lstyle,"CR: Layer %i"%layer,lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_resy,label,lColor,lstyle,"CR: Layer %i"%layer,lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_lefty,label,lColor,lstyle,"1D: Layer %i Odd Type 1"%layer,lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_righty,label,lColor,lstyle,"1D: Layer %i Odd Type 2"%layer,lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_lefty,label,lColor,lstyle,"CR: Layer %i Odd Type 1"%layer,lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_righty,label,lColor,lstyle,"CR: Layer %i Odd Type 2"%layer,lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_eveny,label,lColor,lstyle,"1D: Layer %i Even"%layer,lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_eveny,label,lColor,lstyle,"CR: Layer %i Even"%layer,lTag2,outDir,0,100)



    output1D(h2_probQ,label,kBlue,lstyle,"CR",lTag2,outDir)
    output1D(h1_probQ,label,kBlue,lstyle,"1D",lTag2,outDir)

    output1D(h2_probXY,label,kBlue,lstyle,"CR",lTag2,outDir)
    output1D(h1_probXY,label,kBlue,lstyle,"1D",lTag2,outDir)


