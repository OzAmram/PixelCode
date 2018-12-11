
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

def getOnTrackHist(h, pTree):
    max_resid = 500.
    for e in xrange(pTree.GetEntries()):
        pTree.GetEntry(e)
        #13x21 array of floats
        resid1Dy = pTree.resid1Dy
        resid1Dx = pTree.resid1Dx

        cut = (pTree.isOnEdge2D) and pTree.assocTrack and (pTree.clustymin==0 or pTree.clustymax==415)
        if((abs(resid1Dy) < max_resid) and (abs(resid1Dx) < max_resid) and cut):
                h.Fill(abs(pTree.resid1Dy), pTree.onTrack)
    return h

if __name__ == "__main__":
    if(len(sys.argv) < 3):
        print("Requires input filename and output label \n")
        exit(1)

    inputFile = sys.argv[1]
    label = sys.argv[2]+"_edge"
    if(len(sys.argv) > 3):
        layer = int(sys.argv[3])
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

    pTree = fin.Get("Layer1_Residuals/tree")
    #print_edges(pTree)
    #exit(1)

    cut = "(clustxmin == 0 || clustxmax == 159 || clustymin==0 || clustymax==415) && assocTrack "
    cutx = cut + "&& (clustxmin ==0 || clustxmax ==159)"
    cuty = cut + "&& (clustymin ==0 || clustymax ==415)"

    #max x row is 159, max y col is 415
    cutxleft = cut + "&& (clustxmin  == 0)"
    cutxright = cut + "&& (clustxmax == 159)"

    cutyleft = cut + "&& (clustymin  == 0)"
    cutyright = cut + "&& (clustymax == 415)"

    cuteven = cut + "&& (clustSizeY % 2 == 0)"

    nbins = 60
    range_ = 500.
    xrange_ = 300.

    h1_resy = getHist(pTree, "resid1Dy", "1Dy", cuty, nbins, -range_,range_)
    h2_resy = getHist(pTree, "resid2Dy", "CRy", cuty, nbins, -range_, range_)


    h1_lefty = getHist(pTree, "resid1Dy", "1Dyleft", cutyleft, nbins, -range_, range_ )
    h2_lefty = getHist(pTree, "resid2Dy", "2Dyleft", cutyleft,  nbins, -range_, range_)

    h1_righty = getHist(pTree, "resid1Dy", "1Dyright", cutyright,  nbins, -range_, range_)
    h2_righty = getHist(pTree, "resid2Dy", "2Dyright", cutyright,  nbins, -range_, range_)

    h1_resx = getHist(pTree, "resid1Dx", "1Dx", cutx, nbins, -xrange_, xrange_)
    h2_resx = getHist(pTree, "resid2Dx", "CRx", cutx, nbins, -xrange_, xrange_)

    h1_leftx = getHist(pTree, "resid1Dx", "1Dxleft", cutxleft, nbins, -xrange_, xrange_ )
    h2_leftx = getHist(pTree, "resid2Dx", "2Dxleft", cutxleft,  nbins, -xrange_, xrange_)

    h1_yleftx = getHist(pTree, "resid1Dx", "1Dx_yleft", cutyleft, nbins, -xrange_, xrange_ )
    h2_yleftx = getHist(pTree, "resid2Dx", "2Dx_yleft", cutyleft,  nbins, -xrange_, xrange_)


    h1_rightx = getHist(pTree, "resid1Dx", "1Dyxright", cutxright,  nbins, -xrange_, xrange_)
    h2_rightx = getHist(pTree, "resid2Dx", "2Dyxright", cutxright,  nbins, -xrange_, xrange_)

    h1_yrightx = getHist(pTree, "resid1Dx", "1Dyx_yright", cutyright,  nbins, -xrange_, xrange_)
    h2_yrightx = getHist(pTree, "resid2Dx", "2Dyx_yright", cutyright,  nbins, -xrange_, xrange_)


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

    h2_rightx.GetXaxis().SetTitle("#Delta X [#mum]")
    h2_leftx.GetXaxis().SetTitle("#Delta X [#mum]")
    h1_leftx.GetXaxis().SetTitle("#Delta X [#mum]")
    h1_rightx.GetXaxis().SetTitle("#Delta X [#mum]")

    h2_yrightx.GetXaxis().SetTitle("#Delta X [#mum]")
    h2_yleftx.GetXaxis().SetTitle("#Delta X [#mum]")
    h1_yleftx.GetXaxis().SetTitle("#Delta X [#mum]")
    h1_yrightx.GetXaxis().SetTitle("#Delta X [#mum]")


    lstyle = 0
    lColor = kBlack
    lTag2 = 'Track p_{T}>10 GeV'
    outDir = 'plots/Dec6'

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_resx,label,lColor,lstyle,"1D: Layer %i"%layer,lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_resy,label,lColor,lstyle,"1D: Layer %i"%layer,lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_resx,label,lColor,lstyle,"CR: Layer %i"%layer,lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_resy,label,lColor,lstyle,"CR: Layer %i"%layer,lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_lefty,label,lColor,lstyle,"1D: Layer %i EdgeY ymin=0"%layer,lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_righty,label,lColor,lstyle,"1D: Layer %i EdgeY ymax=415"%layer,lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_lefty,label,lColor,lstyle,"CR: Layer %i EdgeY ymin=0"%layer,lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_righty,label,lColor,lstyle,"CR: Layer %i EdgeY ymax=415"%layer,lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_leftx,label,lColor,lstyle,"1D: Layer %i EdgeX xmin=0"%layer,lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_rightx,label,lColor,lstyle,"1D: Layer %i EdgeX xmax=159"%layer,lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_leftx,label,lColor,lstyle,"CR: Layer %i EdgeX xmin=0"%layer,lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_rightx,label,lColor,lstyle,"CR: Layer %i EdgeX xmax=159"%layer,lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_yleftx,label,lColor,lstyle,"1D: Layer %i EdgeY ymin=0"%layer,lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h1_yrightx,label,lColor,lstyle,"1D: Layer %i EdgeY ymax=415"%layer,lTag2,outDir,0,100)

    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_yleftx,label,lColor,lstyle,"CR: Layer %i EdgeY ymin=0"%layer,lTag2,outDir,0,100)
    lmean,lmeanerr,lsigma,lsigmaerr = make1D(h2_yrightx,label,lColor,lstyle,"CR: Layer %i EdgeY ymax=415"%layer,lTag2,outDir,0,100)


    h_residy_ontrack = TH2F("residy_ontrack","", 20, 0, 500, 2,0,1.2)
    getOnTrackHist(h_residy_ontrack, pTree)
    prof_y_trk = h_residy_ontrack.ProfileX("pfx_trk")
    c = TCanvas("c", "", 0,0, 1200, 800)
    c.cd()
    prof_y_trk.Draw("")
    prof_y_trk.GetYaxis().SetTitle("Fraction onTrack")
    prof_y_trk.GetXaxis().SetTitle("#DeltaY (1D) (#mu m)")
    tag = ROOT.TLatex(0.80,0.82, "Layer %i Edge" %layer)
    tag.SetNDC(); tag.SetTextFont(42); tag.SetTextSize(0.025);
    tag.Draw()
    c.SaveAs("%s/Residuals_%s_%s.png"%(outDir,label, "residy_ontrack"))


    output1D(h2_probQ,label,kBlue,lstyle,"CR",lTag2,outDir)
    output1D(h1_probQ,label,kBlue,lstyle,"1D",lTag2,outDir)

    output1D(h2_probXY,label,kBlue,lstyle,"CR",lTag2,outDir)
    output1D(h1_probXY,label,kBlue,lstyle,"1D",lTag2,outDir)


