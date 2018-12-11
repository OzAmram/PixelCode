
import ROOT
from ROOT import *

from optparse import OptionParser
from plotUtils import *
from array import array
import sys

def outputResid(h, outputName):
    c = TCanvas("c", "", 0,0, 800, 800)
    c.cd()
    h.Draw("hist")
    c.Print(outputName)
    return True



if __name__ == "__main__":
    if(len(sys.argv) != 3):
        print("Requires input filename and output label \n")
        exit(1)

    inputFile = sys.argv[1]
    label = sys.argv[2]
    fin = TFile.Open(inputFile)


    import tdrstyle
    tdrstyle.setTDRStyle()
    #ROOT.gStyle.SetPadTopMargin(0.10)
    #ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.20)
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetPaintTextFormat("1.1f")
    #ROOT.gStyle.SetOptFit(0000)
    ROOT.gROOT.SetBatch()

    if(not fin):
        print("Cant open file %s. Exitting \n" % inputFile)
        exit(1)

    pTree = fin.Get("Layer1_Residuals/tree")

    numPix = 0
    maxPix = 100
    max_row = 13
    max_col = 21
    buff = 2


    pix_height = 100
    pix_width = 150

    max_resid = 300
    outDir = 'plots/Dec4/pix/'
    clustMatrix = array('f', [0.]*max_row *max_col)
    pTree.SetBranchAddress("clustMatrix", clustMatrix)

    
    for e in xrange(pTree.GetEntries()):
        pTree.GetEntry(e)
        #13x21 array of floats
        pTree.clustMatrix
        resid1Dy = pTree.resid1Dy
        resid1Dx = pTree.resid1Dx

        if((abs(resid1Dy) < max_resid) and (abs(resid1Dx) < max_resid) and pTree.hasBadPixels2D):
            #print clustMatrix
            print ("clust size is %i x %i \n" % (pTree.clustSizeX, pTree.clustSizeY))
            if(pTree.clustSizeX > max_row or pTree.clustSizeY > max_col):
                print ("Clust is too big! \n")
                continue
            hist_rows = pTree.clustSizeX + 2*buff
            hist_cols = pTree.clustSizeY + 2*buff
            h_pix = TH2D("h" + str(numPix), "", hist_cols, -buff*pix_width, (hist_cols-buff) * pix_width, hist_rows, -buff*pix_height, (hist_rows-buff)*pix_height)
            for i in range(max_row):
                for j in range(max_col):
                    h_pix.SetBinContent( j+buff+1, i+buff+1, clustMatrix[i*max_col + j])
                    

            #local X becomes y-axis so have to swap coords
            offsetY = pTree.clust_start_y
            offsetX = pTree.clust_start_x
            tys = array('f', [pTree.trackIntersectx - offsetX])
            txs = array('f', [pTree.trackIntersecty - offsetY])
            tyerr = array('f', [ pTree.trackIntersectErrorx])
            txerr = array('f', [pTree.trackIntersectErrory])
            print("Trk intsect is %.0f +/- %.0f %.0f +/- %0.0f \n" %(pTree.trackIntersectx, pTree.trackIntersectErrorx, pTree.trackIntersecty, pTree.trackIntersectErrory))
            print("Clust start is %.0f %.0f \n" %(pTree.clust_start_x, pTree.clust_start_y))
            print("Resid 1D is %.0f %.0f \n" %(resid1Dx, resid1Dy))
            g_trk = TGraphErrors(1,  txs, tys, txerr, tyerr)
            

            ys1d = array('f', [pTree.trackIntersectx + pTree.resid1Dx - offsetX])
            xs1d = array('f', [pTree.trackIntersecty + pTree.resid1Dy - offsetY])
            ys1derr = array('f', [pTree.err1Dx])
            xs1derr = array('f', [pTree.err1Dy])
            g_1d = TGraphErrors(1, xs1d, ys1d, xs1derr, ys1derr)
            g_1d.SetMarkerColor(kGray)
            g_1d.SetMarkerSize(2)
            g_1d.SetMarkerStyle(kFullSquare)

            ys2d = array('f', [pTree.trackIntersectx + pTree.resid2Dx - offsetX])
            xs2d = array('f', [pTree.trackIntersecty + pTree.resid2Dy - offsetY])
            ys2derr = array('f', [pTree.err2Dx])
            xs2derr = array('f', [pTree.err2Dy])
            g_2d = TGraphErrors(1, xs2d, ys2d, xs2derr, ys2derr)
            g_2d.SetMarkerColor(kOrange+3)
            g_2d.SetMarkerStyle(kFullSquare)
            g_2d.SetMarkerSize(2)

            leg = TLegend(0.2,0.2)
            leg.AddEntry(g_trk, "Track Intersect", "p")
            leg.AddEntry(g_1d, "Template Reco 1D", "p")
            leg.AddEntry(g_2d, "Template Reco 2D", "p")
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            

            h_pix.GetXaxis().SetTitle("Local Y (#mu m)")
            h_pix.GetYaxis().SetTitle("Local X (#mu m)")
            c = TCanvas("c", "", 0,0, 800, 800)
            c.cd()
            h_pix.Draw("colz")
            g_trk.Draw("P same")
            g_1d.Draw("P same")
            g_2d.Draw("P same")
            leg.Draw()
            c.Print(outDir + label + '_pix' + str(numPix) + '.png')
            numPix+=1
            if(numPix > maxPix): break



