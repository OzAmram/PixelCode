import os
import glob
import math
import array
import sys
import time
from optparse import OptionParser
import ROOT
from ROOT import *

from plotUtils import *

lFiles_dict ={}
lFiles_dict['1D'] = [line.rstrip() for line in open('files_1D.txt')]
lFiles_dict['2D'] = [line.rstrip() for line in open('files_2D.txt')]
lFiles_dict['2Dforce'] = [line.rstrip() for line in open('files_2Dforce.txt')]

lHists = ['dx_resolution_study_1','dz_resolution_study_1',
          'dx_resolution_study_2','dz_resolution_study_2',
          'dx_resolution_study_3','dz_resolution_study_3']
lHists_dict = {}

def main(options,args):

    if options.fill:
        for iKey,lFiles in lFiles_dict.iteritems():
            for i0,iFile in enumerate(lFiles):
                if i0 > 30: continue
                if 'log' in iFile: continue
                print "root://cmsxrootd.fnal.gov/"+iFile.replace('/eos/uscms','')
                pFile = ROOT.TFile.Open("root://cmsxrootd.fnal.gov/"+iFile.replace('/eos/uscms',''))
                pTree = pFile.Get("FPixResolution_Template/tree")
                for i1,iHist in enumerate(lHists):
                    sTmp = "pTmp_%s_%s"%(iKey,iHist)
                    if 'dx_resolution' in iHist:
                        pTmp = ROOT.TH1F(sTmp,sTmp,100,-150,150);
                    else:
                        pTmp = ROOT.TH1F(sTmp,sTmp,100,-300,300);
                    pTree.Draw("%s>>%s"%(iHist,sTmp),"pt_resolution_study>4")
                    pTmp = ROOT.gDirectory.Get(sTmp)
                    pTmp.SetDirectory(0)
                    if iKey+'_'+iHist in lHists_dict: lHists_dict[iKey+'_'+iHist].Add(pTmp)
                    else: lHists_dict[iKey+'_'+iHist] = pTmp.Clone()
                    lHists_dict[iKey+'_'+iHist].SetDirectory(0)
                pFile.Close()

        pOut=ROOT.TFile.Open('resolution_fpix_2.root','RECREATE')
        for iKey,iHist in lHists_dict.iteritems(): iHist.Write()
        pOut.Close()

    if options.plot:
        lLegend = {'2D': '2D Template',
                   '1D': '1D Template',
                   '2Dforce':' 2D always on'}
        lColor = {'2D': ROOT.kGreen,
                  '1D': ROOT.kBlue,
                  '2Dforce': ROOT.kRed,
                  }
        lStyle = {'2D': 1,
                  '1D': 1,
                  '2Dforce': 2,}
        iFile = ROOT.TFile.Open('resolution_fpix_2.root','read')
        iOFile = ROOT.TFile.Open('resolution_fpix_2_plots.root','recreate')
        for iHist in lHists:
            lSig = {}
            for iKey,lFiles in lFiles_dict.iteritems():
                lTag = ''
                sTmp = "pTmp_%s_%s"%(iKey,iHist)
                lSig[iKey] = iFile.Get(sTmp).Clone()
                if 'dx_resolution' in iHist:
                    lSig[iKey].GetXaxis().SetRangeUser(-150,150)
                    lSig[iKey].GetYaxis().SetTitle('Number of hits / 3. #mum')
                    lSig[iKey].GetXaxis().SetTitle('Residuals x direction (#mum)')
                else:
                    lSig[iKey].GetXaxis().SetRangeUser(-300,300)
                    lSig[iKey].GetYaxis().SetTitle('Number of hits / 6. #mum')
                    lSig[iKey].GetXaxis().SetTitle('Residuals z direction (#mum)')
                lSig[iKey].GetXaxis().SetTitleSize(20)
                lSig[iKey].GetXaxis().SetTitleOffset(1)
                lSig[iKey].GetYaxis().SetTitleSize(0.05)
                lSig[iKey].GetYaxis().SetTitleOffset(1)
                lSig[iKey].SetDirectory(0)
                if 'study_1' in iHist: lTag = 'FPix Disk 1'
                if 'study_2' in iHist: lTag = 'FPix Disk 2'
                if 'study_3' in iHist: lTag = 'FPix Disk 3'
                lTag2 = 'Track p_{T}>4 GeV'
                lmean,lmeanerr,lsigma,lsigmaerr = make1D(lSig[iKey],lLegend[iKey],lColor,lStyle,lTag,lTag2,'plots',0,100)
            #c = makeCanvasComparison(lSig,lLegend,lColor,lStyle,iHist,'plots',options.run,iOFile,False)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option("--run", dest="run", type=float, default = 3243138,help="run number", metavar="run")
    parser.add_option('--plot', action='store_true', dest='plot', default=False, help='plot histograms')
    parser.add_option('--fill', action='store_true', dest='fill', default=False, help='fill histograms')

    (options, args) = parser.parse_args()


    import tdrstyle
    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.10)
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetPaintTextFormat("1.1f")
    ROOT.gStyle.SetOptFit(0000)
    ROOT.gROOT.SetBatch()

    main(options,args)
