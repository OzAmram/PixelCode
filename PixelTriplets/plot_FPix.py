import os
import glob
import math
import array
import sys
import time
from optparse import OptionParser
import ROOT
from ROOT import *


#from plotHelpers import *

lFiles_dict ={}
lFiles_dict['1D'] = [line.rstrip() for line in open('files_1D.txt')]
lFiles_dict['2D'] = [line.rstrip() for line in open('files_2D.txt')]
lFiles_dict['2Dforce'] = [line.rstrip() for line in open('files_2Dforce.txt')]

lHists = ['dx_resolution_study_1','dz_resolution_study_1',
          'dx_resolution_study_2','dz_resolution_study_2',
          'dx_resolution_study_3','dz_resolution_study_3']
lHists_dict = {}
def main(options,args):

    for iKey,lFiles in lFiles_dict.iteritems():
        for i0,iFile in enumerate(lFiles):
            if 'log' in iFile: continue
            print "root://cmsxrootd.fnal.gov/"+iFile.replace('/eos/uscms','')
            pFile = ROOT.TFile.Open("root://cmsxrootd.fnal.gov/"+iFile.replace('/eos/uscms',''))
            pTree = pFile.Get("FPixResolution_Template/tree")
            for i1,iHist in enumerate(lHists):
                sTmp = "pTmp_%s_%s"%(iKey,iHist)
                pTmp = ROOT.TH1F(sTmp,sTmp,100,-10000,10000);
                pTree.Draw("%s>>%s"%(iHist,sTmp))
                pTmp = ROOT.gDirectory.Get(sTmp)
                pTmp.SetDirectory(0)
                if iKey+'_'+iHist in lHists_dict: lHists_dict[iKey+'_'+iHist].Add(pTmp)
                else: lHists_dict[iKey+'_'+iHist] = pTmp.Clone()
                lHists_dict[iKey+'_'+iHist].SetDirectory(0)
            pFile.Close()

    pOut=ROOT.TFile.Open('resolution_fpix.root','RECREATE')
    for iKey,iHist in lHists_dict.iteritems(): iHist.Write()
    pOut.Close()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option("--lumi", dest="lumi", type=float, default = 30,help="luminosity", metavar="lumi")
    parser.add_option('-i','--idir', dest='idir', default = '2017sel8pt475409old/',help='directory with data', metavar='idir')
    parser.add_option('-o','--odir', dest='odir', default = 'plots/',help='directory to write plots', metavar='odir')
    parser.add_option('--jobs', dest='jobs', default=20, type=int, help='#of jobs to process in parallel the trees [default: %default]')
    parser.add_option('--isMuonCR', action='store_true', dest='isMuonCR', default=False, help='run on muon CR')
    parser.add_option('--isData', action='store_true', dest='isData', default=False, help='run on data')
    parser.add_option('--jet', dest='jet', default='AK8', help='jet type')
    parser.add_option('--blind', dest='blind', default=False, help='10th of dataset')
    parser.add_option('--is2016', action='store_true', dest='is2016', default=False, help='run 2016')

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
