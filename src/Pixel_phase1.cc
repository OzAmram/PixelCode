// -*- C++ -*-dz%
//
// Package:  Pixel
// Class:    Pixel
//
// my/Pixel/src/Pixel.cc
//
// Pixel (and strip) triplet residuals
//
// Author: Valere Lambert, UZH 2015
//
// Based off of original Triplet Author:  Daniel Pitzl, DESY
//         Created:  Sat Feb 12 12:12:42 CET 2011
// 
//
// system include files:
#include <memory>
#include <iostream>
#include <iomanip>
#include <cmath>
// ROOT:
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"

// CMS and user include files:
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include <DataFormats/BeamSpot/interface/BeamSpot.h>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <DataFormats/TrackReco/interface/HitPattern.h>

#include <MagneticField/Engine/interface/MagneticField.h>
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"

#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"


#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryBuilder.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include <TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h>
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "DataFormats/TrackerRecHit2D/interface/TkCloner.h"

#include "TrackingTools/DetLayers/interface/DetLayer.h"
//#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelLorentzAngle.h"
#include "CondFormats/DataRecord/interface/SiPixelLorentzAngleRcd.h"

// Helix script from Morris Schwartz
#include "SimpleHelix.h"

// Flag for new tracking rechis, has to be ON for pre7 and later   
#define NEW_TRACKINGRECHITS  // For V71X_pre7 and later 

class Pixel_phase1 : public edm::EDAnalyzer{
public:
  explicit Pixel_phase1(const edm::ParameterSet&);
  ~Pixel_phase1();

private:
  virtual void beginJob() ;
  virtual void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void getResiduals(const edm::Event&, const edm::EventSetup&, std::string detTag);
  std::vector<double> getIntersection(std::vector<double> p1, std::vector<double> p2, double rho, const GeomDet *detHit, std::vector<double> intersection);

  edm::InputTag _triggerSrc;
  std::string _ttrhBuilder;
  HLTConfigProvider HLTConfig;
  bool doBPix;
  bool doFPix;
  int _orbit_beginning;
  int _orbit_end;

  float dx_res_1;
  float dz_res_1;
  float dx_res_2;
  float dz_res_2;
  float dx_res_3;
  float dz_res_3;

  float pt_res;
  float pt_res_refit;
  int cluster_size_res_init;
  int cluster_size_res;
  int hits_track;
  int hits_barrel;
  int hits_endcap;
  int ls_with_measure;
  int ds_with_measure;

  bool isTriplet;
  TTree * tree;

  Double_t layer1dx, layer1dz, layer2dx, layer2dz, layer3dx, layer3dz;
  Double_t trkEta, trkPt;
  Int_t layer1EdgeTypeY, layer2EdgeTypeY, layer3EdgeTypeY;
  Int_t layer1EdgeTypeX, layer2EdgeTypeX, layer3EdgeTypeX;
  bool layer1HasBadPixels, layer2HasBadPixels, layer3HasBadPixels;
  bool layer1OnEdge, layer2OnEdge, layer3OnEdge;
  Int_t pxn1, pxn2, pxn3, pxn4;
  Int_t layer1xmax, layer1xmin, layer1ymin, layer1ymax, layer1SizeY, layer1SizeX;
  Int_t layer2xmax, layer2xmin, layer2ymin, layer2ymax, layer2SizeY, layer2SizeX;
  Int_t layer3xmax, layer3xmin, layer3ymin, layer3ymax, layer3SizeY, layer3SizeX;
  Float_t layer1Charge, layer2Charge, layer3Charge;


  edm::EDGetTokenT<reco::BeamSpot>  t_offlineBeamSpot_;
  edm::EDGetTokenT<reco::VertexCollection> t_offlinePrimaryVertices_ ;
  edm::EDGetTokenT <reco::TrackCollection>  t_generalTracks_;
  edm::EDGetTokenT< edm::View<reco::PFMET>> t_pfMet_;

  // ----------member data ---------------------------

  std::string processName_;
  int run_num = -999;
  int lumi_block = -999;

  //Cluster Probability
  int clusSize_X = -999;
  int clusSize_Y = -999;
  float clusProb_FPix = -999;
  
  // Rechit coordinates in local - 3 hits (xpx1_l, xpy1_l), (xpx2_l, xpy2_l), and (xpx3_l, xpy3_l)
  double xpx1_l = -999;
  double xpy1_l = -999;

  double xpx2_l = -999;
  double xpy2_l = -999;
  
  double xpx3_l = -999;
  double xpy3_l = -999;

  double xpx4_l = -999;
  double xpy4_l = -999;

  //**

  double xpx1_l_r1 = -999;
  double xpy1_l_r1 = -999;

  double xpx2_l_r1 = -999;
  double xpy2_l_r1 = -999;
  
  double xpx3_l_r1 = -999;
  double xpy3_l_r1 = -999;

  double xpx4_l_r1 = -999;
  double xpy4_l_r1 = -999;

  //**

  double xpx1_l_r2 = -999;
  double xpy1_l_r2 = -999;

  double xpx2_l_r2 = -999;
  double xpy2_l_r2 = -999;
  
  double xpx3_l_r2 = -999;
  double xpy3_l_r2 = -999;

  double xpx4_l_r2 = -999;
  double xpy4_l_r2 = -999;

  // Estimated coordinates in local 
  double xl_ideal_1 = -999;
  double yl_ideal_1 = -999;

  double xl_ideal_2 = -999;
  double yl_ideal_2 = -999;

  double xl_ideal_3 = -999;
  double yl_ideal_3 = -999;

  double xl_ideal_4 = -999;
  double yl_ideal_4 = -999;

  double xl_ideal_2_r1 = -999;
  double yl_ideal_2_r1 = -999;

  double xl_ideal_3_r1 = -999;
  double yl_ideal_3_r1 = -999;

  double xl_ideal_4_r1 = -999;
  double yl_ideal_4_r1 = -999;

  double xl_ideal_2_r2 = -999;
  double yl_ideal_2_r2 = -999;

  double xl_ideal_3_r2 = -999;
  double yl_ideal_3_r2 = -999;

  double xl_ideal_4_r2 = -999;
  double yl_ideal_4_r2 = -999;

  // Residuals
  double residual_x_1 = -999;
  double residual_x_2 = -999;
  double residual_x_3 = -999;
  double residual_x_4 = -999;

  double residual_x_2_r1 = -999;
  double residual_x_3_r1 = -999;
  double residual_x_4_r1 = -999;
  double residual_x_2_r2 = -999;
  double residual_x_3_r2 = -999;
  double residual_x_4_r2 = -999;

  double residual_y_1 = -999;
  double residual_y_2 = -999;
  double residual_y_3 = -999;
  double residual_y_4 = -999;

  double residual_y_2_r1 = -999;
  double residual_y_3_r1 = -999;
  double residual_y_4_r1 = -999;
  double residual_y_2_r2 = -999;
  double residual_y_3_r2 = -999;
  double residual_y_4_r2 = -999;

  // Errors on Rechit coordinates in local
  double x_local_error_1 = -999;
  double y_local_error_1 = -999;

  double x_local_error_2 = -999;
  double y_local_error_2 = -999;

  double x_local_error_3 = -999;
  double y_local_error_3 = -999;

  double x_local_error_4 = -999;
  double y_local_error_4 = -999;

  // Helix parameters
  int Pass = -99;
  int nloops = -99;
  int insideD = -99;
  double radius = -99;
  double xcenter = -99;
  double ycenter = -99;
  double dzdphi = -99;
  double z0 = -99;
  double rho = -99;

};

class myCountersPixel_phase1{
   public:
      static int neve;
      static unsigned int prevrun;
};

int myCountersPixel_phase1::neve = 0;
unsigned int myCountersPixel_phase1::prevrun = 0;


Pixel_phase1::Pixel_phase1(const edm::ParameterSet& iConfig)
{
  std::cout << "PxlFPix constructed\n";
  _triggerSrc = iConfig.getParameter<edm::InputTag>("triggerSource");
  _ttrhBuilder = iConfig.getParameter<std::string>("ttrhBuilder");
  doBPix=iConfig.getParameter<bool>("doBPix");
  doFPix=iConfig.getParameter<bool>("doFPix");
  _orbit_beginning=iConfig.getParameter<int>("orbit_beginning");
  _orbit_end=iConfig.getParameter<int>("orbit_end");
  //std::cout<<_triggerSrc<<" "<<_triggerSrc.label()<<" "<<_triggerSrc.process()<<" "
  //	   <<_triggerSrc.instance()<<" "<<std::endl;

  t_offlineBeamSpot_ =    consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  t_offlinePrimaryVertices_ =   consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
  t_generalTracks_= consumes<reco::TrackCollection> (edm::InputTag("generalTracks"));
  t_pfMet_= consumes< edm::View<reco::PFMET>>(edm::InputTag("pfMet"));
  
  edm::Service<TFileService> fsT;
  tree = fsT->make<TTree>("tree", "tree");
  tree->Branch("layer1dx", &layer1dx);
  tree->Branch("layer1dz", &layer1dz);
  tree->Branch("layer2dx", &layer2dx);
  tree->Branch("layer2dz", &layer2dz);
  tree->Branch("layer3dx", &layer3dx);
  tree->Branch("layer3dz", &layer3dz);
  tree->Branch("pxn1", &pxn1);
  tree->Branch("pxn2", &pxn2);
  tree->Branch("pxn3", &pxn3);
  tree->Branch("pxn4", &pxn4);
  tree->Branch("trkEta", &trkEta);
  tree->Branch("trkPt", &trkPt);
  tree->Branch("layer1OnEdge", &layer1OnEdge);
  tree->Branch("layer2OnEdge", &layer2OnEdge);
  tree->Branch("layer3OnEdge", &layer3OnEdge);
  tree->Branch("layer1HasBadPixels", &layer1HasBadPixels);
  tree->Branch("layer2HasBadPixels", &layer2HasBadPixels);
  tree->Branch("layer3HasBadPixels", &layer3HasBadPixels);

  tree->Branch("layer1SizeX", &layer1SizeX);
  tree->Branch("layer1SizeY", &layer1SizeY);
  tree->Branch("layer1xmax", &layer1xmax);
  tree->Branch("layer1xmin", &layer1xmin);
  tree->Branch("layer1ymax", &layer1ymax);
  tree->Branch("layer1ymin", &layer1ymin);

  tree->Branch("layer2SizeX", &layer2SizeX);
  tree->Branch("layer2SizeY", &layer2SizeY);
  tree->Branch("layer2xmax", &layer2xmax);
  tree->Branch("layer2xmin", &layer2xmin);
  tree->Branch("layer2ymax", &layer2ymax);
  tree->Branch("layer2ymin", &layer2ymin);

  tree->Branch("layer3SizeX", &layer3SizeX);
  tree->Branch("layer3SizeY", &layer3SizeY);
  tree->Branch("layer3xmax", &layer3xmax);
  tree->Branch("layer3xmin", &layer3xmin);
  tree->Branch("layer3ymax", &layer3ymax);
  tree->Branch("layer3ymin", &layer3ymin);

  tree->Branch("layer1Charge", &layer1Charge);
  tree->Branch("layer2Charge", &layer2Charge);
  tree->Branch("layer3Charge", &layer3Charge);

}
Pixel_phase1::~Pixel_phase1()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


// member functions:
// method called once each job just before starting event loop

void Pixel_phase1::beginJob()
{

}

//----------------------------------------------------------------------
// method called for each event:

void Pixel_phase1::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}

std::vector<double> Pixel_phase1::getIntersection(std::vector<double> point1, std::vector<double> point2, double rho,  const GeomDet *detHit, std::vector<double> intersection){
  /* 
         Takes two points and the curvature to create a helix, then finds the intersection with a detector plane
	 returns the (x,y,z) local coordinates for the intersection point
  */


  // Create helix from p1, p2 and rho
  SimpleHelix Hel = SimpleHelix(point1, point2, fabs(rho), &Pass);
  
  // Information about the detector plane
  double x_0 = detHit->position().x();
  double y_0 = detHit->position().y();
  double z_0 = detHit->position().z();

  double nX = detHit->surface().normalVector().x();
  double nY = detHit->surface().normalVector().y();
  double nZ = detHit->surface().normalVector().z();

  std::vector<double> plane = {x_0, y_0, z_0, nX, nY, nZ};

  // Find the intersection of the detector plane and helix
  nloops = Hel.SimpleHelix::pposdir(plane, intersection);

  // Convert global intersection point to local
  Surface::GlobalPoint EstimatedPoint(intersection[0],intersection[1],intersection[2]);
  Surface::LocalPoint LocalEstimate = detHit->toLocal(EstimatedPoint);

  // boolean of if the intersection is within the detector plane bounds - currently not used
  insideD = detHit->surface().bounds().inside(LocalEstimate);

  // Get Helix Parameters - currently not used
  Hel.SimpleHelix::parameters(radius, xcenter, ycenter, dzdphi, z0);
  
  std::vector<double> EstLocalCoordinates = {LocalEstimate.x(), LocalEstimate.y(), LocalEstimate.z()};

  return EstLocalCoordinates;
}



void Pixel_phase1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  bool isTriplet;

  if(doFPix && (int)iEvent.orbitNumber() >= (int)_orbit_beginning && (int)iEvent.orbitNumber()<= (int)_orbit_end ){
    std::string detTag = "fpix";
    Pixel_phase1::getResiduals(iEvent, iSetup, detTag);
  }
  if(doBPix && (int)iEvent.orbitNumber() >= (int)_orbit_beginning && (int)iEvent.orbitNumber()<= (int)_orbit_end ){
    std::string detTag = "bpix";
    Pixel_phase1::getResiduals(iEvent, iSetup, detTag);
  }

}

void Pixel_phase1::getResiduals(const edm::Event & iEvent, const edm::EventSetup& iSetup, std::string detTag){
  
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace math;
  
  myCountersPixel_phase1::neve++;
  
  if( myCountersPixel_phase1::prevrun != iEvent.run() ){
    time_t unixZeit = iEvent.time().unixTime();
    cout << "new run " << iEvent.run();
    cout << ", LumiBlock " << iEvent.luminosityBlock();
    cout << " taken " << ctime(&unixZeit); // ctime has endline
    myCountersPixel_phase1::prevrun = iEvent.run();
  }// new run
  
  int idbg = 0;
  if( myCountersPixel_phase1::neve < 2 ) idbg = 1;
  
  int jdbg = 0;
  if( idbg ) {
    cout << endl;
    cout << "run " << iEvent.run();
    cout << ", LumiBlock " << iEvent.luminosityBlock();
    cout << ", event " << iEvent.eventAuxiliary().event();
    time_t unixZeit = iEvent.time().unixTime();
    cout << ", taken " << ctime(&unixZeit); // ctime has endline
  }
  
  
  run_num=iEvent.run();
  lumi_block=iEvent.luminosityBlock();   
  
  //--------------------------------------------------------------------
  // beam spot:
  
  edm::Handle<reco::BeamSpot> rbs;
  iEvent.getByToken( t_offlineBeamSpot_, rbs );

  XYZPoint bsP = XYZPoint(0,0,0);       
  if( rbs.failedToGet() ) return;
  if( ! rbs.isValid() ) return;
  bsP = XYZPoint( rbs->x0(), rbs->y0(), rbs->z0() );

  double bx = rbs->BeamWidthX();
  double by = rbs->BeamWidthY();
 
  //Retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> tTopoH;
  iSetup.get<TrackerTopologyRcd>().get(tTopoH);
  const TrackerTopology *tTopo=tTopoH.product();

  // primary vertices:
  Handle<VertexCollection> vertices;
  iEvent.getByToken( t_offlinePrimaryVertices_,vertices );
  if( vertices.failedToGet() ) return;
  if( !vertices.isValid() ) return;
  int nvertex = vertices->size();
  XYZPoint vtxN = XYZPoint(0,0,0);
  XYZPoint vtxP = XYZPoint(0,0,0);
  double bestNdof = 0;
  double maxSumPt = 0;
  Vertex bestPvx;
  for( VertexCollection::const_iterator iVertex = vertices->begin();
       iVertex != vertices->end(); ++iVertex ) {
    if( ! iVertex->isValid() ) continue;
    else {
      if( iVertex->isFake() ) continue;
      else{
	if( iVertex->ndof() > bestNdof ) {
	  bestNdof = iVertex->ndof();
	  vtxN = XYZPoint( iVertex->x(), iVertex->y(), iVertex->z() );
	}
	if( iVertex->p4().pt() > maxSumPt ) {
	  maxSumPt = iVertex->p4().pt();
	  vtxP = XYZPoint( iVertex->x(), iVertex->y(), iVertex->z() );
	  bestPvx = *iVertex;
	}
      }
    }
  }
  if( maxSumPt < 1 ) return;

  // MET:  
  edm::Handle< edm::View<reco::PFMET> > pfMEThandle;
  iEvent.getByToken(t_pfMet_, pfMEThandle );
  
  // get a fitter to refit TrackCandidates, the same fitter as used in standard reconstruction:
  
  // Fitter                                                                                                                                                                         
  edm::ESHandle<TrajectoryFitter> aFitter;
  iSetup.get<TrajectoryFitter::Record>().get("KFFittingSmootherWithOutliersRejectionAndRK",aFitter);
  std::unique_ptr<TrajectoryFitter> theFitter = aFitter->clone();         // pointer which destroys object when pointer out of scope
  
  // Transient Rechit Builders                                                                                                                                                     
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get( "TransientTrackBuilder", theB );
  
  // Transient rec hits:                                                                                                                                                          
  ESHandle<TransientTrackingRecHitBuilder> hitBuilder;
  iSetup.get<TransientRecHitRecord>().get( _ttrhBuilder, hitBuilder );
  
  // Cloner, New from 71Xpre7                                                                                                                                                   
  const TkTransientTrackingRecHitBuilder * builder = static_cast<TkTransientTrackingRecHitBuilder const *>(hitBuilder.product());
  auto hitCloner = builder->cloner();
  theFitter->setHitCloner(&hitCloner);
  
  // Trackpropagator:
  edm::ESHandle<Propagator> prop;
  iSetup.get<TrackingComponentsRecord>().get( "PropagatorWithMaterial", prop );
  const Propagator* thePropagator = prop.product();
  
  // tracks:
  Handle<TrackCollection> tracks;
  iEvent.getByToken( t_generalTracks_, tracks );
  if( tracks.failedToGet() ) return;
  if( !tracks.isValid() ) return;

  // get tracker geometry:
  edm::ESHandle<TrackerGeometry> pTG;
  iSetup.get<TrackerDigiGeometryRecord>().get( pTG );
  if( ! pTG.isValid() ) {
    cout << "Unable to find TrackerDigiGeometry. Return\n";
    return;
  }
  

  double sumpt = 0;     // total pt of tracks from vtx
  double sumq = 0;      // total charge from vtx

  for( TrackCollection::const_iterator iTrack = tracks->begin();
       iTrack != tracks->end(); ++iTrack ) {

    isTriplet = true;
    double pt = iTrack->pt();

    trkPt = pt;
    trkEta = iTrack->eta();

    if( abs( iTrack->dxy(vtxP) ) > 5*iTrack->dxyError() ) continue; // if trans. IP > 5x its error, skip
    sumpt += pt;
    sumq += iTrack->charge();
    double logpt = log(pt) / log(10);
    const reco::HitPattern& hp = iTrack->hitPattern();    

    const double pi = 4*atan(1);
    const double wt = 180/pi;
    const double twopi = 2*pi;
    const double pihalf = 2*atan(1);
    const double sqrtpihalf = sqrt(pihalf);

    double phi = iTrack->phi();
    double dca = iTrack->d0(); // w.r.t. origin                                                                                                          
    //double dca = -iTrack->dxy(); // dxy = -d0                                                                                                          
    double dip = iTrack->lambda();
    double z0  = iTrack->dz();
    double tet = pihalf - dip;
    //double eta = iTrack->eta();                                                                                                                        

    // beam line at z of track, taking beam tilt into account                                                                                            
    double zBeam = iTrack->dz(bsP);//z0p of track along beam line w.r.t. beam z center                                                                   
    double xBeam = rbs->x0() + rbs->dxdz() * zBeam;//beam at z of track                                                                                  
    double yBeam = rbs->y0() + rbs->dydz() * zBeam;
    double z0p =  zBeam + bsP.z(); // z0p of track along beam line w.r.t. CMS z = 0                                                                      
    XYZPoint blP = XYZPoint( xBeam, yBeam, z0p );//point on beam line at z of track                                                                      

    double bcap = iTrack->dxy(blP);//impact parameter to beam                                                                                            
    double edca = iTrack->dxyError();
    double ebca = sqrt( edca*edca + bx*by );//round beam                                                                                                 
    //double sbca = bcap / ebca;//impact parameter significance                                                                                            

    if( hp.trackerLayersWithMeasurement() < 7 ) continue; // select only tracks that go into strips

    // transient track:    
    TransientTrack tTrack = theB->build(*iTrack);
    TrajectoryStateOnSurface initialTSOS = tTrack.innermostMeasurementState();
    double kap = tTrack.initialFreeState().transverseCurvature();                          // curvature of track
    rho = 1/kap;  
      
    //   Get the Pixel Hits from the track for the triplet
    // rec hits from track extra:
    if( iTrack->extra().isNull() ) continue;//next track
    if( ! iTrack->extra().isAvailable() ) continue;//next track

    uint32_t innerDetId = 0;
    double xPX1 = 0;      // global x hit 1
    double yPX1 = 0;      // global y hit 1
    double zPX1 = 0;      // global z hit 1
    double ePX1 = 0;      // sqrt(covariance of x)
    double fPX1 = 0;      // sqrt(covariance of y)
    
    double xPX2 = 0;
    double yPX2 = 0;
    double zPX2 = 0;
    double ePX2 = 0;
    double fPX2 = 0;
    
    double xPX3 = 0;
    double yPX3 = 0;
    double zPX3 = 0;
    double ePX3 = 0;
    double fPX3 = 0;

    double xPX4 = 0;
    double yPX4 = 0;
    double zPX4 = 0;
    double ePX4 = 0;
    double fPX4 = 0;

    int n1 = 0;
    int n2 = 0;
    int n3 = 0;
    int n4 = 0;
    const GeomDet * det1 = NULL;  // Detector for first hit
    const GeomDet * det2 = NULL;
    const GeomDet * det3 = NULL;
    const GeomDet * det4 = NULL;    

    double xPX2_r1 = 0;
    double yPX2_r1 = 0;
    double zPX2_r1 = 0;
    
    double xPX3_r1 = 0;
    double yPX3_r1 = 0;
    double zPX3_r1 = 0;

    double xPX4_r1 = 0;
    double yPX4_r1 = 0;
    double zPX4_r1 = 0;

    int n2_r1 = 0;
    int n3_r1 = 0;
    int n4_r1 = 0;

    const GeomDet * det2_r1 = NULL;
    const GeomDet * det3_r1 = NULL;
    const GeomDet * det4_r1 = NULL;    

    //**********
    
    double xPX2_r2 = 0;
    double yPX2_r2 = 0;
    double zPX2_r2 = 0;
    
    double xPX3_r2 = 0;
    double yPX3_r2 = 0;
    double zPX3_r2 = 0;

    double xPX4_r2 = 0;
    double yPX4_r2 = 0;
    double zPX4_r2 = 0;

    int n2_r2 = 0;
    int n3_r2 = 0;
    int n4_r2 = 0;

    const GeomDet * det2_r2 = NULL;
    const GeomDet * det3_r2 = NULL;
    const GeomDet * det4_r2 = NULL;    

    edm::OwnVector<TrackingRecHit> recHitVector;                     // for seed
    std::vector<TransientTrackingRecHit::RecHitPointer> myTTRHvec;
    Trajectory::RecHitContainer coTTRHvec;                           // for fit, constant
    
    // loop over recHits on this track:

    for( trackingRecHit_iterator irecHit = iTrack->recHitsBegin();
	 irecHit != iTrack->recHitsEnd(); ++irecHit ) {
      DetId detId = (*irecHit)->geographicalId();                          // get detector 
      uint32_t subDet = detId.subdetId();                                  // get subdetector
      recHitVector.push_back( (*irecHit)->clone() );
            
      // build transient hit: 
      auto tmprh = (*irecHit)->cloneForFit(*builder->geometry()->idToDet((**irecHit).geographicalId()));
      auto transRecHit = hitCloner.makeShared(tmprh, initialTSOS);
      myTTRHvec.push_back( transRecHit );
      coTTRHvec.push_back( transRecHit );
      
      if( ! (*irecHit)->isValid() ) continue;

      if( subDet == PixelSubdetector::PixelEndcap) {
	int idisk=tTopo->pxfDisk(detId);
	const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>( &*(*irecHit) );
	edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const & clust = pixhit->cluster();
	if( idisk == 1 ) {
	  layer1HasBadPixels = pixhit->hasBadPixels();
	  layer1OnEdge = pixhit->isOnEdge();
	  layer1Charge = clust->charge();
	  layer1SizeX = clust->sizeX();
	  layer1SizeY = clust->sizeY();
	  layer1xmin = clust->minPixelRow();
	  layer1xmax = clust->maxPixelRow();
	  layer1ymin = clust->minPixelCol();
	  layer1ymax = clust->maxPixelCol();
	}
	else if(idisk ==2){
	  layer2HasBadPixels = pixhit->hasBadPixels();
	  layer2OnEdge = pixhit->isOnEdge();
	  layer2Charge = clust->charge();
	  layer2SizeX = clust->sizeX();
	  layer2SizeY = clust->sizeY();
	  layer2xmin = clust->minPixelRow();
	  layer2xmax = clust->maxPixelRow();
	  layer2ymin = clust->minPixelCol();
	  layer2ymax = clust->maxPixelCol();
	}
	else if(idisk ==3){
	  layer3HasBadPixels = pixhit->hasBadPixels();
	  layer3OnEdge = pixhit->isOnEdge();
	  layer3Charge = clust->charge();
	  layer3SizeX = clust->sizeX();
	  layer3SizeY = clust->sizeY();
	  layer3xmin = clust->minPixelRow();
	  layer3xmax = clust->maxPixelRow();
	  layer3ymin = clust->minPixelCol();
	  layer3ymax = clust->maxPixelCol();
        }
      }

      double xloc = transRecHit->localPosition().x();       // 1st meas coord
      double yloc = transRecHit->localPosition().y();       // 2nd meas coord or zero
      double vxloc = transRecHit->localPositionError().xx();//covariance
      double vyloc = transRecHit->localPositionError().yy();//covariance
      double gX = transRecHit->globalPosition().x();
      double gY = transRecHit->globalPosition().y();
      double gZ = transRecHit->globalPosition().z();

      if( transRecHit->canImproveWithTrack() ) {//use z from track to apply alignment
	TrajectoryStateOnSurface propTSOS = thePropagator->propagate( initialTSOS, transRecHit->det()->surface() );
	if( propTSOS.isValid() ){
	  auto preciseHit = hitCloner.makeShared(tmprh,propTSOS); //pre7  

	  xloc = preciseHit->localPosition().x();// 1st meas coord
	  yloc = preciseHit->localPosition().y();// 2nd meas coord or zero

	  vxloc = preciseHit->localPositionError().xx();//covariance
	  vyloc = preciseHit->localPositionError().yy();//covariance
	  
	  gX = preciseHit->globalPosition().x();
	  gY = preciseHit->globalPosition().y();
	  gZ = preciseHit->globalPosition().z();
	}
      }//canImprove
      
      // ============================================================== 
      if(detTag == "fpix"){
	
	// PXB:      
	if( subDet == PixelSubdetector::PixelBarrel ) {

	  int ilay=tTopo->pxbLayer(detId);
	  
	  if( ilay == 2 ) {
	    
	    n1++;
	    xPX1 = gX;
	    yPX1 = gY;
	    zPX1 = gZ;
	    xpx1_l = xloc;
	    xpy1_l = yloc;	        
  	    ePX1 = sqrt( vxloc );
	    fPX1 = sqrt( vyloc );

	    det1 = transRecHit->det();

	  }// BPIX1
	  
	}//PXB1
    
	if( subDet == PixelSubdetector::PixelEndcap) {

	  int idisk=tTopo->pxfDisk(detId); 
	  int blade  = tTopo->pxfBlade(detId);     // Phase 1: Inner blades 1-22, Outer blades 23-56
	  int ring = 1 + (blade>22);               // Phase 1: Inner: 1, Outer: 2
	  
	  if( idisk == 1 ){
	    n2++;
	    xPX2 = gX; // precise hit in CMS global coordinates
	    yPX2 = gY;
	    zPX2 = gZ;
	    xpx2_l = xloc; // precise hit in local coordinates (w.r.t. sensor center)
	    xpy2_l = yloc;
	    ePX2 = sqrt( vxloc );
	    fPX2 = sqrt( vyloc );
	        	        
	    det2 = transRecHit->det();
	        
	    const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>( &*(*irecHit) );
	    
	    if( pixhit->hasFilledProb() ){
	      clusProb_FPix = pixhit->clusterProbability(0);                                      
	    }
	   
	    edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const & clust = pixhit->cluster();

	    if( clust.isNonnull() ) {
	      clusSize_Y = clust->sizeY();
	      clusSize_X = clust->sizeX();
	    }
	    
	    if( ring == 2 ){
	      n2_r2++;
	      xPX2_r2 = gX; // precise hit in CMS global coordinates
	      yPX2_r2 = gY;
	      zPX2_r2 = gZ;
	      xpx2_l_r2 = xloc; // precise hit in local coordinates (w.r.t. sensor center)
	      xpy2_l_r2 = yloc;
	      //ePX2_r2 = sqrt( vxloc );
	      //fPX2_r2 = sqrt( vyloc );
	      
	      det2_r2 = transRecHit->det();
	      
	      const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>( &*(*irecHit) );
	      edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const & clust = pixhit->cluster();
	    }

	    if( ring == 1 ){
	      n2_r1++;
	      xPX2_r1 = gX; // precise hit in CMS global coordinates
	      yPX2_r1 = gY;
	      zPX2_r1 = gZ;
	      xpx2_l_r1 = xloc; // precise hit in local coordinates (w.r.t. sensor center)
	      xpy2_l_r1 = yloc;
	      
	      det2_r1 = transRecHit->det();
	      
	      const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>( &*(*irecHit) );
	      edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const & clust = pixhit->cluster();
	    }


	  }//PXF1
	    
	  if( idisk == 2 ){
	        
	    n3++;
	    xPX3 = gX;
	    yPX3 = gY;
	    zPX3 = gZ;
	    xpx3_l = xloc;
	    xpy3_l = yloc;
	    ePX3 = sqrt( vxloc );
	    fPX3 = sqrt( vyloc );

	    det3 = transRecHit->det();


	    if( ring == 2 ){
	      n3_r2++;
	      xPX3_r2 = gX; // precise hit in CMS global coordinates
	      yPX3_r2 = gY;
	      zPX3_r2 = gZ;
	      xpx3_l_r2 = xloc; // precise hit in local coordinates (w.r.t. sensor center)
	      xpy3_l_r2 = yloc;
	      
	      det3_r2 = transRecHit->det();
	      
	      const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>( &*(*irecHit) );
	      edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const & clust = pixhit->cluster();
	    }


	    if( ring == 1 ){
	      n3_r1++;
	      xPX3_r1 = gX; // precise hit in CMS global coordinates
	      yPX3_r1 = gY;
	      zPX3_r1 = gZ;
	      xpx3_l_r1 = xloc; // precise hit in local coordinates (w.r.t. sensor center)
	      xpy3_l_r1 = yloc;
	      
	      det3_r1 = transRecHit->det();
	      
	      const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>( &*(*irecHit) );
	      edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const & clust = pixhit->cluster();
	    }


	  }//PXF2
	  	  
	  if( idisk == 3 ){
	        
	    n4++;
	    xPX4 = gX;
	    yPX4 = gY;
	    zPX4 = gZ;
	    xpx4_l = xloc;
	    xpy4_l = yloc;
	    ePX4 = sqrt( vxloc );
	    fPX4 = sqrt( vyloc );

	    det4 = transRecHit->det();

	    if( ring == 2 ){
	      n4_r2++;
	      xPX4_r2 = gX; // precise hit in CMS global coordinates
	      yPX4_r2 = gY;
	      zPX4_r2 = gZ;
	      xpx4_l_r2 = xloc; // precise hit in local coordinates (w.r.t. sensor center)
	      xpy4_l_r2 = yloc;
	      
	      det4_r2 = transRecHit->det();
	      
	      const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>( &*(*irecHit) );
	      edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const & clust = pixhit->cluster();
	    }

	    if( ring == 1 ){
	      n4_r1++;
	      xPX4_r1 = gX; // precise hit in CMS global coordinates
	      yPX4_r1 = gY;
	      zPX4_r1 = gZ;
	      xpx4_l_r1 = xloc; // precise hit in local coordinates (w.r.t. sensor center)
	      xpy4_l_r1 = yloc;

	      det4_r1 = transRecHit->det();
	      
	      const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>( &*(*irecHit) );
	      edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const & clust = pixhit->cluster();
	    }

	    }//PXF3

	}//PXF
	
      }//doFPix
      else{
	cout << "detector tag not specified"<< endl;
      }
    
    }//loop rechits
  
    //------------------------------------------------------------------------
    // refit the track:
    
    PTrajectoryStateOnDet PTraj = trajectoryStateTransform::persistentState( initialTSOS, innerDetId );
    const TrajectorySeed seed( PTraj, recHitVector, alongMomentum );
    
    //if( idbg ) cout << "  have seed\n";

    std::vector<Trajectory> refitTrajectoryCollection = theFitter->fit( seed, coTTRHvec, initialTSOS );
    
    if( refitTrajectoryCollection.size() > 0 ) { // should be either 0 or 1            
      const Trajectory& refitTrajectory = refitTrajectoryCollection.front();
      // Trajectory.measurements:
      Trajectory::DataContainer refitTMs = refitTrajectory.measurements();

      pt_res_refit = refitTrajectory.geometricalInnermostState().globalMomentum().perp();

      // trajectory residuals:
      
      for( Trajectory::DataContainer::iterator iTM = refitTMs.begin(); iTM != refitTMs.end(); iTM++ ) {
	if( ! iTM->recHit()->isValid() ) continue;
	DetId detId = iTM->recHit()->geographicalId();
	uint32_t subDet = detId.subdetId();
	
	// enum SubDetector{ PixelBarrel=1, PixelEndcap=2 };
	// enum SubDetector{ TIB=3, TID=4, TOB=5, TEC=6 };
	
	double xHit = iTM->recHit()->localPosition().x(); // primary measurement direction
	double yHit = iTM->recHit()->localPosition().y(); // always 0 in strips
	
	double dx = xHit - iTM->predictedState().localPosition().x();
	double dy = yHit - iTM->predictedState().localPosition().y();
	//double vxh = iTM->recHit()->localPositionError().xx();//covariance Not used
	//double vxt = iTM->predictedState().localError().positionError().xx();// ditto
	
	//if( subDet == 1 && idbg ){//1=PXB
	//if( subDet == 4 && idbg ){4=TID
	  //cout << "  predictdStateResid = " << dx*1E4 << " um";
	  //cout << ", eh = " << sqrt(vxh)*1E4 << " um";
	  //cout << ", et = " << sqrt(vxt)*1E4 << " um";
	  //cout << endl;
	//}

	TrajectoryStateOnSurface combinedPredictedState =
	  TrajectoryStateCombiner().combine( iTM->forwardPredictedState(), iTM->backwardPredictedState() );
	
	if( ! combinedPredictedState.isValid() ) continue;//skip hit
	
	if( jdbg ) cout << "  have combinedPredictedState\n";
	
	double xptch;
	double yptch;	
	
	if( subDet <  3 ){//1,2=pixel
	  PixelTopology & pixelTopol = (PixelTopology&) iTM->recHit()->detUnit()->topology();
	  xptch = pixelTopol.pitch().first;
	  yptch = pixelTopol.pitch().second;
	}
	else {//strip
	  StripTopology & stripTopol = (StripTopology&) iTM->recHit()->detUnit()->topology();
	  xptch = stripTopol.localPitch( combinedPredictedState.localPosition() );
	  yptch = stripTopol.localStripLength( combinedPredictedState.localPosition() );
	}
      
	dx = xHit - combinedPredictedState.localPosition().x(); //x = primary measurement
	dy = yHit - combinedPredictedState.localPosition().y(); //
	//vxh = iTM->recHit()->localPositionError().xx();//covariance Unused
	//vxt = combinedPredictedState.localError().positionError().xx();// ditto
	
	// angles of incidence:
	// local z = upwards = normal vector
	// local x = primary measurement direction
	// local y = secondary measurement direction
	
	// use Topology. no effect in PXB, essential in TID, TEC
	
	const Topology* theTopology = &(iTM->recHit()->detUnit()->topology() );
	    
	// MeasurementPoint [pitch] (like channel number)
	
	// TODO: Use the measurementPosition(point, trackdir) version of this function in order to take bows into account!
	MeasurementPoint hitMeasurement = theTopology->measurementPosition( iTM->recHit()->localPosition() );
	
	// TID and TEC have trapezoidal detectors:
	// translation from channel number into local x depends on local y
	// track prediction has local x,y => can convert into proper channel number MeasurementPoint:
	
	// TODO: Use the measurementPosition(point, trackdir) version of this function in order to take bows into account!
	MeasurementPoint combinedPredictedMeasurement = theTopology->measurementPosition( combinedPredictedState.localPosition() );

	dx = hitMeasurement.x() - combinedPredictedMeasurement.x(); //in units of pitch  (pitch = size of pixel or strip)
	dy = hitMeasurement.y() - combinedPredictedMeasurement.y(); //in units of pitch
	dx = dx * xptch;//convert back into [cm] using local pitch
	dy = dy * yptch;//[cm]
	
	// use topology: needed for TEC	
	double xx = hitMeasurement.x();
	double yy;
	if( subDet < 3 ) // pixel is 2D
	  yy = hitMeasurement.y();
	else // strips are 1D
	  yy = combinedPredictedMeasurement.y();
	
	MeasurementPoint mp( xx, yy );
	
	//2012 StripTopology & stripTopol = (StripTopology&) iTM->recHit()->detUnit()->topology();
	
	Surface::LocalPoint lp = theTopology->localPosition( mp );
	const GeomDet * myGeomDet = iTM->recHit()->det(); // makes no difference in TEC
	Surface::GlobalPoint gp = myGeomDet->toGlobal( lp );
	
	double gX = gp.x();
	double gY = gp.y();
	double gZ = gp.z();
	
	double lX = lp.x();
	double lY = lp.y();
      

	//overwrite PXB global coordinates once more, using topology:

	// =================================================
	// Modify BPIX Points with topology refit
	// =================================================
	if(detTag == "bpix"){
	  if( subDet == PixelSubdetector::PixelBarrel ) {
	        
	    int ilay=tTopo->pxbLayer(detId);
 
	    if( ilay == 1 ) {
	      xPX1 = gX;
	      yPX1 = gY;
	      zPX1 = gZ;

	      xpx1_l=lX;
	      xpy1_l=lY;
	      det1 = iTM->recHit()->det();
	    }// layer 1
	    else if( ilay == 2 ) {
	      xPX2 = gX;
	      yPX2 = gY;
	      zPX2 = gZ;
	            
	      xpx2_l=lX;
	      xpy2_l=lY;
	      det2 = iTM->recHit()->det();
	    }// layer 2
	    else if( ilay == 3 ) {
	      xPX3 = gX;
	      yPX3 = gY;
	      zPX3 = gZ;
	            
	      xpx3_l=lX;
	      xpy3_l=lY;
	      det3 = iTM->recHit()->det();
	    }// layer 3
	    else if( ilay == 4) {
	      xPX4 = gX;
	      yPX4 = gY;
	      zPX4 = gZ;
	      
	      xpx4_l=lX;
	      xpy4_l=lY;
	      det4 = iTM->recHit()->det();
	      }// layer 4
	      
	        
	  }// barrel
	  }// doBPix

	// ================================================
	// Modify FPIX Points with topology refit
	// ================================================
	else if(detTag =="fpix"){

	  if( subDet == PixelSubdetector::PixelBarrel ) {

	    int ilay=tTopo->pxbLayer(detId);

	    if( ilay == 2 ) {
	      xPX1 = gX;
	      yPX1 = gY;
	      zPX1 = gZ;

	      xpx1_l=lX;
	      xpy1_l=lY;
	      det1 = iTM->recHit()->det();
	    }// layer 1
	  }// barrel	  
	  else if( subDet == PixelSubdetector::PixelEndcap) {
	    
	    int idisk=tTopo->pxfDisk(detId);
	    int blade  = tTopo->pxfBlade(detId);     // Phase 1: Inner blades 1-22, Outer blades 23-56
	    int ring = 1 + (blade>22);               // Phase 1: Inner: 1, Outer: 2

	    if( idisk == 1 ) {
	      xPX2 = gX;
	      yPX2 = gY;
	      zPX2 = gZ;
	      
	      xpx2_l=lX;
	      xpy2_l=lY;
	      det2 = iTM->recHit()->det();
	    
	      if( ring == 2 ){
		xpx2_l_r2 = lX; // precise hit in local coordinates (w.r.t. sensor center)
		xpy2_l_r2 = lY;
	      }

	      if( ring == 1 ){
		xpx2_l_r1 = lX; // precise hit in local coordinates (w.r.t. sensor center)
		xpy2_l_r1 = lY;
	    }

	    }// disk 1
	    else if( idisk == 2 ) {
	      xPX3 = gX;
	      yPX3 = gY;
	      zPX3 = gZ;
	      
	      xpx3_l=lX;
	      xpy3_l=lY;
	      det3 = iTM->recHit()->det();


	      if( ring == 2 ){
		xPX3_r2 = gX; // precise hit in CMS global coordinates
		yPX3_r2 = gY;
		zPX3_r2 = gZ;
		xpx3_l_r2 = lX; // precise hit in local coordinates (w.r.t. sensor center)
		xpy3_l_r2 = lY;
		
		det3_r2 = iTM->recHit()->det();
	      }

	      if( ring == 1 ){
		xPX3_r1 = gX; // precise hit in CMS global coordinates
		yPX3_r1 = gY;
		zPX3_r1 = gZ;
		xpx3_l_r1 = lX; // precise hit in local coordinates (w.r.t. sensor center)
		xpy3_l_r1 = lY;
		
		det3_r1 = iTM->recHit()->det(); 
	    }


	    }// disk 2
	    else if( idisk == 3 ) {
	      xPX4 = gX;
	      yPX4 = gY;
	      zPX4 = gZ;
	      
	      xpx4_l=lX;
	      xpy4_l=lY;
	      det4 = iTM->recHit()->det();

	      if( ring == 2 ){
		xPX4_r2 = gX; // precise hit in CMS global coordinates
		yPX4_r2 = gY;
		zPX4_r2 = gZ;
		xpx4_l_r2 = lX; // precise hit in local coordinates (w.r.t. sensor center)
		xpy4_l_r2 = lY;
		
		det4_r2 = iTM->recHit()->det();
	      }

	      if( ring == 1 ){
		xPX4_r1 = gX; // precise hit in CMS global coordinates
		yPX4_r1 = gY;
		zPX4_r1 = gZ;
		xpx4_l_r1 = lX; // precise hit in local coordinates (w.r.t. sensor center)
		xpy4_l_r1 = lY;
		
		det4_r1 = iTM->recHit()->det();
	    }


	      }// disk 3
	    
	  }// endcaps
	}// doFPix
	else {
	  cout << "detTag not specified"<<endl;
	}
	
      }//loop iTM
      
    }//refitted trajectory
  
    //------------------------------------------------------------------------
    // 1-2-3 pixel triplet:
    
    pxn1 = n1;
    pxn2 = n2;
    pxn3 = n3;
    pxn4 = n4;

    if( n1*n2*n3 > 0 ) {

      {// let's open a scope, so we can redefine the variables further down
	
	// create points to feed helix
	std::vector<double> p1 = {xPX1, yPX1, zPX1};
	std::vector<double> p2 = {xPX2, yPX2, zPX2};
	std::vector<double> p3 = {xPX3, yPX3, zPX3};
	
	std::vector<double> intersection1 = {};
	std::vector<double> intersection2 = {};
	std::vector<double> intersection3 = {};
	
	// Create helix from two points and curvature, return the intersection point in local coordinates
	std::vector<double> IntersectionPointLocal_1 = Pixel_phase1::getIntersection(p2, p3, rho, det1, intersection1); 
	std::vector<double> IntersectionPointLocal_2 = Pixel_phase1::getIntersection(p1, p3, rho, det2, intersection2);
	std::vector<double> IntersectionPointLocal_3 = Pixel_phase1::getIntersection(p1, p2, rho, det3, intersection3);
	
	// Intersection point in local coordinates
	xl_ideal_1 = IntersectionPointLocal_1[0];
	yl_ideal_1 = IntersectionPointLocal_1[1];
	
	xl_ideal_2 = IntersectionPointLocal_2[0];
	yl_ideal_2 = IntersectionPointLocal_2[1];
	
	xl_ideal_3 = IntersectionPointLocal_3[0];
	yl_ideal_3 = IntersectionPointLocal_3[1];
	
	// Residuals with rechit and intersection point
	residual_x_1= (xpx1_l - xl_ideal_1)*1E4;
	residual_y_1= (xpy1_l - yl_ideal_1)*1E4;
	
	residual_x_2= (xpx2_l - xl_ideal_2)*1E4;
	residual_y_2= (xpy2_l - yl_ideal_2)*1E4;
	
	residual_x_3= (xpx3_l - xl_ideal_3)*1E4;
	residual_y_3= (xpy3_l - yl_ideal_3)*1E4;
	
	if(n2_r1>0){
	  std::vector<double> p2_r1 = {xPX2_r1, yPX2_r1, zPX2_r1};
	  std::vector<double> intersection2_r1 = {};
	  std::vector<double> IntersectionPointLocal_2_r1 = Pixel_phase1::getIntersection(p1, p3, rho, det2_r1, intersection2_r1);
          xl_ideal_2_r1 = IntersectionPointLocal_2_r1[0];
          yl_ideal_2_r1 = IntersectionPointLocal_2_r1[1];
          residual_x_2_r1= (xpx2_l_r1 - xl_ideal_2_r1)*1E4;
          residual_y_2_r1= (xpy2_l_r1 - yl_ideal_2_r1)*1E4;
        }
	
	if(n2_r2>0){
	  std::vector<double> p2_r2 = {xPX2_r2, yPX2_r2, zPX2_r2};
	  std::vector<double> intersection2_r2 = {};
	  std::vector<double> IntersectionPointLocal_2_r2 = Pixel_phase1::getIntersection(p1, p3, rho, det2_r2, intersection2_r2);
          xl_ideal_2_r2 = IntersectionPointLocal_2_r2[0];
          yl_ideal_2_r2 = IntersectionPointLocal_2_r2[1];
          residual_x_2_r2= (xpx2_l_r2 - xl_ideal_2_r2)*1E4;
          residual_y_2_r2= (xpy2_l_r2 - yl_ideal_2_r2)*1E4;
        }

	// Local errors for rechit
	x_local_error_1 = ePX1*1E4;
	y_local_error_1 = fPX1*1E4;
	
	x_local_error_2 = ePX2*1E4;
	y_local_error_2 = fPX2*1E4;
	
	x_local_error_3 = ePX3*1E4;
	y_local_error_3 = fPX3*1E4;
      
	// Fill Histograms for FPIX
	if(detTag == "fpix"){
	  
	  if( pt > 0.8) {
            isTriplet = true;
          }

  	  dx_res_1 = residual_x_2;
          dz_res_1 = residual_y_2;

	  //actual residuals!!
	  if(pt>4){
	    printf("layer1dx %f for x %f ideal: %f",residual_x_2,xpx2_l,xl_ideal_1 );
	    layer1dx = residual_x_2;
	    layer1dz = residual_y_2;
	  }
	  
	}
      

      }//triplet 
    }// three hits: 1-2-3
    
    if( n2*n3*n4 > 0 ) {
      {// let's open a scope, so we can redefine the variables further down
       
	// create points to feed helix
	std::vector<double> p2 = {xPX2, yPX2, zPX2};
	std::vector<double> p3 = {xPX3, yPX3, zPX3};
	std::vector<double> p4 = {xPX4, yPX4, zPX4};
	
	std::vector<double> intersection2 = {};
	std::vector<double> intersection3 = {};
	std::vector<double> intersection4 = {};
	
	// Create helix from two points and curvature, return the intersection point in local coordinates
	std::vector<double> IntersectionPointLocal_2 = Pixel_phase1::getIntersection(p3, p4, rho, det2, intersection2); 
	std::vector<double> IntersectionPointLocal_3 = Pixel_phase1::getIntersection(p2, p4, rho, det3, intersection3);
	std::vector<double> IntersectionPointLocal_4 = Pixel_phase1::getIntersection(p2, p3, rho, det4, intersection4);

	if(n3_r1>0){
	  std::vector<double> p3_r1 = {xPX3_r1, yPX3_r1, zPX3_r1};
	  std::vector<double> intersection3_r1 = {};	
	  std::vector<double> IntersectionPointLocal_3_r1 = Pixel_phase1::getIntersection(p2, p4, rho, det3_r1, intersection3_r1);
	  xl_ideal_3_r1 = IntersectionPointLocal_3_r1[0];
	  yl_ideal_3_r1 = IntersectionPointLocal_3_r1[1];
	  residual_x_3_r1= (xpx3_l_r1 - xl_ideal_3_r1)*1E4;
	  residual_y_3_r1= (xpy3_l_r1 - yl_ideal_3_r1)*1E4;
	}

	if(n4_r1>0){
	  std::vector<double> p4_r1 = {xPX4_r1, yPX4_r1, zPX4_r1};
	  std::vector<double> intersection4_r1 = {};
	  std::vector<double> IntersectionPointLocal_4_r1 = Pixel_phase1::getIntersection(p2, p3, rho, det4_r1, intersection4_r1);
	  xl_ideal_4_r1 = IntersectionPointLocal_4_r1[0];
	  yl_ideal_4_r1 = IntersectionPointLocal_4_r1[1];
	  residual_x_4_r1= (xpx4_l_r1 - xl_ideal_4_r1)*1E4;
	  residual_y_4_r1= (xpy4_l_r1 - yl_ideal_4_r1)*1E4;
	}


	if(n3_r2>0){
	  std::vector<double> p3_r2 = {xPX3_r2, yPX3_r2, zPX3_r2};
	  std::vector<double> intersection3_r2 = {};	
	  std::vector<double> IntersectionPointLocal_3_r2 = Pixel_phase1::getIntersection(p2, p4, rho, det3_r2, intersection3_r2);
	  xl_ideal_3_r2 = IntersectionPointLocal_3_r2[0];
	  yl_ideal_3_r2 = IntersectionPointLocal_3_r2[1];
	  residual_x_3_r2= (xpx3_l_r2 - xl_ideal_3_r2)*1E4;
	  residual_y_3_r2= (xpy3_l_r2 - yl_ideal_3_r2)*1E4;
	}

	if(n4_r2>0){
	  std::vector<double> p4_r2 = {xPX4_r2, yPX4_r2, zPX4_r2};
	  std::vector<double> intersection4_r2 = {};
	  std::vector<double> IntersectionPointLocal_4_r2 = Pixel_phase1::getIntersection(p2, p3, rho, det4_r2, intersection4_r2);
	  xl_ideal_4_r2 = IntersectionPointLocal_4_r2[0];
	  yl_ideal_4_r2 = IntersectionPointLocal_4_r2[1];
	  residual_x_4_r2= (xpx4_l_r2 - xl_ideal_4_r2)*1E4;
	  residual_y_4_r2= (xpy4_l_r2 - yl_ideal_4_r2)*1E4;
	}

	// Intersection point in local coordinates
	xl_ideal_2 = IntersectionPointLocal_2[0];
	yl_ideal_2 = IntersectionPointLocal_2[1];
	
	xl_ideal_3 = IntersectionPointLocal_3[0];
	yl_ideal_3 = IntersectionPointLocal_3[1];
	
	xl_ideal_4 = IntersectionPointLocal_4[0];
	yl_ideal_4 = IntersectionPointLocal_4[1];
	
	// Residuals with rechit and intersection point
	residual_x_2= (xpx2_l - xl_ideal_2)*1E4;
	residual_y_2= (xpy2_l - yl_ideal_2)*1E4;
	
	residual_x_3= (xpx3_l - xl_ideal_3)*1E4;
	residual_y_3= (xpy3_l - yl_ideal_3)*1E4;
	
	residual_x_4= (xpx4_l - xl_ideal_4)*1E4;
	residual_y_4= (xpy4_l - yl_ideal_4)*1E4;
	
	// Local errors for rechit
	x_local_error_2 = ePX2*1E4;
	y_local_error_2 = fPX2*1E4;
	
	x_local_error_3 = ePX3*1E4;
	y_local_error_3 = fPX3*1E4;
	
	x_local_error_4 = ePX4*1E4;
	y_local_error_4 = fPX4*1E4;
      
	// Fill Histograms for FPIX
	if(detTag == "fpix"){
	  
	  if( refitTrajectoryCollection.size() > 0 ){
            const Trajectory& refitTrajectory = refitTrajectoryCollection.front();
            pt_res_refit = refitTrajectory.geometricalInnermostState().globalMomentum().perp();
          }
          else{
            pt_res_refit = -9999.;
          }

	  if( pt > 0.8) {
            isTriplet = true;
          }

	  dx_res_2 = residual_x_3;
          dz_res_2 = residual_y_3;

	  dx_res_3 = residual_x_4;
	  dz_res_3 = residual_y_4;

	  if(pt>4){
	    layer2dx = residual_x_3;
            layer2dz = residual_y_3;

            layer3dx = residual_x_4;
            layer3dz = residual_y_4;
	  }
	}
      
	else{}

      }//triplet 
    }// three hits: 2-3-4

  }// loop over tracks
  tree->Fill();
  
}//event
//----------------------------------------------------------------------
// method called just after ending the event loop:
//
void Pixel_phase1::endJob() {
  
  std::cout << "end of job after " << myCountersPixel_phase1::neve << " events.\n";
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(Pixel_phase1);

