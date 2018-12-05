//

// system include files:
#include <memory>
#include <iostream>
#include <iomanip>
#include <cmath>

// ROOT:
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TProfile.h"

// CMS and user include files:
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include "FWCore/Framework/interface/MakerMacros.h"

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

// To convert detId to subdet/layer number:
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
//#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
//#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
//#include "DataFormats/SiStripDetId/interface/TECDetId.h"
//#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
//#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h" //GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"



#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
//#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"


#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryBuilder.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include <TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h>
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
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

#include "FWCore/Framework/interface/ConsumesCollector.h"


// Flag for new tracking rechis, has to be ON for pre7 and later
//#define NEW_TRACKINGRECHITS  // For V71X_pre7 and later

//
// class declaration:
//
namespace {

    class myCounters{
        public:
            static int neve;
            static unsigned int prevrun;
    };

    int myCounters::neve = 0;
    unsigned int myCounters::prevrun = 0;
    //

} // namespace

class Pixel : public edm::EDAnalyzer{

    public:
        explicit Pixel(const edm::ParameterSet& // ,edm::ConsumesCollector&&
                );
        ~Pixel();

    private:
        virtual void beginJob() ;
        virtual void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
        virtual void analyze(const edm::Event&, const edm::EventSetup&// ,edm::ConsumesCollector&&
                );
        virtual void endJob() ;

        edm::InputTag _triggerSrc;
        std::string _ttrhBuilder;
        HLTConfigProvider HLTConfig;
        bool singleParticleMC;
        int _OC_beginning;
        int _OC_end;
        edm::EDGetTokenT<reco::BeamSpot>  t_offlineBeamSpot_;
        edm::EDGetTokenT<reco::VertexCollection> t_offlinePrimaryVertices_ ;
        edm::EDGetTokenT<edm::TriggerResults> t_triggerSrc_ ;
        edm::EDGetTokenT<reco::TrackCollection>  t_generalTracks_;
        edm::EDGetTokenT<edm::View<reco::PFMET>> t_pfMet_;
        // ----------member data:

        TTree *tree;
        Double_t layer1dx, layer1dz, layer2dx, layer2dz, layer3dx, layer3dz, layer4dx, layer4dz;
        Double_t trkEta, trkPt;
        Int_t layer1EdgeTypeY, layer2EdgeTypeY, layer3EdgeTypeY, layer4EdgeTypeY;
        Int_t layer1EdgeTypeX, layer2EdgeTypeX, layer3EdgeTypeX, layer4EdgeTypeX;
        bool layer1Used2D, layer2Used2D, layer3Used2D, layer4Used2D;
        bool layer1OnEdge, layer2OnEdge, layer3OnEdge, layer4OnEdge;
        Int_t pxn1, pxn2, pxn3, pxn4;
        Int_t layer1xmax, layer1xmin, layer1ymin, layer1ymax, layer1SizeY, layer1SizeX;
        Int_t layer2xmax, layer2xmin, layer2ymin, layer2ymax, layer2SizeY, layer2SizeX;
        Int_t layer3xmax, layer3xmin, layer3ymin, layer3ymax, layer3SizeY, layer3SizeX;
        Int_t layer4xmax, layer4xmin, layer4ymin, layer4ymax, layer4SizeY, layer4SizeX;
        Float_t layer1Charge, layer2Charge, layer3Charge, layer4Charge;

};

//
// constants, enums and typedefs:
//

//
// static data member definitions:
//

//
// constructor:
//
Pixel::Pixel(const edm::ParameterSet& iConfig// , edm::ConsumesCollector && ic
        )
{
    std::cout << "Pixel constructed\n";
    _triggerSrc = iConfig.getParameter<edm::InputTag>("triggerSource");
    _ttrhBuilder = iConfig.getParameter<std::string>("ttrhBuilder");
    singleParticleMC  = iConfig.getUntrackedParameter<bool>("singleParticleMC",false);
    std::cout<<_triggerSrc<<" "<<_triggerSrc.label()<<" "<<_triggerSrc.process()<<" "
        <<_triggerSrc.instance()<<" "<<std::endl;



    //Definition of parameters
    t_triggerSrc_ = consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggerSource"));
    t_offlineBeamSpot_ =    consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
    t_offlinePrimaryVertices_ =   consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
    t_generalTracks_= consumes<reco::TrackCollection> (edm::InputTag("generalTracks"));//"generalTracks"));
    t_pfMet_= consumes< edm::View<reco::PFMET>>(edm::InputTag("pfMet"));
    _OC_beginning=iConfig.getParameter<int>("orbit_beginning");
    _OC_end=iConfig.getParameter<int>("orbit_end");
    // tok_caloHH_ = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits", "HcalHits"));
    edm::Service<TFileService> fsT;
    tree = fsT->make<TTree>("tree", "tree");
    tree->Branch("layer1dx", &layer1dx);
    tree->Branch("layer1dz", &layer1dz);
    tree->Branch("layer2dx", &layer2dx);
    tree->Branch("layer2dz", &layer2dz);
    tree->Branch("layer3dx", &layer3dx);
    tree->Branch("layer3dz", &layer3dz);
    tree->Branch("layer4dx", &layer4dx);
    tree->Branch("layer4dz", &layer4dz);
    tree->Branch("pxn1", &pxn1);
    tree->Branch("pxn2", &pxn2);
    tree->Branch("pxn3", &pxn3);
    tree->Branch("pxn4", &pxn4);
    tree->Branch("trkEta", &trkEta);
    tree->Branch("trkPt", &trkPt);
    tree->Branch("layer1EdgeTypeY", &layer1EdgeTypeY);
    tree->Branch("layer2EdgeTypeY", &layer2EdgeTypeY);
    tree->Branch("layer3EdgeTypeY", &layer3EdgeTypeY);
    tree->Branch("layer4EdgeTypeY", &layer4EdgeTypeY);
    tree->Branch("layer1EdgeTypeX", &layer1EdgeTypeX);
    tree->Branch("layer2EdgeTypeX", &layer2EdgeTypeX);
    tree->Branch("layer3EdgeTypeX", &layer3EdgeTypeX);
    tree->Branch("layer4EdgeTypeX", &layer4EdgeTypeX);
    tree->Branch("layer1OnEdge", &layer1OnEdge);
    tree->Branch("layer2OnEdge", &layer2OnEdge);
    tree->Branch("layer3OnEdge", &layer3OnEdge);
    tree->Branch("layer4OnEdge", &layer4OnEdge);
    tree->Branch("layer1Used2D", &layer1Used2D);
    tree->Branch("layer2Used2D", &layer2Used2D);
    tree->Branch("layer3Used2D", &layer3Used2D);
    tree->Branch("layer4Used2D", &layer4Used2D);

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

    tree->Branch("layer4SizeX", &layer4SizeX);
    tree->Branch("layer4SizeY", &layer4SizeY);
    tree->Branch("layer4xmax", &layer4xmax);
    tree->Branch("layer4xmin", &layer4xmin);
    tree->Branch("layer4ymax", &layer4ymax);
    tree->Branch("layer4ymin", &layer4ymin);

    tree->Branch("layer1Charge", &layer1Charge);
    tree->Branch("layer2Charge", &layer2Charge);
    tree->Branch("layer3Charge", &layer3Charge);
    tree->Branch("layer4Charge", &layer4Charge);
}
//
// destructor:
//
Pixel::~Pixel()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions:
// method called once each job just before starting event loop
//
void Pixel::beginJob()
{
}

//----------------------------------------------------------------------
// method called for each event:

void Pixel::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
    const int run = iRun.run();


}

void Pixel::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup // , edm::ConsumesCollector && ic
        ){
    if((int)iEvent.orbitNumber() >= (int)_OC_beginning && (int)iEvent.orbitNumber() <= (int)_OC_end ){
        using namespace std;
        using namespace edm;
        using namespace reco;
        using namespace math;

        const double pi = 4*atan(1);
        const double wt = 180/pi;
        const double twopi = 2*pi;
        const double pihalf = 2*atan(1);
        //const double sqrtpihalf = sqrt(pihalf);

        myCounters::neve++;
        //FOR NOW only do 1/3 events to get jobs to finish

        /*
        if(myCounters::neve % 3 != 1){
            //printf("Skipping event %i \n", myCounters::neve);
            return;
        }
        */



        if( myCounters::prevrun != iEvent.run() ){

            time_t unixZeit = iEvent.time().unixTime();
            cout << "new run " << iEvent.run();
            cout << ", LumiBlock " << iEvent.luminosityBlock();
            cout << " taken " << ctime(&unixZeit); // ctime has endline

            myCounters::prevrun = iEvent.run();

        }// new run

        int idbg = 0;  // printout for the first few events
        if( myCounters::neve < 1 ) idbg = 1;

        int jdbg = 0; // special printout

        if( idbg ) {
            cout << endl;
            cout << "run " << iEvent.run();
            cout << ", LumiBlock " << iEvent.luminosityBlock();
            cout << ", event " << iEvent.eventAuxiliary().event();
            time_t unixZeit = iEvent.time().unixTime();
            cout << ", taken " << ctime(&unixZeit); // ctime has endline
        }


        //----------------- HLT -------------------------------------------

        // // Trigger information, do only if the container is defined
        if(_triggerSrc.label()!="" && idbg==1 ) {
            edm::Handle<edm::TriggerResults> triggerResults;

            iEvent.getByToken( t_triggerSrc_ , triggerResults);

            //   iEvent.getByLabel(_triggerSrc, triggerResults);
            assert(triggerResults->size() == HLTConfig.size());

            const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
            for(unsigned int i = 0; i < triggerResults->size(); ++i){
                std::string triggerName = triggerNames.triggerName(i);
                // this does not work in 74X
                //std::pair<int, int> prescale = HLTConfig.prescaleValues(iEvent, iSetup, triggerName);
                // this compiles&runs but I am not sure how to use it?
                // std::pair<std::vector<std::pair<std::basic_string<char>, int> >, int>
                // 	prescale = HLTConfig.prescaleValuesInDetail(iEvent, iSetup, triggerName);
                // std::cout << i << ": " << triggerName << ", " << prescale.second << std::endl;
            }
        }

        //--------------------------------------------------------------------
        // beam spot:
        edm::Handle<reco::BeamSpot> rbs;
        iEvent.getByToken( t_offlineBeamSpot_, rbs );

        XYZPoint bsP = XYZPoint(0,0,0);
        int ibs = 0;

        if( rbs.failedToGet() ) return;
        if( ! rbs.isValid() ) return;

        ibs = 1;
        bsP = XYZPoint( rbs->x0(), rbs->y0(), rbs->z0() );

        double bx = rbs->BeamWidthX();
        double by = rbs->BeamWidthY();

        if( idbg ){
            cout << "beam spot x " << rbs->x0();
            cout << ", y " << rbs->y0();
            cout << ", z " << rbs->z0();
            cout << endl;
        }

        //--------------------------------------------------------------------
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

        // need vertex global point for tracks
        // from #include "DataFormats/GeometryVector/interface/GlobalPoint.h"
        // Global points are three-dimensional by default
        // typedef Global3DPoint  GlobalPoint;
        // typedef Point3DBase< float, GlobalTag> Global3DPoint;

        XYZPoint vtxN = XYZPoint(0,0,0);
        XYZPoint vtxP = XYZPoint(0,0,0);

        double bestNdof = 0;
        double maxSumPt = 0;
        Vertex bestPvx;

        for( VertexCollection::const_iterator iVertex = vertices->begin();
                iVertex != vertices->end(); ++iVertex ) {

            if( ! iVertex->isValid() )
                continue;

            else {

                if( iVertex->isFake() )
                    continue;

                else{
                    if( idbg ){
                        cout << "vertex";
                        cout << ": x " << iVertex->x();
                        cout << ", y " << iVertex->y();
                        cout << ", z " << iVertex->z();
                        cout << ", ndof " << iVertex->ndof();
                        cout << ", sumpt " << iVertex->p4().pt();
                        cout << endl;
                    }


                    if( iVertex->ndof() > bestNdof ) {
                        bestNdof = iVertex->ndof();
                        vtxN = XYZPoint( iVertex->x(), iVertex->y(), iVertex->z() );
                    }


                    if( iVertex->p4().pt() > maxSumPt ) {
                        maxSumPt = iVertex->p4().pt();
                        vtxP = XYZPoint( iVertex->x(), iVertex->y(), iVertex->z() );
                        bestPvx = *iVertex;
                    }
                }// non-fake
            }//valid
        } // loop over vertices


        //if( maxSumPt < 1 ) return;
        if(!singleParticleMC && maxSumPt < 1 ) return;

        if( maxSumPt < 1 ) vtxP = vtxN;


        double xBS = 0;
        double yBS = 0;
        if( ibs ) {
            xBS = bsP.x();
            yBS = bsP.y();
        }
        else {
            xBS = vtxP.x();
            yBS = vtxP.y();
        }

        //--------------------------------------------------------------------
        // get a fitter to refit TrackCandidates, the same fitter as used in standard reconstruction:
        // Fitter = cms.string('KFFittingSmootherWithOutliersRejectionAndRK'),
        // KalmanFilter
        // RungeKutta


        //#ifdef NEW_TRACKINGRECHITS

        // Fitter
        edm::ESHandle<TrajectoryFitter> aFitter;
        iSetup.get<TrajectoryFitter::Record>().get("KFFittingSmootherWithOutliersRejectionAndRK",aFitter);
        std::unique_ptr<TrajectoryFitter> theFitter = aFitter->clone();

        //----------------------------------------------------------------------------
        // Transient Rechit Builders
        edm::ESHandle<TransientTrackBuilder> theB;
        iSetup.get<TransientTrackRecord>().get( "TransientTrackBuilder", theB );

        // Transient rec hits:
        ESHandle<TransientTrackingRecHitBuilder> hitBuilder;
        iSetup.get<TransientRecHitRecord>().get( _ttrhBuilder, hitBuilder );

        // Cloner, New from 71Xpre7
        const TkTransientTrackingRecHitBuilder * builder =
            static_cast<TkTransientTrackingRecHitBuilder const *>(hitBuilder.product());
        //dynamic_cast<TkTransientTrackingRecHitBuilder const *>(hitBuilder.product());
        auto hitCloner = builder->cloner();
        //static_cast<TkTransientTrackingRecHitBuilder const *>(hitBuilder.product())->cloner();

        theFitter->setHitCloner(&hitCloner);
        //#endif

        //--------------------------------------------------------------------
        // TrackPropagator:
        edm::ESHandle<Propagator> prop;
        iSetup.get<TrackingComponentsRecord>().get( "PropagatorWithMaterial", prop );
        const Propagator* thePropagator = prop.product();


        //--------------------------------------------------------------------
        // tracks:
        Handle<reco::TrackCollection> tracks;
        iEvent.getByToken( t_generalTracks_, tracks );
        //iEvent.getByLabel( "generalTracks", tracks );

        if( tracks.failedToGet() ) return;
        if( !tracks.isValid() ) return;
        // cout << "  tracks " << tracks->size();
        if( idbg ){
            cout << "  tracks " << tracks->size();
            cout << endl;
        }


        //----------------------------------------------------------------------------
        // get tracker geometry:

        edm::ESHandle<TrackerGeometry> pTG;
        iSetup.get<TrackerDigiGeometryRecord>().get( pTG );

        if( ! pTG.isValid() ) {
            cout << "Unable to find TrackerDigiGeometry. Return\n";
            return;
        }

        //----------------------------------------------------------------------------
        // Tracks:

        double sumpt = 0;
        double sumq = 0;
        Surface::GlobalPoint origin = Surface::GlobalPoint(0,0,0);
        for( TrackCollection::const_iterator iTrack = tracks->begin();
                iTrack != tracks->end(); ++iTrack ) {

            // cpt = cqRB = 0.3*R[m]*B[T] = 1.14*R[m] for B=3.8T
            // D = 2R = 2*pt/1.14
            // calo: D = 1.3 m => pt = 0.74 GeV/c

            double pt = iTrack->pt();
            double pp = iTrack->p();

            trkPt = pt;
            trkEta = iTrack->eta();

            //if( pt < 0.75 ) continue;// curls up
            //if( pt < 1.75 ) continue;// want sharper image

            //float tmp = abs(iTrack->dxy(vtxP))/iTrack->dxyError();
            //cout<<pt<<" "<<abs(iTrack->dxy(vtxP))<<" "<<iTrack->dxyError()<<" "<<tmp<<endl;


            if(!singleParticleMC &&
                    (abs( iTrack->dxy(vtxP) ) > 5*iTrack->dxyError()) ) continue; // not prompt


            sumpt += pt;
            sumq += iTrack->charge();

            //double logpt = log(pt) / log(10);


            const reco::HitPattern& hp = iTrack->hitPattern();

            if( idbg ) {
                cout << endl;
                cout << "Track " << distance( tracks->begin(), iTrack );
                cout << ": pt " << iTrack->pt();
                cout << ", eta " << iTrack->eta();
                cout << ", phi " << iTrack->phi()*wt;
                cout << setprecision(1);
                cout << ", dxyv " << iTrack->dxy(vtxP)*1E4 << " um";
                cout << ", dzv " << iTrack->dz(vtxP)*1E1 << " mm";
                cout << setprecision(4);
                //cout << ", hits " << hp.numberOfHits(HitPattern::TRACK_HITS);/// assuming HitCategory =TRACK_HITS = 0
                cout << ", valid " << hp.numberOfValidTrackerHits();
                cout << endl;
            }

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

            xBS = xBeam;//improve
            yBS = yBeam;//beam tilt taken into account

            double bcap = iTrack->dxy(blP);//impact parameter to beam
            double edca = iTrack->dxyError();
            double ebca = sqrt( edca*edca + bx*by );//round beam
            //double sbca = bcap / ebca;//impact parameter significance

            if( hp.trackerLayersWithMeasurement() < 7 ) continue; // select only tracks which go into the strips



            // transient track:
            TransientTrack tTrack = theB->build(*iTrack);

            TrajectoryStateOnSurface initialTSOS = tTrack.innermostMeasurementState();

            double kap = tTrack.initialFreeState().transverseCurvature();

            TrajectoryStateClosestToPoint tscp = tTrack.trajectoryStateClosestToPoint( origin );

            if( tscp.isValid() ) {

                kap = tscp.perigeeParameters().transverseCurvature();
                phi = tscp.perigeeParameters().phi();
                dca = tscp.perigeeParameters().transverseImpactParameter();
                tet = tscp.perigeeParameters().theta();
                z0  = tscp.perigeeParameters().longitudinalImpactParameter();
                dip = pihalf - tet;


            }//tscp valid

            double cf = cos(phi);
            double sf = sin(phi);
            //double xdca =  dca * sf;
            //double ydca = -dca * cf;

            //double tt = tan(tet);

            double rinv = -kap; // Karimaki
            double rho = 1/kap;
            double erd = 1.0 - kap*dca;
            double drd = dca * ( 0.5*kap*dca - 1.0 ); // 0.5 * kap * dca**2 - dca;
            double hkk = 0.5*kap*kap;

            // track w.r.t. beam (cirmov):

            double dp = -xBS*sf + yBS*cf + dca;
            double dl = -xBS*cf - yBS*sf;
            double sa = 2*dp + rinv * ( dp*dp + dl*dl );
            double dcap = sa / ( 1 + sqrt(1 + rinv*sa) );// distance to beam
            double ud = 1 + rinv*dca;
            double phip = atan2( -rinv*xBS + ud*sf, rinv*yBS + ud*cf );//direction

            // track at R(PXB1), from FUNPHI, FUNLEN:

            double R1 = 3; //34.4; //3.0; // PXB1-4.4 in run 1 -3.0 in phase1

            double s1 = 0;
            double fpos1 = phi - pihalf;

            if( R1 >= abs(dca) ) {

                // sin(delta phi):

                double sindp = ( 0.5*kap * (R1*R1 + dca*dca) - dca ) / (R1*erd);
                fpos1 = phi + asin(sindp); // phi position at R1

                // sin(alpha):

                double sina = R1*kap * sqrt( 1.0 - sindp*sindp );

                // s = alpha / kappa:

                if( sina >= 1.0 )
                    s1 = pi / kap;
                else{
                    if( sina <= -1.0 )
                        s1 = -pi / kap;
                    else
                        s1 = asin(sina) / kap;//always positive
                }

                // Check direction: limit to half-turn

                if( hkk * ( R1*R1 - dca*dca ) > erd ) s1 = pi/abs(kap) - s1; // always positive

            }// R1 > dca

            if( fpos1 > pi ) fpos1 -= twopi;
            else if( fpos1 < -pi ) fpos1 += twopi;

            double zR1 = z0 + s1*tan(dip); // z at R1

            //--------------------------------------------------------------------------
            // loop over tracker detectors:

            double xcrss[99];
            double ycrss[99];
            double zcrss[99];
            int ncrss = 0;

            for( TrackerGeometry::DetContainer::const_iterator idet = pTG->dets().begin();
                    idet != pTG->dets().end(); ++idet ) {

                DetId mydetId = (*idet)->geographicalId();
                uint32_t mysubDet = mydetId.subdetId();

                if( mysubDet != PixelSubdetector::PixelBarrel ) continue;

                /*
                   cout << ": PXB layer " << PXBDetId(mydetId).layer();
                   cout << ", ladder " << PXBDetId(mydetId).ladder();
                   cout << ", module " << PXBDetId(mydetId).module();
                   cout << ", at R1 " << (*idet)->position().perp();
                   cout << ", F " << (*idet)->position().barePhi()*wt;
                   cout << ", z " << (*idet)->position().z();
                   cout << endl;
                   */

                if(// PXBDetId(mydetId).layer()
                        tTopo->pxbLayer(mydetId) == 1 ) {

                    double dz = zR1 - (*idet)->position().z();

                    if( abs(dz) > 4.0 ) continue;

                    double df = fpos1 - (*idet)->position().barePhi();//track phi at R1 vs ladder phi

                    if( df > pi ) df -= twopi;
                    else if( df < -pi ) df += twopi;

                    if( abs(df)*wt > 36.0 ) continue;//coarse matching track to ladder

                    // normal vector: includes alignment (varies from module to module along z on one ladder)
                    // neighbouring ladders alternate with inward/outward orientation

                    /*
                       cout << "normal";
                       cout << ": x " << (*idet)->surface().normalVector().x();
                       cout << ", y " << (*idet)->surface().normalVector().y();
                       cout << ", z " << (*idet)->surface().normalVector().z();
                       cout << ", f " << (*idet)->surface().normalVector().barePhi()*wt;
                       cout << endl;
                       */

                    double phiN = (*idet)->surface().normalVector().barePhi();//normal vector

                    double phidet = phiN - pihalf;// orientation of sensor plane in x-y
                    double ux = cos(phidet);// vector in sensor plane
                    double uy = sin(phidet);
                    double x = (*idet)->position().x();
                    double y = (*idet)->position().y();

                    // intersect helix with line: FUNRXY (in FUNPHI) from V. Blobel
                    // factor f for intersect point (x + f*ux, y + f*uy)

                    double a =                                 kap * ( ux*ux + uy*uy ) * 0.5;
                    double b =       erd * ( ux*sf - uy*cf ) + kap * ( ux*x + uy*y );
                    double c = drd + erd * (  x*sf -  y*cf ) + kap * (  x*x +  y*y ) * 0.5;
                    double dis = b*b - 4.0*a*c;
                    double f = 0;

                    if( dis > 0 ) {

                        dis = sqrt(dis);
                        double f1 = 0;
                        double f2 = 0;

                        if( b < 0 ) {
                            f1 = 0.5 * ( dis - b ) / a;
                            f2 = 2.0 * c / ( dis - b );
                        }
                        else{
                            f1 = -0.5 * ( dis + b ) / a;
                            f2 = -2.0 * c / ( dis + b );
                        }

                        f = f1;
                        if( abs(f2) < abs(f1) ) f = f2;

                    }//dis

                    xcrss[ncrss] = x + f*ux;
                    ycrss[ncrss] = y + f*uy;
                    double r = sqrt( xcrss[ncrss]*xcrss[ncrss] + ycrss[ncrss]*ycrss[ncrss] );

                    double s = 0;

                    if( r >= abs(dca) ) {
                        double sindp = ( 0.5*kap * ( r*r + dca*dca ) - dca ) / (r*erd);
                        double sina = r*kap * sqrt( 1.0 - sindp*sindp );
                        if( sina >= 1.0 )
                            s = pi / kap;
                        else {
                            if( sina <= -1.0 )
                                s = -pi / kap;
                            else
                                s = asin(sina) / kap;
                        }
                        if( hkk * ( r*r - dca*dca ) > erd ) s = pi/abs(kap) - s;
                    }

                    zcrss[ncrss] = z0 + s*tan(dip); // z at r
                    ncrss++;

                }//PXB1

            }//idet

            //--------------------------------------------------------------------------
            // rec hits from track extra:

            if( iTrack->extra().isNull() ) continue;//next track
            if( ! iTrack->extra().isAvailable() ) continue;//next track


            double rmin = 99.9;
            uint32_t innerDetId = 0;

            double xPXB1 = 0;
            double yPXB1 = 0;
            double zPXB1 = 0;
            double uPXB1 = 0;
            //double vPXB1 = 0;
            //double fPXB1 = 0;


            double xPXB2 = 0;
            double yPXB2 = 0;
            double zPXB2 = 0;
            double uPXB2 = 0;
            double vPXB2 = 0;

            double xPXB3 = 0;
            double yPXB3 = 0;
            double zPXB3 = 0;
            double uPXB3 = 0;

            double xPXB4 = 0;
            double yPXB4 = 0;
            double zPXB4 = 0;
            //double uPXB4 = 0;


            double vPXB3 = 0;

            int n1 = 0;
            int n2 = 0;
            int n3 = 0;
            int n4 = 0;
            double phiN1 = 0;
            double phiN2 = 0;
            double phiN3 = 0;
            double phiN4 = 0;
            //double clch1 = 0;
            //double clch3 = 0;
            //int ncol1 = 0;
            //int ncol2 = 0;
            //int ncol3 = 0;
            //double etaX1 = 0;
            //double etaX3 = 0;
            //double cogp1 = 0;
            double cogp2 = 0;
            double cogp3 = 0;
            double xmid2 = 0;
            double ymid2 = 0;
            const GeomDet * det2 = NULL;
            //int ilad2 = 0;
            int zmin2 = 0;
            int zmax2 = 0;

            double xmid3 = 0;
            double ymid3 = 0;
            const GeomDet * det3 = NULL;
            //int ilad3 = 0;
            //int xmin3 = 0;
            //int xmax3 = 0;
            //int zmin3 = 0;
            //int zmax3 = 0;

            //int nTIB1 = 0;

            // std::cout << " imod 0 "<< imod<< std::endl;

            edm::OwnVector<TrackingRecHit> recHitVector; // for seed

            std::vector<TransientTrackingRecHit::RecHitPointer> myTTRHvec;

            Trajectory::RecHitContainer coTTRHvec; // for fit, constant

            // loop over recHits on this track:

            for( trackingRecHit_iterator irecHit = iTrack->recHitsBegin();
                    irecHit != iTrack->recHitsEnd(); ++irecHit ) {

                DetId detId = (*irecHit)->geographicalId();
                uint32_t subDet = detId.subdetId();

                // enum Detector { Tracker=1, Muon=2, Ecal=3, Hcal=4, Calo=5 };
                if( detId.det() != 1 ){
                    cout << "rec hit ID = " << detId.det() << " not in tracker!?!?\n";
                    continue;
                }

                int ilay = tTopo->pxbLayer(detId);//PXBDetId(detId).layer();
                recHitVector.push_back( (*irecHit)->clone() );

                int xmin = 0;
                int xmax = 0;
                int ymin = 0;
                int ymax = 0;

                double cogp = 0;

                //int icol = 0;
                //int irow = 0;
                int ncol = 0;
                int nrow = 0;
                double clch = 0;

                //bool halfmod = 0;

                double Q_f_X = 0.0;//first
                double Q_l_X = 0.0;//last
                double Q_m_X = 0.0;//middle
                double etaX = 0;

                double Q_f_Y = 0.0;//first
                double Q_l_Y = 0.0;//last
                double Q_m_Y = 0.0;//middle


                if( (*irecHit)->isValid() ) {

                    // enum SubDetector{ PixelBarrel=1, PixelEndcap=2 };
                    // enum SubDetector{ TIB=3, TID=4, TOB=5, TEC=6 };

                    //if( idbg ) cout << "  hit in " << subDet << endl;


                    // cast to SiPixelRecHit:
                    // TrackingRecHit -> RecHit2DLocalPos -> BaseSiTrackerRecHit2DLocalPos -> SiPixelRecHit

                    // Examine the detector
                    if( subDet == 1 ){ // PXB

                        int ilay = tTopo->pxbLayer(detId);//PXBDetId(detId).layer();
                        int ilad = tTopo->pxbLadder(detId);//PXBDetId(detId).ladder();

                        if( idbg ) {
                            cout << "  layer  " << tTopo->pxbLayer(detId);// PXBDetId(detId).layer();
                            cout << ", ladder " << tTopo->pxbLadder(detId);// PXBDetId(detId).ladder();
                            cout << ", module " << tTopo->pxbModule(detId);//PXBDetId(detId).module();
                            cout << endl;
                        }

                        const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>( &*(*irecHit) );
                        // cout << transRecHit->localPositionError().xx()<< endl;


                        // pixel cluster:
                        // TrackingRecHit -> RecHit2DLocalPos -> BaseSiTrackerRecHit2DLocalPos -> SiPixelRecHit -> SiPixelCluster
                        edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const & clust = pixhit->cluster();

                        if(ilay == 1){
                            //layer1EdgeTypeY = pixhit->edgeTypeY_;
                            //layer1EdgeTypeX = pixhit->edgeTypeX_;
                            //layer1Used2D = pixhit->filled_from_2d;
                            layer1OnEdge = pixhit->isOnEdge();
                            layer1Charge = clust->charge();
                            layer1SizeX = clust->sizeX();
                            layer1SizeY = clust->sizeY();
                            layer1xmin = clust->minPixelRow();
                            layer1xmax = clust->maxPixelRow();
                            layer1ymin = clust->minPixelCol();
                            layer1ymax = clust->maxPixelCol();
                        }
                        else if(ilay ==2){
                            //layer2EdgeTypeY = pixhit->edgeTypeY_;
                            //layer2EdgeTypeX = pixhit->edgeTypeX_;
                            //layer2Used2D = pixhit->filled_from_2d;
                            layer2OnEdge = pixhit->isOnEdge();
                            layer2Charge = clust->charge();
                            layer2SizeX = clust->sizeX();
                            layer2SizeY = clust->sizeY();
                            layer2xmin = clust->minPixelRow();
                            layer2xmax = clust->maxPixelRow();
                            layer2ymin = clust->minPixelCol();
                            layer2ymax = clust->maxPixelCol();
                        }
                        else if(ilay ==3){
                            //layer3EdgeTypeY = pixhit->edgeTypeY_;
                            //layer3EdgeTypeX = pixhit->edgeTypeX_;
                            //layer3Used2D = pixhit->filled_from_2d;
                            layer3OnEdge = pixhit->isOnEdge();
                            layer3Charge = clust->charge();
                            layer3SizeX = clust->sizeX();
                            layer3SizeY = clust->sizeY();
                            layer3xmin = clust->minPixelRow();
                            layer3xmax = clust->maxPixelRow();
                            layer3ymin = clust->minPixelCol();
                            layer3ymax = clust->maxPixelCol();
                        }
                        else if(ilay ==4){
                            //layer4EdgeTypeY = pixhit->edgeTypeY_;
                            //layer4EdgeTypeX = pixhit->edgeTypeX_;
                            //layer4Used2D = pixhit->filled_from_2d;
                            layer4OnEdge = pixhit->isOnEdge();
                            layer4Charge = clust->charge();
                            layer4SizeX = clust->sizeX();
                            layer4SizeY = clust->sizeY();
                            layer4xmin = clust->minPixelRow();
                            layer4xmax = clust->maxPixelRow();
                            layer4ymin = clust->minPixelCol();
                            layer4ymax = clust->maxPixelCol();
                        }


                        if( clust.isNonnull() ) {
                            /*if( ilay ==4) {
                              PXB4_clusSizeX = clust->sizeX();
                              PXB4_clusSizeY = clust->sizeY();
                              }*/
                            if( idbg ) {
                                cout << setprecision(0);
                                cout << "  charge " << clust->charge();
                                cout << setprecision(4);
                                cout << ", cols " << clust->minPixelCol() << " - " << clust->maxPixelCol(); //0..416 = 8*52
                                cout << " = " << clust->size();
                                cout << ", rows " << clust->minPixelRow() << " - " << clust->maxPixelRow();//0..159 left and right
                                cout << " = " << clust->sizeX();
                                cout << endl;
                            }

                            // Fetch the pixels vector from the cluster:
                            const vector<SiPixelCluster::Pixel> & pixelsVec = clust->pixels();

                            // Obtain boundaries in index units:
                            xmin = clust->minPixelRow();
                            xmax = clust->maxPixelRow();
                            ymin = clust->minPixelCol();
                            ymax = clust->maxPixelCol();

                            // cluster matrix:

                            int QQ[9][99];
                            for( int ix = 0; ix < 9; ++ix ){
                                for( int jz = 0; jz < 99; ++jz ){
                                    QQ[ix][jz] = 0;
                                }
                            }
                            double xsum = 0;
                            double qsum = 0;

                            // loop over the pixels:
                            int isize = pixelsVec.size();
                            for( int i = 0;  i < isize; ++i ) {

                                int ix = pixelsVec[i].x - xmin;
                                int jz = pixelsVec[i].y - ymin;
                                if( ix > 8 ) ix = 8;
                                if( jz > 98 ) jz = 98;
                                QQ[ix][jz] = pixelsVec[i].adc;

                                double pix_adc = pixelsVec[i].adc;
                                qsum += pix_adc;
                                xsum += pix_adc * pixelsVec[i].x;

                                // X projection:
                                if( pixelsVec[i].x == xmin )
                                    Q_f_X += pix_adc;
                                else{
                                    if( pixelsVec[i].x == xmax )
                                        Q_l_X += pix_adc;
                                    else
                                        Q_m_X += pix_adc;
                                }


                                // Y projection:
                                if( pixelsVec[i].y == ymin )
                                    Q_f_Y += pix_adc;
                                else{
                                    if( pixelsVec[i].y == ymax )
                                        Q_l_Y += pix_adc;
                                    else
                                        Q_m_Y += pix_adc;
                                }

                            }//loop over pixels

                            etaX = ( Q_f_X - Q_l_X ) / ( Q_f_X + Q_l_X + Q_m_X );

                            cogp = xsum / qsum;

                            clch = clust->charge();//electrons
                            /*icol = clust->minPixelCol();
                              irow = clust->minPixelRow();*/
                            ncol = clust->sizeY();
                            nrow = clust->sizeX();

                            if( ncol > 5 && idbg ){
                                cout << setprecision(1);
                                cout.setf(ios::showpoint);
                                cout.setf(ios::uppercase);
                                //cout.setf(ios::scientific);
                                cout.setf(ios::fixed);
                                cout << "  dip " << setw(5) << dip*wt;
                                cout << setprecision(4);
                                cout << ", layer  " << tTopo->pxbLayer(detId); //PXBDetId(detId).layer();
                                cout << ", ladder " << tTopo->pxbLadder(detId);//PXBDetId(detId).ladder();
                                cout << ", module " << tTopo->pxbModule(detId);//PXBDetId(detId).module();
                                cout << ", x " << xmin << " - " << xmax;
                                cout << ", z " << ymin << " - " << ymax;
                                cout << endl;
                                for( int ix = 0; ix < min( 9, nrow ); ++ix ){
                                    cout << "    row " << setw(3) << ix + xmin;
                                    for( int jz = 0; jz < min( 99, ncol ); ++jz ){
                                        cout << setw(6) << QQ[ix][jz] / 100;
                                    }
                                    cout << endl;
                                }
                            }//long ncol
                            //
                            //dip -7.E+01, layer  1, ladder 18, module 2, x 148 - 149, z 141 - 148
                            //row 148     0     0     0     0     0     0    75    78
                            //row 149   191   164   111   166    98    96   150     0
                            //
                            //dip 8.E+01, layer  1, ladder 9, module 8, x 21 - 22, z 368 - 375
                            //row  21     0     0     0     0     0     0     0    84
                            //row  22   259   171   113   115   144   106   161     0
                            //
                            //dip -8.E+01, layer  1, ladder 8, module 1, x 8 - 10, z 279 - 290
                            //row   8     0     0     0     0     0     0     0     0     0     0    70    99
                            //row   9     0     0    91   129   131   112    97    97   146   134   107     0
                            //row  10    55   107    56     0     0     0     0     0     0     0     0     0
                            //
                            //dip 8.E+01, layer  1, ladder 2, module 7, x 139 - 140, z 135 - 144
                            //row 139    94   107   173    52     0     0     0     0     0     0
                            //row 140     0     0     0    70   121   128   116   132    99    91
                            //
                            //dip -7.E+01, layer  2, ladder 10, module 2, x 139 - 140, z 408 - 415
                            //row 139    50   397   107   100    71     0     0     0
                            //row 140     0     0     0     0    48   130   246   239


                        }//clust nonNull

                    }//PXB

                }//valid

                // build transient hit:
                //#ifdef NEW_TRACKINGRECHITS
                // for pre7
                auto tmprh =
                    (*irecHit)->cloneForFit(*builder->geometry()->idToDet((**irecHit).geographicalId()));
                auto transRecHit =
                    hitCloner.makeShared(tmprh, initialTSOS);
                //#else

                myTTRHvec.push_back( transRecHit );
                coTTRHvec.push_back( transRecHit );

                if( ! (*irecHit)->isValid() ) continue;

                double xloc = transRecHit->localPosition().x();// 1st meas coord
                double yloc = transRecHit->localPosition().y();// 2nd meas coord or zero
                //double zloc = transRecHit->localPosition().z();// up, always zero, unused

                double vxloc = transRecHit->localPositionError().xx();//covariance
                double vyloc = transRecHit->localPositionError().yy();//covariance

                double gX = transRecHit->globalPosition().x();
                double gY = transRecHit->globalPosition().y();
                double gZ = transRecHit->globalPosition().z();

                if( transRecHit->canImproveWithTrack() ) {//use z from track to apply alignment

                    //if( idbg ) cout << "  try to improve\n";

                    TrajectoryStateOnSurface propTSOS =
                        thePropagator->propagate( initialTSOS, transRecHit->det()->surface() );

                    if( propTSOS.isValid() ){

                        //if( idbg ) cout << "  have propTSOS\n";

                        //#ifdef NEW_TRACKINGRECHITS

                        auto preciseHit = hitCloner.makeShared(tmprh,propTSOS); //pre7

                        //#else
                        xloc = preciseHit->localPosition().x();// 1st meas coord
                        yloc = preciseHit->localPosition().y();// 2nd meas coord or zero
                        // zloc = preciseHit->localPosition().z();// up, always zero

                        vxloc = preciseHit->localPositionError().xx();//covariance
                        //vyloc = preciseHit->localPositionError().yy();//covariance

                        if( idbg && subDet==1 ) {
                            cout << "  improved hit in lay " << tTopo->pxbLayer(detId);//PXBDetId(detId).layer();
                            cout << setprecision(4);
                            cout << ", xloc from " << transRecHit->localPosition().x();
                            cout << " to " << preciseHit->localPosition().x();
                            cout << ", yloc from " << transRecHit->localPosition().y();
                            cout << " to " << preciseHit->localPosition().y();
                            cout << endl;
                            cout<<sqrt(transRecHit->localPositionError().xx())*1E4<<" ";
                            cout<<sqrt(transRecHit->localPositionError().yy())*1E4<<" ";
                            cout<<sqrt(preciseHit->localPositionError().xx())*1E4<<" ";
                            cout<<sqrt(preciseHit->localPositionError().yy())*1E4<<" ";
                            cout << endl;
                        }

                        gX = preciseHit->globalPosition().x();
                        gY = preciseHit->globalPosition().y();
                        gZ = preciseHit->globalPosition().z();

                    }//valid propTSOS
                    else{
                        if( idbg ) cout << "  propTSOS not valid\n";
                    }
                }//canImprove

                double gF = atan2( gY, gX );
                double gR = sqrt( gX*gX + gY*gY );


                if( subDet == PixelSubdetector::PixelBarrel ||
                        subDet == StripSubdetector::TIB ||
                        subDet == StripSubdetector::TOB ) { // barrel

                }

                if( gR < rmin ) {
                    rmin = gR;
                    innerDetId = detId.rawId();
                }

                double phiN = transRecHit->det()->surface().normalVector().barePhi();//normal vector

                double xmid = transRecHit->det()->position().x();
                double ymid = transRecHit->det()->position().y();

                // PXB:

                if( subDet == PixelSubdetector::PixelBarrel ) {

                    double xpix = fmod( xloc+0.82, 0.01 );// xpix = 0..0.01


                    double df = phiN - gF;//normal vector vs position vector: inwards or outwards

                    // take care of cut at +180/-180:

                    if( df > pi ) df -= twopi;
                    else if( df < -pi ) df += twopi;

                    // outward/inward have different Lorentz drift:


                    int ilay = tTopo->pxbLayer(detId);//PXBDetId(detId).layer();
                    int ilad = tTopo->pxbLadder(detId);//PXBDetId(detId).ladder();
                    int imod = tTopo->pxbModule(detId);//PXBDetId(detId).module();

                    if( idbg ) {
                        cout << "  xloc " << xloc;
                        cout << ", cogp " << cogp;
                        double cogx = (cogp + 0.5 - 80) * 0.01 - 0.0054;
                        if( cogp < 79 ) cogx -= 0.01; // big pix
                        if( cogp > 80 ) cogx += 0.01; // big pix
                        cout << ", cogx " << cogx;
                        cout << ", dx = " << cogx - xloc;
                        cout << endl;
                    }


                    if( ilay == 1 ) {
                        n1++;
                        xPXB1 = gX;
                        yPXB1 = gY;
                        zPXB1 = gZ;
                        uPXB1 = xloc;
                        //vPXB1 = yloc;
                        //fPXB1 = sqrt( vyloc );
                        phiN1 = phiN;
                        //clch1 = clch;
                        //ncol1 = ncol;
                        //etaX1 = etaX;
                        //cogp1 = cogp;




                        // if(      ilad ==  5 ) halfmod = 1;
                        // else if( ilad ==  6 ) halfmod = 1;
                        // else if( ilad == 15 ) halfmod = 1;
                        // else if( ilad == 16 ) halfmod = 1;



                        // my crossings:

                        for( int icrss = 0; icrss < ncrss; ++icrss ){

                            double fcrss = atan2( ycrss[icrss], xcrss[icrss] );
                            double df = gF - fcrss;
                            if( df > pi ) df -= twopi;
                            else if( df < -pi ) df += twopi;
                            double du = gR*df;
                            double dz = gZ - zcrss[icrss];


                        }//crss

                    }//PXB1

                    if( ilay == 2 ){

                        n2++;
                        xPXB2 = gX; // precise hit in CMS global coordinates
                        yPXB2 = gY;
                        zPXB2 = gZ;
                        uPXB2 = xloc; // precise hit in local coordinates (w.r.t. sensor center)
                        vPXB2 = yloc;
                        phiN2 = phiN;
                        //ncol2 = ncol;
                        cogp2 = cogp;
                        xmid2 = xmid; // middle of sensor in global CMS coordinates
                        ymid2 = ymid;
                        //ilad2 = ilad;
                        zmin2 = ymin;
                        zmax2 = ymax;

                        det2 = transRecHit->det();


                    }//PXB2

                    if( ilay == 3 ){

                        n3++;
                        xPXB3 = gX;
                        yPXB3 = gY;
                        zPXB3 = gZ;
                        uPXB3 = xloc;
                        vPXB3 = yloc;
                        phiN3 = phiN;
                        //clch3 = clch;
                        //ncol3 = ncol;
                        //etaX3 = etaX;
                        cogp3 = cogp;
                        xmid3 = xmid; // middle of sensor in global CMS coordinates
                        ymid3 = ymid;
                        //ilad3 = ilad;
                        //xmax3 = xmax;
                        //zmin3 = ymin;
                        //zmax3 = ymax;

                        det3 = transRecHit->det();

                    }//PXB3
                    if( ilay == 4 ){

                        n4++;
                        xPXB4 = gX;
                        yPXB4 = gY;
                        zPXB4 = gZ;
                        //uPXB4 = xloc;
                        //vPXB4 = yloc;
                        phiN4 = phiN;
                        //clch4 = clch;
                        //ncol4 = ncol;
                        //etaX4 = etaX;
                        //cogp4 = cogp;


                    }//PXB4

                }//PXB

            }//loop rechits

            //if( pt < 0.75 ) continue;// curls up

            //------------------------------------------------------------------------
            // refit the track:

            //edm::RefToBase<TrajectorySeed> seed = iTrack->seedRef(); // not present in RECO

            //if( idbg ) cout << "  prepare refit\n";

            PTrajectoryStateOnDet PTraj = trajectoryStateTransform::persistentState( initialTSOS, innerDetId );
            const TrajectorySeed seed( PTraj, recHitVector, alongMomentum );

            //if( idbg ) cout << "  have seed\n";

            std::vector<Trajectory> refitTrajectoryCollection = theFitter->fit( seed, coTTRHvec, initialTSOS );

            if( refitTrajectoryCollection.size() > 0 ) { // should be either 0 or 1

                const Trajectory& refitTrajectory = refitTrajectoryCollection.front();

                // Trajectory.measurements:

                Trajectory::DataContainer refitTMs = refitTrajectory.measurements();

                if( idbg ) {
                    cout << "  refitTrajectory has " << refitTMs.size() <<" hits in subdet";
                }

                // hits in subDet:

                if( idbg ) {

                    for( Trajectory::DataContainer::iterator iTM = refitTMs.begin();
                            iTM != refitTMs.end(); iTM++ ) {

                        TransientTrackingRecHit::ConstRecHitPointer iTTRH = iTM->recHit();
                        if( iTTRH->hit()->isValid() ){
                            cout << "  " << iTTRH->geographicalId().subdetId();
                        }
                    }
                    cout << endl;

                    cout << "         pt " << refitTrajectory.geometricalInnermostState().globalMomentum().perp();
                    cout << ", eta " << refitTrajectory.geometricalInnermostState().globalMomentum().eta();
                    cout << ", phi " << refitTrajectory.geometricalInnermostState().globalMomentum().barePhi()*wt;
                    cout << ", at R " << refitTrajectory.geometricalInnermostState().globalPosition().perp();
                    cout << ", z " << refitTrajectory.geometricalInnermostState().globalPosition().z();
                    cout << ", phi " << refitTrajectory.geometricalInnermostState().globalPosition().barePhi()*wt;
                    cout << endl;

                }//dbg

                // trajectory residuals:

                for( Trajectory::DataContainer::iterator iTM = refitTMs.begin();
                        iTM != refitTMs.end(); iTM++ ) {

                    if( ! iTM->recHit()->isValid() ) continue;

                    DetId detId = iTM->recHit()->geographicalId();

                    uint32_t subDet = detId.subdetId();

                    // enum SubDetector{ PixelBarrel=1, PixelEndcap=2 };
                    // enum SubDetector{ TIB=3, TID=4, TOB=5, TEC=6 };

                    double xHit = iTM->recHit()->localPosition().x(); // primary measurement direction
                    double yHit = iTM->recHit()->localPosition().y(); // always 0 in strips
                    /*
                       int ilay = 0;
                       if( detId.subdetId() == 1 ){
                       ilay = PXBDetId( detId ).layer();
                       }

                       if( subDet == 1 && idbg ){//1=PXB
                       cout << "  PXB layer " << ilay << endl;
                       }
                       */

                    double dx = xHit - iTM->predictedState().localPosition().x();
                    double dy = yHit - iTM->predictedState().localPosition().y();
                    double vxh = iTM->recHit()->localPositionError().xx();//covariance
                    double vxt = iTM->predictedState().localError().positionError().xx();//

                    if( subDet == 4 && idbg ){//4=TID
                        cout << "  predictdStateResid = " << dx*1E4 << " um";
                        cout << ", eh = " << sqrt(vxh)*1E4 << " um";
                        cout << ", et = " << sqrt(vxt)*1E4 << " um";
                        cout << endl;
                    }

                    TrajectoryStateOnSurface combinedPredictedState =
                        TrajectoryStateCombiner().combine( iTM->forwardPredictedState(), iTM->backwardPredictedState() );

                    if( ! combinedPredictedState.isValid() ) continue;//skip hit

                    if( jdbg ) cout << "  have combinedPredictedState\n";

                    //double R = combinedPredictedState.globalPosition().perp();
                    double F = combinedPredictedState.globalPosition().barePhi();
                    double Z = combinedPredictedState.globalPosition().z();

                    double xptch = 0;
                    double yptch = 0;

                    if( subDet <  3 ){//1,2=pixel
                        PixelTopology & pixelTopol = (PixelTopology&) iTM->recHit()->detUnit()->topology();
                        xptch = pixelTopol.pitch().first;
                        yptch = pixelTopol.pitch().second;
                    }

                    dx = xHit - combinedPredictedState.localPosition().x(); //x = primary measurement
                    dy = yHit - combinedPredictedState.localPosition().y(); //
                    vxh = iTM->recHit()->localPositionError().xx();//covariance
                    vxt = combinedPredictedState.localError().positionError().xx();//

                    // angles of incidence:
                    // local z = upwards = normal vector
                    // local x = primary measurement direction
                    // local y = secondary measurement direction

                    double alf_inc = atan2( combinedPredictedState.localDirection().x(), combinedPredictedState.localDirection().z() );
                    double bet_inc = atan2( combinedPredictedState.localDirection().y(), combinedPredictedState.localDirection().z() );

                    double phiinc = alf_inc;
                    if( phiinc > pihalf ) phiinc -= pi;
                    else if( phiinc < -pihalf ) phiinc += pi;

                    if( bet_inc > pihalf ) bet_inc -= pi;
                    else if( bet_inc < -pihalf ) bet_inc += pi;

                    if( subDet == 1 && idbg ){//1=PXB
                        cout << setprecision(1);
                        cout << "  combinedStateResid = " << dx*1E4 << " um";
                        cout << ", eh = " << sqrt(vxh)*1E4 << " um";
                        cout << ", et = " << sqrt(vxt)*1E4 << " um";
                        cout << ", dy = " << dy*1E4 << " um";
                        cout << setprecision(4);
                        cout << ", track at x " << combinedPredictedState.localPosition().x();
                        cout << ", y " << combinedPredictedState.localPosition().y();
                        cout << endl;
                    }

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

                    dx = hitMeasurement.x() - combinedPredictedMeasurement.x(); //in units of pitch
                    dy = hitMeasurement.y() - combinedPredictedMeasurement.y(); //in units of pitch
                    dx = dx * xptch;//convert back into [cm] using local pitch
                    dy = dy * yptch;//[cm]

                    if( jdbg ) cout << "  have combinedPredictedMeasurement\n";

                    if( subDet == 1 && idbg ){//1=PXB
                        cout << setprecision(1);
                        cout << "  topologyStateResid = " << dx*1E4 << " um";
                        cout << setprecision(4);
                        cout << ", hit at x " << hitMeasurement.x();
                        cout << ", y " << hitMeasurement.y();
                        cout << ", track at x " << combinedPredictedMeasurement.x();
                        cout << ", y " << combinedPredictedMeasurement.y();
                        cout << endl;
                    }

                    // plot resids:



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
                    //2012 Surface::LocalPoint lp = stripTopol.localPosition( mp ); // makes no difference in TEC

                    //const GeomDet * myGeomDet = pTG->idToDet( detId );
                    const GeomDet * myGeomDet = iTM->recHit()->det(); // makes no difference in TEC
                    Surface::GlobalPoint gp = myGeomDet->toGlobal( lp );

                    double gX = gp.x();
                    double gY = gp.y();
                    double gZ = gp.z();

                    //double phiN = iTM->recHit()->det()->surface().normalVector().barePhi();//normal vector

                    // if( idbg && subDet == StripSubdetector::TEC ) {
                    //   cout << "  TEC hit ";
                    //   cout << ", side " << tTopo->tecSide(detId);
                    //   cout << ", wheel " << tTopo->tecWheel(detId);
                    //   cout << ", ring " << tTopo->tecRing(detId);
                    //   cout << ", petal " << tTopo->tecPetalNumber(detId);
                    //   cout << ", order " << tTopo->tecOrder(detId);//1=back or 2=front
                    //   cout << ", module " << tTopo->tecModule(detId);
                    //   cout << ", stereo " << tTopo->tecIsStereo(detId);
                    //   cout << ", comp " << myGeomDet->components().size(); //always 0
                    //   cout << endl;
                    // }

                    //2012: overwrite PXB global coordinates once more, using topology:

                    if( subDet == PixelSubdetector::PixelBarrel ) {

                        int ilay = tTopo->pxbLayer(detId); //PXBDetId(detId).layer();

                        if( ilay == 1 ) {
                            xPXB1 = gX;
                            yPXB1 = gY;
                            zPXB1 = gZ;
                        }

                        else if( ilay == 2 ) {
                            xPXB2 = gX;
                            yPXB2 = gY;
                            zPXB2 = gZ;
                        }

                        else if( ilay == 3 ) {
                            xPXB3 = gX;
                            yPXB3 = gY;
                            zPXB3 = gZ;
                        }
                        // no layer 4!
                    }//PXB


                }//loop iTM

            }//refitted trajectory

            //if( pt < 0.75 ) continue;// curls up

            //------------------------------------------------------------------------
            // refit once more, leaving out hit in 2nd PXB:

            // double  clusProb_pxl = -99;
            pxn1 = n1;
            pxn2 = n2;
            pxn3 = n3;
            pxn4 = n4;
            if( n2 > 0 ){

                Trajectory::RecHitContainer nyTTRHvec; // for fit

                for( vector<TransientTrackingRecHit::RecHitPointer>::iterator iTTRH = myTTRHvec.begin();
                        iTTRH != myTTRHvec.end(); ++iTTRH ) {

                    if( idbg == 9 ) {
                        cout << "  hit " << distance( (vector<TransientTrackingRecHit::RecHitPointer>::iterator)myTTRHvec.begin(), iTTRH );
                        if( (*iTTRH)->hit()->isValid() ){
                            cout << ": subdet " << (*iTTRH)->hit()->geographicalId().subdetId();
                            cout << ", weight " << (*iTTRH)->weight(); // default weight is 1

                            cout << endl;
                        }
                        else cout << " not valid\n";
                    }

                    int iuse = 1;
                    if( (*iTTRH)->hit()->isValid() ){
                        if( (*iTTRH)->hit()->geographicalId().subdetId() == 1 ){ //PXB
                            if( PXBDetId( (*iTTRH)->hit()->geographicalId() ).layer() == 2 ) iuse = 0; // skip PXB2: unbiased track
                        }
                    }
                    if( iuse ) nyTTRHvec.push_back( *iTTRH );
                }//copy

                if( idbg ) {
                    cout << "  all hits " << myTTRHvec.size();
                    cout << ", without PXB2 " << nyTTRHvec.size();
                    cout << endl;
                }


                std::vector<Trajectory> refitTrajectoryVec2 = theFitter->fit( seed, nyTTRHvec, initialTSOS );

                if( refitTrajectoryVec2.size() > 0 ) { // should be either 0 or 1

                    const Trajectory& refitTrajectory2 = refitTrajectoryVec2.front();

                    // Trajectory.measurements:

                    Trajectory::DataContainer refitTMvec2 = refitTrajectory2.measurements();

                    if( idbg ) {
                        cout << "  refitTrajectory2 has " << refitTMvec2.size() <<" hits in subdet";
                    }

                    for( Trajectory::DataContainer::iterator iTM = refitTMvec2.begin();
                            iTM != refitTMvec2.end(); iTM++ ) {

                        TransientTrackingRecHit::ConstRecHitPointer iTTRH = iTM->recHit();

                        if( iTTRH->hit()->isValid() ){
                            // const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>( iTTRH->hit() );
                            // if (pixhit->hasFilledProb()){
                            //   cout << "PXL layer 2:  clust prob "<<pixhit->clusterProbability(0)<<endl;
                            //   clusProb_pxl = pixhit->clusterProbability(0);
                            // PXB2_clusterProb
                            // }else     cout << "PXL layer 2: NO  clust prob "<< endl;

                            if( idbg ) cout << "  " << iTTRH->geographicalId().subdetId();

                        }
                    }
                    if( idbg ) cout << endl;

                    if( idbg ) {
                        cout << "  ndof " << refitTrajectory2.ndof();
                        cout << ", found " << refitTrajectory2.foundHits();
                        cout << ", missing " << refitTrajectory2.lostHits();
                        cout << ", chi2 " << refitTrajectory2.chiSquared();
                        cout << ", /ndof " << refitTrajectory2.chiSquared() / refitTrajectory2.ndof();
                        cout << endl;
                        //}

                        //if( idbg ) {
                        cout << "         pt " << refitTrajectory2.geometricalInnermostState().globalMomentum().perp();
                        cout << ", eta " << refitTrajectory2.geometricalInnermostState().globalMomentum().eta();
                        cout << ", phi " << refitTrajectory2.geometricalInnermostState().globalMomentum().barePhi()*wt;
                        cout << ", at R " << refitTrajectory2.geometricalInnermostState().globalPosition().perp();
                        cout << ", z " << refitTrajectory2.geometricalInnermostState().globalPosition().z();
                        cout << ", phi " << refitTrajectory2.geometricalInnermostState().globalPosition().barePhi()*wt;
                        cout << endl;
                }

                // now use it!

                }//refit2
                else
                    if( idbg ) cout << "  no refit\n";

            }// n2 > 0
            //------------------------------------------------------------------------
            //------------------------------------------------------------------------
            // 1-2-4 pixel triplet:

            if( n1*n2*n4 > 0 ) {

                { // let's open a scope, so we can redefine the variables further down

                    if( jdbg ) cout << "  triplet 1+4 -> 2\n";

                    double f2 = atan2( yPXB2, xPXB2 );//position angle in layer 2

                    double ax = xPXB4 - xPXB1;
                    double ay = yPXB4 - yPXB1;
                    double aa = sqrt( ax*ax + ay*ay ); // from 1 to 4

                    double xmid = 0.5 * ( xPXB1 + xPXB4 );
                    double ymid = 0.5 * ( yPXB1 + yPXB4 );
                    double bx = xPXB2 - xmid;
                    double by = yPXB2 - ymid;
                    double bb = sqrt( bx*bx + by*by ); // from mid point to point 2


                    // Author: Johannes Gassner (15.11.1996)
                    // Make track from 2 space points and kappa (cmz98/ftn/csmktr.f)
                    // Definition of the Helix :

                    // x( t ) = X0 + KAPPA^-1 * SIN( PHI0 + t )
                    // y( t ) = Y0 - KAPPA^-1 * COS( PHI0 + t )          t > 0
                    // z( t ) = Z0 + KAPPA^-1 * TAN( DIP ) * t

                    // Center of the helix in the xy-projection:

                    // X0 = + ( DCA - KAPPA^-1 ) * SIN( PHI0 )
                    // Y0 = - ( DCA - KAPPA^-1 ) * COS( PHI0 )

                    // Point 1 must be in the inner layer, 4 in the outer:

                    double r1 = sqrt( xPXB1*xPXB1 + yPXB1*yPXB1 );
                    double r4 = sqrt( xPXB4*xPXB4 + yPXB4*yPXB4 );

                    //	cout << "!!!warn r1 = " << r1 << ", r4 = " << r4 << endl;

                    if( r4-r1 < 2.0 ) cout << "warn r1 = " << r1 << ", r4 = " << r4 << endl;

                    // Calculate the centre of the helix in xy-projection that
                    // transverses the two spacepoints. The points with the same
                    // distance from the two points are lying on a line.
                    // LAMBDA is the distance between the point in the middle of
                    // the two spacepoints and the centre of the helix.

                    // we already have kap and rho = 1/kap

                    double lam = sqrt( -0.25 +
                            rho*rho / ( ( xPXB1 - xPXB4 )*( xPXB1 - xPXB4 ) + ( yPXB1 - yPXB4 )*( yPXB1 - yPXB4 ) ) );

                    // There are two solutions, the sign of kap gives the information
                    // which of them is correct:

                    if( kap > 0 ) lam = -lam;

                    // ( X0, Y0 ) is the centre of the circle
                    // that describes the helix in xy-projection:

                    double x0 =  0.5*( xPXB1 + xPXB4 ) + lam * ( -yPXB1 + yPXB4 );
                    double y0 =  0.5*( yPXB1 + yPXB4 ) + lam * (  xPXB1 - xPXB4 );

                    // Calculate theta:

                    double num = ( yPXB4 - y0 ) * ( xPXB1 - x0 ) - ( xPXB4 - x0 ) * ( yPXB1 - y0 );
                    double den = ( xPXB1 - x0 ) * ( xPXB4 - x0 ) + ( yPXB1 - y0 ) * ( yPXB4 - y0 );
                    double tandip = kap * ( zPXB4 - zPXB1 ) / atan( num / den );
                    double udip = atan(tandip);
                    //double utet = pihalf - udip;

                    // To get phi0 in the right interval one must distinguish
                    // two cases with positve and negative kap:

                    double uphi;
                    if( kap > 0 ) uphi = atan2( -x0,  y0 );
                    else          uphi = atan2(  x0, -y0 );

                    // The distance of the closest approach DCA depends on the sign
                    // of kappa:

                    double udca;
                    if( kap > 0 ) udca = rho - sqrt( x0*x0 + y0*y0 );
                    else          udca = rho + sqrt( x0*x0 + y0*y0 );

                    // angle from first hit to dca point:

                    double dphi = atan( ( ( xPXB1 - x0 ) * y0 - ( yPXB1 - y0 ) * x0 )
                            / ( ( xPXB1 - x0 ) * x0 + ( yPXB1 - y0 ) * y0 ) );

                    double uz0 = zPXB1 + tandip * dphi * rho;


                    // interpolate to middle hit:
                    // cirmov
                    // we already have rinv = -kap

                    double cosphi = cos(uphi);
                    double sinphi = sin(uphi);
                    double dp = -xPXB2*sinphi + yPXB2*cosphi + udca;
                    double dl = -xPXB2*cosphi - yPXB2*sinphi;
                    double sa = 2*dp + rinv * ( dp*dp + dl*dl );
                    double dca2 = sa / ( 1 + sqrt(1 + rinv*sa) );// distance to hit 2
                    double ud = 1 + rinv*udca;
                    double phi2 = atan2( -rinv*xPXB2 + ud*sinphi, rinv*yPXB2 + ud*cosphi );//direction

                    double phiinc = phi2 - phiN2;//angle of incidence in rphi w.r.t. normal vector

                    // phiN alternates inward/outward
                    // reduce phiinc:

                    if( phiinc > pihalf ) phiinc -= pi;
                    else if( phiinc < -pihalf ) phiinc += pi;

                    // arc length:

                    double xx = xPXB2 + dca2 * sin(phi2); // point on track
                    double yy = yPXB2 - dca2 * cos(phi2);

                    double vx = xx - xmid2;//from module center
                    double vy = yy - ymid2;
                    double vv = sqrt( vx*vx + vy*vy );

                    double f0 = uphi;//
                    double kx = kap*xx;
                    double ky = kap*yy;
                    double kd = kap*udca;

                    // Solve track equation for s:

                    double dx = kx - (kd-1)*sin(f0);
                    double dy = ky + (kd-1)*cos(f0);
                    double ks = atan2( dx, -dy ) - f0;// turn angle

                    // Limit to half-turn:

                    if(      ks >  pi ) ks = ks - twopi;
                    else if( ks < -pi ) ks = ks + twopi;

                    double s = ks*rho; // signed
                    double uz2 = uz0 + s*tandip; // track z at R2
                    double dz2 = zPXB2 - uz2;

                    Surface::GlobalPoint gp2( xx, yy, uz2 );
                    Surface::LocalPoint lp2 = det2->toLocal( gp2 );

                    if( idbg ) {
                        std::cout <<"**** local  Point coord ****" <<std::endl;
                        std::cout <<"gp2.x() "<< gp2.x() <<std::endl;
                        std::cout <<"gp2.y() "<< gp2.y() <<std::endl;
                        std::cout <<"gp2.z() "<< gp2.z() <<std::endl;

                        std::cout <<"lp2.x() "<< lp2.x() <<std::endl;
                        std::cout <<"lp2.y() "<< lp2.y() <<std::endl;
                        std::cout <<"lp2.z() "<< lp2.z() <<std::endl;

                        std::cout <<"  uPXB2 = xloc;  precise hit in local coordinates (w.r.t. sensor center)"<< uPXB2 <<std::endl;
                        std::cout <<"  vPXB2 = yloc;precise hit in local coordinates (w.r.t. sensor center)" << vPXB2 <<std::endl;

                    }
                    // local x = phi
                    // local y = z in barrel
                    // local z = radial in barrel (thickness)

                    //double xpix = fmod( uPXB2 + 0.82, 0.01 ); // xpix = 0..0.01 reconstructed
                    //double xpx2 = fmod( uPXB2 + 0.82, 0.02 ); // xpix = 0..0.02 reconstructed
                    //double xpx1 = fmod( uPXB1 + 0.82, 0.01 ); // xpix = 0..0.01 reconstructed
                    //double xpx4 = fmod( uPXB4 + 0.82, 0.01 ); // xpix = 0..0.01 reconstructed

                    //double dpix = fmod( uPXB2 + dca2 + 0.82, 0.01 ); // dpix = 0..0.01 predicted

                    double vpix = fmod( vv, 0.01 ); // vpix = 0..0.01 predicted
                    if( uPXB2 < 0 ) vpix = -vpix; // vv is unsigned distance from module center

                    //double lpix = fmod( lp2.x() + 0.82, 0.01 ); // lpix = 0..0.01 predicted
                    //double tpix = fmod( lp2.x() + 0.82, 0.02 ); // tpix = 0..0.02 predicted

                    //double zpix = fmod( lp2.y() + 3.24, 0.015 ); // zpix = 0..0.015 predicted
                    //double spix = fmod( lp2.y() + 3.24, 0.03  ); // spix = 0..0.03  predicted

                    //int smin = zmin2%52; // 0..51 column along z
                    //int smax = zmax2%52; // 0..51 column along z

                    double cogx = (cogp2 + 0.5 - 80) * 0.01 - 0.0054; // Lorentz shift
                    if( cogp2 < 79 ) cogx -= 0.01; // big pix
                    if( cogp2 > 80 ) cogx += 0.01; // big pix

                    double mpix = fmod( cogx + 0.82, 0.01 ); // mpix = 0..0.01 from cluster COG
                    //double cogdx = cogx - lp2.x(); // residual

                    // hybrid method:

                    //double hybx = uPXB2; // template
                    //if( mpix*1E4 < 20 ) hybx = cogx; // COG
                    //if( mpix*1E4 > 75 ) hybx = cogx;
                    //double hpix = fmod( hybx + 0.82, 0.01 ); // hpix = 0..0.01 from cluster hybrid method
                    //double hybdx = hybx - lp2.x(); // residual

                    /*bool halfmod = 0;
                      if(      ilad2 ==  8 ) halfmod = 1;
                      else if( ilad2 ==  9 ) halfmod = 1;
                      else if( ilad2 == 24 ) halfmod = 1;
                      else if( ilad2 == 25 ) halfmod = 1;
                      */


                    ////check for bias dot flip, unflipped, zplus,zminus,xplu,xminus combination





                    // //////// Modules:1-8


                    // profile of abs(dca) gives mean abs(dca):
                    // mean of abs(Gauss) = 0.7979 * RMS = 1/sqrt(pi/2)
                    // => rms = sqrt(pi/2) * mean of abs (sqrt(pi/2) = 1.2533)
                    // point resolution = 1/sqrt(3/2) * triplet middle residual width
                    // => sqrt(pi/2)*sqrt(2/3) = sqrt(pi/3) = 1.0233, almost one


                    // pt bins:
                    // uPXB2 = local x
                    // vPXB2 = local z
                    // low pt: material
                }

            }


            if (n3*n2*n4>0){

                //------------------------------------------------------------------------
                // triplet 3+4 -> 2:

                {
                    if( jdbg ) cout << "  triplet 3+4 -> 2\n";

                    double f2 = atan2( yPXB2, xPXB2 );//position angle in layer 2

                    double ax = xPXB4 - xPXB3;
                    double ay = yPXB4 - yPXB3;
                    double aa = sqrt( ax*ax + ay*ay ); // from 2 to 3

                    double xmid = 0.5 * ( xPXB3 + xPXB4 );
                    double ymid = 0.5 * ( yPXB3 + yPXB4 );
                    double bx = xPXB2 - xmid;
                    double by = yPXB2 - ymid;
                    double bb = sqrt( bx*bx + by*by ); // from mid point to point 2

                    // Calculate the centre of the helix in xy-projection that
                    // transverses the two spacepoints. The points with the same
                    // distance from the two points are lying on a line.
                    // LAMBDA is the distance between the point in the middle of
                    // the two spacepoints and the centre of the helix.

                    // we already have kap and rho = 1/kap

                    double lam = sqrt( -0.25 +
                            rho*rho / ( ( xPXB3 - xPXB4 )*( xPXB3 - xPXB4 ) + ( yPXB3 - yPXB4 )*( yPXB3 - yPXB4 ) ) );

                    // There are two solutions, the sign of kap gives the information
                    // which of them is correct:

                    if( kap > 0 ) lam = -lam;

                    // ( X0, Y0 ) is the centre of the circle
                    // that describes the helix in xy-projection:

                    double x0 =  0.5*( xPXB3 + xPXB4 ) + lam * ( -yPXB3 + yPXB4 );
                    double y0 =  0.5*( yPXB3 + yPXB4 ) + lam * (  xPXB3 - xPXB4 );

                    // Calculate theta:

                    double num = ( yPXB4 - y0 ) * ( xPXB3 - x0 ) - ( xPXB4 - x0 ) * ( yPXB3 - y0 );
                    double den = ( xPXB3 - x0 ) * ( xPXB4 - x0 ) + ( yPXB3 - y0 ) * ( yPXB4 - y0 );
                    double tandip = kap * ( zPXB4 - zPXB3 ) / atan( num / den );
                    double udip = atan(tandip);
                    //double utet = pihalf - udip;

                    // To get phi0 in the right interval one must distinguish
                    // two cases with positve and negative kap:

                    double uphi;
                    if( kap > 0 ) uphi = atan2( -x0,  y0 );
                    else          uphi = atan2(  x0, -y0 );

                    // The distance of the closest approach DCA depends on the sign
                    // of kappa:

                    double udca;
                    if( kap > 0 ) udca = rho - sqrt( x0*x0 + y0*y0 );
                    else          udca = rho + sqrt( x0*x0 + y0*y0 );

                    // angle from first hit to dca point:

                    double dphi = atan( ( ( xPXB3 - x0 ) * y0 - ( yPXB3 - y0 ) * x0 )
                            / ( ( xPXB3 - x0 ) * x0 + ( yPXB3 - y0 ) * y0 ) );

                    double uz0 = zPXB3 + tandip * dphi * rho;


                    // extrapolate to inner hit:
                    // cirmov
                    // we already have rinv = -kap

                    double cosphi = cos(uphi);
                    double sinphi = sin(uphi);
                    double dp = -xPXB2*sinphi + yPXB2*cosphi + udca;
                    double dl = -xPXB2*cosphi - yPXB2*sinphi;
                    double sa = 2*dp + rinv * ( dp*dp + dl*dl );
                    double dca2 = sa / ( 1 + sqrt(1 + rinv*sa) );// distance to hit 1
                    double ud = 1 + rinv*udca;
                    double phi2 = atan2( -rinv*xPXB2 + ud*sinphi, rinv*yPXB2 + ud*cosphi );//direction

                    double phiinc = phi2 - phiN2;//angle of incidence in rphi w.r.t. normal vector

                    // phiN alternates inward/outward
                    // reduce phiinc:

                    if( phiinc > pihalf ) phiinc -= pi;
                    else if( phiinc < -pihalf ) phiinc += pi;

                    // arc length:

                    double xx = xPXB2 + dca2 * sin(phi2); // point on track
                    double yy = yPXB2 - dca2 * cos(phi2);

                    double f0 = uphi;//
                    double kx = kap*xx;
                    double ky = kap*yy;
                    double kd = kap*udca;

                    // Solve track equation for s:

                    double dx = kx - (kd-1)*sin(f0);
                    double dy = ky + (kd-1)*cos(f0);
                    double ks = atan2( dx, -dy ) - f0;// turn angle

                    //---  Limit to half-turn:

                    if(      ks >  pi ) ks = ks - twopi;
                    else if( ks < -pi ) ks = ks + twopi;

                    double s = ks*rho;// signed
                    double uz2 = uz0 + s*tandip; //track z at R2
                    double dz2 = zPXB2 - uz2;



                    // residual profiles: alignment check


                    // profile of abs(dca) gives mean abs(dca):
                    // mean of abs(Gauss) = 0.7979 * RMS = 1/sqrt(pi/2)
                    // => rms = sqrt(pi/2) * mean of abs (sqrt(pi/2) = 1.2533)
                    // point resolution = 1/sqrt(3/2) * triplet middle residual width
                    // => sqrt(pi/2)*sqrt(2/3) = sqrt(pi/3) = 1.0233, almost one

                }

                // 2 + 3 => 4
                {
                    if( jdbg ) cout << "  triplet 2+3 -> 4\n";

                    double f4 = atan2( yPXB4, xPXB4 );//position angle in layer 4

                    double ax = xPXB3 - xPXB2;
                    double ay = yPXB3 - yPXB2;
                    double aa = sqrt( ax*ax + ay*ay ); // from 2 to 3

                    double xmid = 0.5 * ( xPXB2 + xPXB3 );
                    double ymid = 0.5 * ( yPXB2 + yPXB3 );
                    double bx = xPXB4 - xmid;
                    double by = yPXB4 - ymid;
                    double bb = sqrt( bx*bx + by*by ); // from mid point to point 4

                    // Calculate the centre of the helix in xy-projection that
                    // transverses the two spacepoints. The points with the same
                    // distance from the two points are lying on a line.
                    // LAMBDA is the distance between the point in the middle of
                    // the two spacepoints and the centre of the helix.

                    // we already have kap and rho = 1/kap

                    double lam = sqrt( -0.25 +
                            rho*rho / ( ( xPXB2 - xPXB3 )*( xPXB2 - xPXB3 ) + ( yPXB2 - yPXB3 )*( yPXB2 - yPXB3 ) ) );

                    // There are two solutions, the sign of kap gives the information
                    // which of them is correct:

                    if( kap > 0 ) lam = -lam;

                    // ( X0, Y0 ) is the centre of the circle
                    // that describes the helix in xy-projection:

                    double x0 =  0.5*( xPXB2 + xPXB3 ) + lam * ( -yPXB2 + yPXB3 );
                    double y0 =  0.5*( yPXB2 + yPXB3 ) + lam * (  xPXB2 - xPXB3 );

                    // Calculate theta:

                    double num = ( yPXB3 - y0 ) * ( xPXB2 - x0 ) - ( xPXB3 - x0 ) * ( yPXB2 - y0 );
                    double den = ( xPXB2 - x0 ) * ( xPXB3 - x0 ) + ( yPXB2 - y0 ) * ( yPXB3 - y0 );
                    double tandip = kap * ( zPXB3 - zPXB2 ) / atan( num / den );
                    double udip = atan(tandip);
                    //double utet = pihalf - udip;

                    // To get phi0 in the right interval one must distinguish
                    // two cases with positve and negative kap:

                    double uphi;
                    if( kap > 0 ) uphi = atan2( -x0,  y0 );
                    else          uphi = atan2(  x0, -y0 );

                    // The distance of the closest approach DCA depends on the sign
                    // of kappa:

                    double udca;
                    if( kap > 0 ) udca = rho - sqrt( x0*x0 + y0*y0 );
                    else          udca = rho + sqrt( x0*x0 + y0*y0 );

                    // angle from first hit to dca point:

                    double dphi = atan( ( ( xPXB2 - x0 ) * y0 - ( yPXB2 - y0 ) * x0 )
                            / ( ( xPXB2 - x0 ) * x0 + ( yPXB2 - y0 ) * y0 ) );

                    double uz0 = zPXB2 + tandip * dphi * rho;


                    // extrapolate to outer hit:
                    // cirmov
                    // we already have rinv = -kap

                    double cosphi = cos(uphi);
                    double sinphi = sin(uphi);
                    double dp = -xPXB4*sinphi + yPXB4*cosphi + udca;
                    double dl = -xPXB4*cosphi - yPXB4*sinphi;
                    double sa = 2*dp + rinv * ( dp*dp + dl*dl );
                    double dca4 = sa / ( 1 + sqrt(1 + rinv*sa) );// distance to hit 4
                    double ud = 1 + rinv*udca;
                    double phi4 = atan2( -rinv*xPXB4 + ud*sinphi, rinv*yPXB4 + ud*cosphi );//track direction

                    double phiinc = phi4 - phiN4;//angle of incidence in rphi w.r.t. normal vector

                    // phiN alternates inward/outward
                    // reduce phiinc:

                    if( phiinc > pihalf ) phiinc -= pi;
                    else if( phiinc < -pihalf ) phiinc += pi;

                    // arc length:

                    double xx = xPXB4 + dca4 * sin(phi4); // point on track
                    double yy = yPXB4 - dca4 * cos(phi4);

                    double f0 = uphi;//
                    double kx = kap*xx;
                    double ky = kap*yy;
                    double kd = kap*udca;

                    // Solve track equation for s:

                    double dx = kx - (kd-1)*sin(f0);
                    double dy = ky + (kd-1)*cos(f0);
                    double ks = atan2( dx, -dy ) - f0;// turn angle

                    //---  Limit to half-turn:

                    if(      ks >  pi ) ks = ks - twopi;
                    else if( ks < -pi ) ks = ks + twopi;

                    double s = ks*rho;// signed
                    double uz4 = uz0 + s*tandip; //track z at R4
                    double dz4 = zPXB4 - uz4;

                    layer4dx = dca4 * 1E4;
                    layer4dz = dz4 * 1E4;


                }


                if( n3*n2*n4 > 0 ) {




                    { // let's open a scope, so we can redefine the variables further down

                        if( jdbg ) cout << "  triplet 2+4 -> 3\n";

                        double f3 = atan2( yPXB3, xPXB3 );//position angle in layer 3

                        double ax = xPXB4 - xPXB2;
                        double ay = yPXB4 - yPXB2;
                        double aa = sqrt( ax*ax + ay*ay ); // from 2 to 4

                        double xmid = 0.5 * ( xPXB2 + xPXB4 );
                        double ymid = 0.5 * ( yPXB2 + yPXB4 );
                        double bx = xPXB3 - xmid;
                        double by = yPXB3 - ymid;
                        double bb = sqrt( bx*bx + by*by ); // from mid point to point 3


                        // Author: Johannes Gassner (15.11.1996)
                        // Make track from 3 space points and kappa (cmz98/ftn/csmktr.f)
                        // Definition of the Helix :

                        // x( t ) = X0 + KAPPA^-1 * SIN( PHI0 + t )
                        // y( t ) = Y0 - KAPPA^-1 * COS( PHI0 + t )          t > 0
                        // z( t ) = Z0 + KAPPA^-1 * TAN( DIP ) * t

                        // Center of the helix in the xy-projection:

                        // X0 = + ( DCA - KAPPA^-1 ) * SIN( PHI0 )
                        // Y0 = - ( DCA - KAPPA^-1 ) * COS( PHI0 )

                        // Point 1 must be in the inner layer, 4 in the outer:

                        double r2 = sqrt( xPXB2*xPXB2 + yPXB2*yPXB2 );
                        double r4 = sqrt( xPXB4*xPXB4 + yPXB4*yPXB4 );

                        //	cout << "!!!warn r2 = " << r2 << ", r4 = " << r4 << endl;

                        if( r4-r2 < 2.0 ) cout << "warn r2 = " << r2 << ", r4 = " << r4 << endl;

                        // Calculate the centre of the helix in xy-projection that
                        // transverses the two spacepoints. The points with the same
                        // distance from the two points are lying on a line.
                        // LAMBDA is the distance between the point in the middle of
                        // the two spacepoints and the centre of the helix.

                        // we already have kap and rho = 1/kap

                        double lam = sqrt( -0.25 +
                                rho*rho / ( ( xPXB2 - xPXB4 )*( xPXB2 - xPXB4 ) + ( yPXB2 - yPXB4 )*( yPXB2 - yPXB4 ) ) );

                        // There are two solutions, the sign of kap gives the information
                        // which of them is correct:

                        if( kap > 0 ) lam = -lam;

                        // ( X0, Y0 ) is the centre of the circle
                        // that describes the helix in xy-projection:

                        double x0 =  0.5*( xPXB2 + xPXB4 ) + lam * ( -yPXB2 + yPXB4 );
                        double y0 =  0.5*( yPXB2 + yPXB4 ) + lam * (  xPXB2 - xPXB4 );

                        // Calculate theta:

                        double num = ( yPXB4 - y0 ) * ( xPXB2 - x0 ) - ( xPXB4 - x0 ) * ( yPXB2 - y0 );
                        double den = ( xPXB2 - x0 ) * ( xPXB4 - x0 ) + ( yPXB2 - y0 ) * ( yPXB4 - y0 );
                        double tandip = kap * ( zPXB4 - zPXB2 ) / atan( num / den );
                        double udip = atan(tandip);
                        //double utet = pihalf - udip;

                        // To get phi0 in the right interval one must distinguish
                        // two cases with positve and negative kap:

                        double uphi;
                        if( kap > 0 ) uphi = atan2( -x0,  y0 );
                        else          uphi = atan2(  x0, -y0 );

                        // The distance of the closest approach DCA depends on the sign
                        // of kappa:

                        double udca;
                        if( kap > 0 ) udca = rho - sqrt( x0*x0 + y0*y0 );
                        else          udca = rho + sqrt( x0*x0 + y0*y0 );

                        // angle from first hit to dca point:

                        double dphi = atan( ( ( xPXB2 - x0 ) * y0 - ( yPXB2 - y0 ) * x0 )
                                / ( ( xPXB2 - x0 ) * x0 + ( yPXB2 - y0 ) * y0 ) );

                        double uz0 = zPXB2 + tandip * dphi * rho;


                        // interpolate to middle hit:
                        // cirmov
                        // we already have rinv = -kap

                        double cosphi = cos(uphi);
                        double sinphi = sin(uphi);
                        double dp = -xPXB3*sinphi + yPXB3*cosphi + udca;
                        double dl = -xPXB3*cosphi - yPXB3*sinphi;
                        double sa = 2*dp + rinv * ( dp*dp + dl*dl );
                        double dca3 = sa / ( 1 + sqrt(1 + rinv*sa) );// distance to hit 3
                        double ud = 1 + rinv*udca;
                        double phi3 = atan2( -rinv*xPXB3 + ud*sinphi, rinv*yPXB3 + ud*cosphi );//direction

                        double phiinc = phi3 - phiN3;//angle of incidence in rphi w.r.t. normal vector

                        // phiN alternates inward/outward
                        // reduce phiinc:

                        if( phiinc > pihalf ) phiinc -= pi;
                        else if( phiinc < -pihalf ) phiinc += pi;

                        // arc length:

                        double xx = xPXB3 + dca3 * sin(phi3); // point on track
                        double yy = yPXB3 - dca3 * cos(phi3);

                        double vx = xx - xmid3;//from module center
                        double vy = yy - ymid3;
                        double vv = sqrt( vx*vx + vy*vy );

                        double f0 = uphi;//
                        double kx = kap*xx;
                        double ky = kap*yy;
                        double kd = kap*udca;

                        // Solve track equation for s:

                        double dx = kx - (kd-1)*sin(f0);
                        double dy = ky + (kd-1)*cos(f0);
                        double ks = atan2( dx, -dy ) - f0;// turn angle

                        // Limit to half-turn:

                        if(      ks >  pi ) ks = ks - twopi;
                        else if( ks < -pi ) ks = ks + twopi;

                        double s = ks*rho; // signed
                        double uz3 = uz0 + s*tandip; // track z at R3
                        double dz3 = zPXB3 - uz3;

                        Surface::GlobalPoint gp3( xx, yy, uz3 );
                        Surface::LocalPoint lp3 = det3->toLocal( gp3 );

                        if( idbg ) {
                            std::cout <<"**** local  Point coord ****" <<std::endl;
                            std::cout <<"gp3.x() "<< gp3.x() <<std::endl;
                            std::cout <<"gp3.y() "<< gp3.y() <<std::endl;
                            std::cout <<"gp3.z() "<< gp3.z() <<std::endl;

                            std::cout <<"lp3.x() "<< lp3.x() <<std::endl;
                            std::cout <<"lp3.y() "<< lp3.y() <<std::endl;
                            std::cout <<"lp3.z() "<< lp3.z() <<std::endl;

                            std::cout <<"  uPXB3 = xloc;  precise hit in local coordinates (w.r.t. sensor center)"<< uPXB3 <<std::endl;
                            std::cout <<"  vPXB3 = yloc;precise hit in local coordinates (w.r.t. sensor center)" << vPXB3 <<std::endl;

                        }
                        // local x = phi
                        // local y = z in barrel
                        // local z = radial in barrel (thickness)

                        //double xpix = fmod( uPXB3 + 0.82, 0.01 ); // xpix = 0..0.01 reconstructed
                        //double xpx3 = fmod( uPXB3 + 0.82, 0.02 ); // xpix = 0..0.02 reconstructed
                        //double xpx2 = fmod( uPXB2 + 0.82, 0.01 ); // xpix = 0..0.01 reconstructed
                        //double xpx4 = fmod( uPXB4 + 0.82, 0.01 ); // xpix = 0..0.01 reconstructed

                        //double dpix = fmod( uPXB3 + dca3 + 0.82, 0.01 ); // dpix = 0..0.01 predicted

                        double vpix = fmod( vv, 0.01 ); // vpix = 0..0.01 predicted
                        if( uPXB3 < 0 ) vpix = -vpix; // vv is unsigned distance from module center

                        double lpix = fmod( lp3.x() + 0.82, 0.01 ); // lpix = 0..0.01 predicted
                        //double tpix = fmod( lp3.x() + 0.82, 0.02 ); // tpix = 0..0.02 predicted

                        //double zpix = fmod( lp3.y() + 3.24, 0.015 ); // zpix = 0..0.015 predicted
                        //double spix = fmod( lp3.y() + 3.24, 0.03  ); // spix = 0..0.03  predicted

                        //int smin = zmin3%52; // 0..51 column along z
                        //int smax = zmax3%52; // 0..51 column along z

                        double cogx = (cogp3 + 0.5 - 80) * 0.01 - 0.0054; // Lorentz shift
                        if( cogp3 < 79 ) cogx -= 0.01; // big pix
                        if( cogp3 > 80 ) cogx += 0.01; // big pix

                        double mpix = fmod( cogx + 0.82, 0.01 ); // mpix = 0..0.01 from cluster COG
                        double cogdx = cogx - lp3.x(); // residual

                        layer3dx = dca3 * 1E4;
                        layer3dz = dz3 * 1E4;

                        // hybrid method:

                        //double hybx = uPXB3; // template
                        //if( mpix*1E4 < 20 ) hybx = cogx; // COG
                        //if( mpix*1E4 > 75 ) hybx = cogx;
                        //double hpix = fmod( hybx + 0.82, 0.01 ); // hpix = 0..0.01 from cluster hybrid method
                        //double hybdx = hybx - lp3.x(); // residual

                        /*bool halfmod = 0;
                          if(      ilad3 ==  8 ) halfmod = 1;
                          else if( ilad3 ==  9 ) halfmod = 1;
                          else if( ilad3 == 24 ) halfmod = 1;
                          else if( ilad3 == 25 ) halfmod = 1;
                          */

                        // residual profiles: alignment check


                    }//pt > 4, ! halfmod
                }

                //triplet 1+3 -> 2




                // triplet 1+2 -> 4:
                if( n1*n2*n4 > 0 ) {
                    {
                        if( jdbg ) cout << "  triplet 1+2 -> 4\n";

                        double f4 = atan2( yPXB4, xPXB4 );//position angle in layer 4

                        double ax = xPXB2 - xPXB1;
                        double ay = yPXB2 - yPXB1;
                        double aa = sqrt( ax*ax + ay*ay ); // from 1 to 2

                        double xmid = 0.5 * ( xPXB1 + xPXB2 );
                        double ymid = 0.5 * ( yPXB1 + yPXB2 );
                        double bx = xPXB4 - xmid;
                        double by = yPXB4 - ymid;
                        double bb = sqrt( bx*bx + by*by ); // from mid point to point 3

                        // Calculate the centre of the helix in xy-projection that
                        // transverses the two spacepoints. The points with the same
                        // distance from the two points are lying on a line.
                        // LAMBDA is the distance between the point in the middle of
                        // the two spacepoints and the centre of the helix.

                        // we already have kap and rho = 1/kap

                        double lam = sqrt( -0.25 +
                                rho*rho / ( ( xPXB1 - xPXB2 )*( xPXB1 - xPXB2 ) + ( yPXB1 - yPXB2 )*( yPXB1 - yPXB2 ) ) );

                        // There are two solutions, the sign of kap gives the information
                        // which of them is correct:

                        if( kap > 0 ) lam = -lam;

                        // ( X0, Y0 ) is the centre of the circle
                        // that describes the helix in xy-projection:

                        double x0 =  0.5*( xPXB1 + xPXB2 ) + lam * ( -yPXB1 + yPXB2 );
                        double y0 =  0.5*( yPXB1 + yPXB2 ) + lam * (  xPXB1 - xPXB2 );

                        // Calculate theta:

                        double num = ( yPXB2 - y0 ) * ( xPXB1 - x0 ) - ( xPXB2 - x0 ) * ( yPXB1 - y0 );
                        double den = ( xPXB1 - x0 ) * ( xPXB2 - x0 ) + ( yPXB1 - y0 ) * ( yPXB2 - y0 );
                        double tandip = kap * ( zPXB2 - zPXB1 ) / atan( num / den );
                        double udip = atan(tandip);
                        //double utet = pihalf - udip;

                        // To get phi0 in the right interval one must distinguish
                        // two cases with positve and negative kap:

                        double uphi;
                        if( kap > 0 ) uphi = atan2( -x0,  y0 );
                        else          uphi = atan2(  x0, -y0 );

                        // The distance of the closest approach DCA depends on the sign
                        // of kappa:

                        double udca;
                        if( kap > 0 ) udca = rho - sqrt( x0*x0 + y0*y0 );
                        else          udca = rho + sqrt( x0*x0 + y0*y0 );

                        // angle from first hit to dca point:

                        double dphi = atan( ( ( xPXB1 - x0 ) * y0 - ( yPXB1 - y0 ) * x0 )
                                / ( ( xPXB1 - x0 ) * x0 + ( yPXB1 - y0 ) * y0 ) );

                        double uz0 = zPXB1 + tandip * dphi * rho;


                        // extrapolate to outer hit:
                        // cirmov
                        // we already have rinv = -kap

                        double cosphi = cos(uphi);
                        double sinphi = sin(uphi);
                        double dp = -xPXB4*sinphi + yPXB4*cosphi + udca;
                        double dl = -xPXB4*cosphi - yPXB4*sinphi;
                        double sa = 2*dp + rinv * ( dp*dp + dl*dl );
                        double dca4 = sa / ( 1 + sqrt(1 + rinv*sa) );// distance to hit 4
                        double ud = 1 + rinv*udca;
                        double phi4 = atan2( -rinv*xPXB4 + ud*sinphi, rinv*yPXB4 + ud*cosphi );//track direction

                        double phiinc = phi4 - phiN4;//angle of incidence in rphi w.r.t. normal vector

                        // phiN alternates inward/outward
                        // reduce phiinc:

                        if( phiinc > pihalf ) phiinc -= pi;
                        else if( phiinc < -pihalf ) phiinc += pi;

                        // arc length:

                        double xx = xPXB4 + dca4 * sin(phi4); // point on track
                        double yy = yPXB4 - dca4 * cos(phi4);

                        double f0 = uphi;//
                        double kx = kap*xx;
                        double ky = kap*yy;
                        double kd = kap*udca;

                        // Solve track equation for s:

                        double dx = kx - (kd-1)*sin(f0);
                        double dy = ky + (kd-1)*cos(f0);
                        double ks = atan2( dx, -dy ) - f0;// turn angle

                        //---  Limit to half-turn:

                        if(      ks >  pi ) ks = ks - twopi;
                        else if( ks < -pi ) ks = ks + twopi;

                        double s = ks*rho;// signed
                        double uz4 = uz0 + s*tandip; //track z at R4
                        double dz4 = zPXB4 - uz4;

                    }
                    ////END PXB4
                }


                //------------------------------------------------------------------------
                // 1-2-3 pixel triplet:
                if( n1*n2*n3 > 0 ) {

                    { // let's open a scope, so we can redefine the variables further down

                        if( jdbg ) cout << "  triplet 1+3 -> 2\n";

                        double f2 = atan2( yPXB2, xPXB2 );//position angle in layer 2

                        double ax = xPXB3 - xPXB1;
                        double ay = yPXB3 - yPXB1;
                        double aa = sqrt( ax*ax + ay*ay ); // from 1 to 3

                        double xmid = 0.5 * ( xPXB1 + xPXB3 );
                        double ymid = 0.5 * ( yPXB1 + yPXB3 );
                        double bx = xPXB2 - xmid;
                        double by = yPXB2 - ymid;
                        double bb = sqrt( bx*bx + by*by ); // from mid point to point 2


                        // Author: Johannes Gassner (15.11.1996)
                        // Make track from 2 space points and kappa (cmz98/ftn/csmktr.f)
                        // Definition of the Helix :

                        // x( t ) = X0 + KAPPA^-1 * SIN( PHI0 + t )
                        // y( t ) = Y0 - KAPPA^-1 * COS( PHI0 + t )          t > 0
                        // z( t ) = Z0 + KAPPA^-1 * TAN( DIP ) * t

                        // Center of the helix in the xy-projection:

                        // X0 = + ( DCA - KAPPA^-1 ) * SIN( PHI0 )
                        // Y0 = - ( DCA - KAPPA^-1 ) * COS( PHI0 )

                        // Point 1 must be in the inner layer, 3 in the outer:

                        double r1 = sqrt( xPXB1*xPXB1 + yPXB1*yPXB1 );
                        double r3 = sqrt( xPXB3*xPXB3 + yPXB3*yPXB3 );

                        //	cout << "!!!warn r1 = " << r1 << ", r3 = " << r3 << endl;

                        if( r3-r1 < 2.0 ) cout << "warn r1 = " << r1 << ", r3 = " << r3 << endl;

                        // Calculate the centre of the helix in xy-projection that
                        // transverses the two spacepoints. The points with the same
                        // distance from the two points are lying on a line.
                        // LAMBDA is the distance between the point in the middle of
                        // the two spacepoints and the centre of the helix.

                        // we already have kap and rho = 1/kap

                        double lam = sqrt( -0.25 +
                                rho*rho / ( ( xPXB1 - xPXB3 )*( xPXB1 - xPXB3 ) + ( yPXB1 - yPXB3 )*( yPXB1 - yPXB3 ) ) );

                        // There are two solutions, the sign of kap gives the information
                        // which of them is correct:

                        if( kap > 0 ) lam = -lam;

                        // ( X0, Y0 ) is the centre of the circle
                        // that describes the helix in xy-projection:

                        double x0 =  0.5*( xPXB1 + xPXB3 ) + lam * ( -yPXB1 + yPXB3 );
                        double y0 =  0.5*( yPXB1 + yPXB3 ) + lam * (  xPXB1 - xPXB3 );

                        // Calculate theta:

                        double num = ( yPXB3 - y0 ) * ( xPXB1 - x0 ) - ( xPXB3 - x0 ) * ( yPXB1 - y0 );
                        double den = ( xPXB1 - x0 ) * ( xPXB3 - x0 ) + ( yPXB1 - y0 ) * ( yPXB3 - y0 );
                        double tandip = kap * ( zPXB3 - zPXB1 ) / atan( num / den );
                        double udip = atan(tandip);
                        //double utet = pihalf - udip;

                        // To get phi0 in the right interval one must distinguish
                        // two cases with positve and negative kap:

                        double uphi;
                        if( kap > 0 ) uphi = atan2( -x0,  y0 );
                        else          uphi = atan2(  x0, -y0 );

                        // The distance of the closest approach DCA depends on the sign
                        // of kappa:

                        double udca;
                        if( kap > 0 ) udca = rho - sqrt( x0*x0 + y0*y0 );
                        else          udca = rho + sqrt( x0*x0 + y0*y0 );

                        // angle from first hit to dca point:

                        double dphi = atan( ( ( xPXB1 - x0 ) * y0 - ( yPXB1 - y0 ) * x0 )
                                / ( ( xPXB1 - x0 ) * x0 + ( yPXB1 - y0 ) * y0 ) );

                        double uz0 = zPXB1 + tandip * dphi * rho;


                        // interpolate to middle hit:
                        // cirmov
                        // we already have rinv = -kap

                        double cosphi = cos(uphi);
                        double sinphi = sin(uphi);
                        double dp = -xPXB2*sinphi + yPXB2*cosphi + udca;
                        double dl = -xPXB2*cosphi - yPXB2*sinphi;
                        double sa = 2*dp + rinv * ( dp*dp + dl*dl );
                        double dca2 = sa / ( 1 + sqrt(1 + rinv*sa) );// distance to hit 2
                        double ud = 1 + rinv*udca;
                        double phi2 = atan2( -rinv*xPXB2 + ud*sinphi, rinv*yPXB2 + ud*cosphi );//direction

                        double phiinc = phi2 - phiN2;//angle of incidence in rphi w.r.t. normal vector

                        // phiN alternates inward/outward
                        // reduce phiinc:

                        if( phiinc > pihalf ) phiinc -= pi;
                        else if( phiinc < -pihalf ) phiinc += pi;

                        // arc length:

                        double xx = xPXB2 + dca2 * sin(phi2); // point on track
                        double yy = yPXB2 - dca2 * cos(phi2);

                        double vx = xx - xmid2;//from module center
                        double vy = yy - ymid2;
                        double vv = sqrt( vx*vx + vy*vy );

                        double f0 = uphi;//
                        double kx = kap*xx;
                        double ky = kap*yy;
                        double kd = kap*udca;

                        // Solve track equation for s:

                        double dx = kx - (kd-1)*sin(f0);
                        double dy = ky + (kd-1)*cos(f0);
                        double ks = atan2( dx, -dy ) - f0;// turn angle

                        // Limit to half-turn:

                        if(      ks >  pi ) ks = ks - twopi;
                        else if( ks < -pi ) ks = ks + twopi;

                        double s = ks*rho; // signed
                        double uz2 = uz0 + s*tandip; // track z at R2
                        double dz2 = zPXB2 - uz2;

                        Surface::GlobalPoint gp2( xx, yy, uz2 );
                        Surface::LocalPoint lp2 = det2->toLocal( gp2 );

                        if( idbg ) {
                            std::cout <<"**** local  Point coord ****" <<std::endl;
                            std::cout <<"gp2.x() "<< gp2.x() <<std::endl;
                            std::cout <<"gp2.y() "<< gp2.y() <<std::endl;
                            std::cout <<"gp2.z() "<< gp2.z() <<std::endl;

                            std::cout <<"lp2.x() "<< lp2.x() <<std::endl;
                            std::cout <<"lp2.y() "<< lp2.y() <<std::endl;
                            std::cout <<"lp2.z() "<< lp2.z() <<std::endl;

                            std::cout <<"  uPXB2 = xloc;  precise hit in local coordinates (w.r.t. sensor center)"<< uPXB2 <<std::endl;
                            std::cout <<"  vPXB2 = yloc;precise hit in local coordinates (w.r.t. sensor center)" << vPXB2 <<std::endl;

                        }
                        // local x = phi
                        // local y = z in barrel
                        // local z = radial in barrel (thickness)

                        double xpix = fmod( uPXB2 + 0.82, 0.01 ); // xpix = 0..0.01 reconstructed
                        double xpx2 = fmod( uPXB2 + 0.82, 0.02 ); // xpix = 0..0.02 reconstructed
                        double xpx1 = fmod( uPXB1 + 0.82, 0.01 ); // xpix = 0..0.01 reconstructed
                        double xpx3 = fmod( uPXB3 + 0.82, 0.01 ); // xpix = 0..0.01 reconstructed

                        //double dpix = fmod( uPXB2 + dca2 + 0.82, 0.01 ); // dpix = 0..0.01 predicted

                        double vpix = fmod( vv, 0.01 ); // vpix = 0..0.01 predicted
                        if( uPXB2 < 0 ) vpix = -vpix; // vv is unsigned distance from module center

                        double lpix = fmod( lp2.x() + 0.82, 0.01 ); // lpix = 0..0.01 predicted
                        double tpix = fmod( lp2.x() + 0.82, 0.02 ); // tpix = 0..0.02 predicted

                        double zpix = fmod( lp2.y() + 3.24, 0.015 ); // zpix = 0..0.015 predicted
                        double spix = fmod( lp2.y() + 3.24, 0.03  ); // spix = 0..0.03  predicted

                        int smin = zmin2%52; // 0..51 column along z
                        int smax = zmax2%52; // 0..51 column along z

                        double cogx = (cogp2 + 0.5 - 80) * 0.01 - 0.0054; // Lorentz shift
                        if( cogp2 < 79 ) cogx -= 0.01; // big pix
                        if( cogp2 > 80 ) cogx += 0.01; // big pix

                        double mpix = fmod( cogx + 0.82, 0.01 ); // mpix = 0..0.01 from cluster COG
                        double cogdx = cogx - lp2.x(); // residual

                        // hybrid method:

                        double hybx = uPXB2; // template
                        if( mpix*1E4 < 20 ) hybx = cogx; // COG
                        if( mpix*1E4 > 75 ) hybx = cogx;
                        //double hpix = fmod( hybx + 0.82, 0.01 ); // hpix = 0..0.01 from cluster hybrid method
                        double hybdx = hybx - lp2.x(); // residual

                        // bool halfmod = 0;
                        // if(      ilad2 ==  8 ) halfmod = 1;
                        // else if( ilad2 ==  9 ) halfmod = 1;
                        // else if( ilad2 == 24 ) halfmod = 1;
                        // else if( ilad2 == 25 ) halfmod = 1;
                        layer2dx = dca2*1E4;
                        layer2dz = dz2*1E4;


                    }//pt > 4, ! halfmod

                }//triplet 1+3 -> 2

                //------------------------------------------------------------------------
                // Triplet 1+3 -> 2 with refitted track:


                //------------------------------------------------------------------------
                // triplet 3+2 -> 1:

                {
                    if( jdbg ) cout << "  triplet 3+2 -> 1\n";

                    double f1 = atan2( yPXB1, xPXB1 );//position angle in layer 1

                    double ax = xPXB3 - xPXB2;
                    double ay = yPXB3 - yPXB2;
                    double aa = sqrt( ax*ax + ay*ay ); // from 2 to 3

                    double xmid = 0.5 * ( xPXB2 + xPXB3 );
                    double ymid = 0.5 * ( yPXB2 + yPXB3 );
                    double bx = xPXB1 - xmid;
                    double by = yPXB1 - ymid;
                    double bb = sqrt( bx*bx + by*by ); // from mid point to point 1

                    // Calculate the centre of the helix in xy-projection that
                    // transverses the two spacepoints. The points with the same
                    // distance from the two points are lying on a line.
                    // LAMBDA is the distance between the point in the middle of
                    // the two spacepoints and the centre of the helix.

                    // we already have kap and rho = 1/kap

                    double lam = sqrt( -0.25 +
                            rho*rho / ( ( xPXB2 - xPXB3 )*( xPXB2 - xPXB3 ) + ( yPXB2 - yPXB3 )*( yPXB2 - yPXB3 ) ) );

                    // There are two solutions, the sign of kap gives the information
                    // which of them is correct:

                    if( kap > 0 ) lam = -lam;

                    // ( X0, Y0 ) is the centre of the circle
                    // that describes the helix in xy-projection:

                    double x0 =  0.5*( xPXB2 + xPXB3 ) + lam * ( -yPXB2 + yPXB3 );
                    double y0 =  0.5*( yPXB2 + yPXB3 ) + lam * (  xPXB2 - xPXB3 );

                    // Calculate theta:

                    double num = ( yPXB3 - y0 ) * ( xPXB2 - x0 ) - ( xPXB3 - x0 ) * ( yPXB2 - y0 );
                    double den = ( xPXB2 - x0 ) * ( xPXB3 - x0 ) + ( yPXB2 - y0 ) * ( yPXB3 - y0 );
                    double tandip = kap * ( zPXB3 - zPXB2 ) / atan( num / den );
                    double udip = atan(tandip);
                    //double utet = pihalf - udip;

                    // To get phi0 in the right interval one must distinguish
                    // two cases with positve and negative kap:

                    double uphi;
                    if( kap > 0 ) uphi = atan2( -x0,  y0 );
                    else          uphi = atan2(  x0, -y0 );

                    // The distance of the closest approach DCA depends on the sign
                    // of kappa:

                    double udca;
                    if( kap > 0 ) udca = rho - sqrt( x0*x0 + y0*y0 );
                    else          udca = rho + sqrt( x0*x0 + y0*y0 );

                    // angle from first hit to dca point:

                    double dphi = atan( ( ( xPXB2 - x0 ) * y0 - ( yPXB2 - y0 ) * x0 )
                            / ( ( xPXB2 - x0 ) * x0 + ( yPXB2 - y0 ) * y0 ) );

                    double uz0 = zPXB2 + tandip * dphi * rho;

                    // extrapolate to inner hit:
                    // cirmov
                    // we already have rinv = -kap

                    double cosphi = cos(uphi);
                    double sinphi = sin(uphi);
                    double dp = -xPXB1*sinphi + yPXB1*cosphi + udca;
                    double dl = -xPXB1*cosphi - yPXB1*sinphi;
                    double sa = 2*dp + rinv * ( dp*dp + dl*dl );
                    double dca1 = sa / ( 1 + sqrt(1 + rinv*sa) );// distance to hit 1
                    double ud = 1 + rinv*udca;
                    double phi1 = atan2( -rinv*xPXB1 + ud*sinphi, rinv*yPXB1 + ud*cosphi );//direction

                    double phiinc = phi1 - phiN1;//angle of incidence in rphi w.r.t. normal vector

                    // phiN alternates inward/outward
                    // reduce phiinc:

                    if( phiinc > pihalf ) phiinc -= pi;
                    else if( phiinc < -pihalf ) phiinc += pi;

                    // arc length:

                    double xx = xPXB1 + dca1 * sin(phi1); // point on track
                    double yy = yPXB1 - dca1 * cos(phi1);

                    double f0 = uphi;//
                    double kx = kap*xx;
                    double ky = kap*yy;
                    double kd = kap*udca;

                    // Solve track equation for s:

                    double dx = kx - (kd-1)*sin(f0);
                    double dy = ky + (kd-1)*cos(f0);
                    double ks = atan2( dx, -dy ) - f0;// turn angle

                    //---  Limit to half-turn:

                    if(      ks >  pi ) ks = ks - twopi;
                    else if( ks < -pi ) ks = ks + twopi;

                    double s = ks*rho;// signed
                    double uz1 = uz0 + s*tandip; //track z at R1
                    double dz1 = zPXB1 - uz1;

                    layer1dx = dca1*1E4;
                    layer1dz = dz1*1E4;

                }






            }//3 PXB layers

            //------------------------------------------------------------------------






            tree->Fill();

        }// loop over tracks

    }
}//event

//----------------------------------------------------------------------
// method called just after ending the event loop:
//
void Pixel::endJob() {

    std::cout << "end of job after " << myCounters::neve << " events.\n";

}

//define this as a plug-in
DEFINE_FWK_MODULE(Pixel);
