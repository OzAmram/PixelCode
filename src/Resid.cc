/*
// CMS and user include files:
#include "data_test/pixel/src/Resid.h"

Resid::Resid(const edm::ParameterSet& iConfig// , edm::ConsumesCollector && ic
        )
{
    std::cout << "Resid constructed\n";
    _triggerSrc = iConfig.getParameter<edm::InputTag>("triggerSource");
    _ttrhBuilder = iConfig.getParameter<std::string>("ttrhBuilder");
    std::cout<<_triggerSrc<<" "<<_triggerSrc.label()<<" "<<_triggerSrc.process()<<" "
        <<_triggerSrc.instance()<<" "<<std::endl;



    //Definition of parameters
    t_triggerSrc_ = consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggerSource"));
    t_offlineBeamSpot_ =    consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
    t_offlinePrimaryVertices_ =   consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
    t_generalTracks_= consumes<reco::TrackCollection> (edm::InputTag("generalTracks"));//"generalTracks"));
    _OC_beginning=iConfig.getParameter<int>("orbit_beginning");
    _OC_end=iConfig.getParameter<int>("orbit_end");
    // tok_caloHH_ = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits", "HcalHits"));
    edm::Service<TFileService> fsT;
    tree = fsT->make<TTree>("tree", "tree");
    tree->Branch("trackIntersectx", &trackIntersectx);
    tree->Branch("trackIntersecty", &trackIntersecty);
    tree->Branch("resid1Dx", &resid1Dx);
    tree->Branch("resid1Dy", &resid1Dy);
    tree->Branch("resid2Dx", &resid2Dx);
    tree->Branch("resid2Dy", &resid2Dy);
    tree->Branch("trkEta", &trkEta);
    tree->Branch("trkPt", &trkPt);

    tree->Branch("clustSizeX", &clustSizeX);
    tree->Branch("clustSizeY", &clustSizeY);
    tree->Branch("clustxmax", &clustxmax);
    tree->Branch("clustxmin", &clustxmin);
    tree->Branch("clustymax", &clustymax);
    tree->Branch("clustymin", &clustymin);

    tree->Branch("probQ1D", &probQ1D);
    tree->Branch("probXY1D", &probXY1D);
    tree->Branch("qBin1D", &qBin1D);
    tree->Branch("hasBadPixels1D", &hasBadPixels1D);
    tree->Branch("isOnEdge1D", &isOnEdge1D);
    tree->Branch("hasFilledProb1D", &hasFilledProb1D);

    tree->Branch("probQ2D", &probQ2D);
    tree->Branch("probXY2D", &probXY2D);
    tree->Branch("qBin2D", &qBin2D);
    tree->Branch("hasBadPixels2D", &hasBadPixels2D);
    tree->Branch("isOnEdge2D", &isOnEdge2D);
    tree->Branch("hasFilledProb2D", &hasFilledProb2D);

    tree->Branch("hit_type", &hit_type);
    tree->Branch("onTrack", &onTrack);

    // Tokens
    rawDataErrorToken_ = consumes<edm::DetSetVector<SiPixelRawDataError>> (edm::InputTag("siPixelDigis"));
    primaryVerticesToken_  = consumes<reco::VertexCollection> (edm::InputTag("offlinePrimaryVertices"));
    triggerResultsToken_ = consumes<edm::TriggerResults>(triggerTag_);
    pileupSummaryToken_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"));
    clustersToken_ = consumes<edmNew::DetSetVector<SiPixelCluster>>(edm::InputTag("siPixelClusters"));
    trajTrackCollectionToken_ = consumes<TrajTrackAssociationCollection>
        (iConfig.getParameter<edm::InputTag>("trajectoryInput"));
    measurementTrackerEventToken_ = consumes<MeasurementTrackerEvent>
        (edm::InputTag("MeasurementTrackerEvent"));
    pixelDigiCollectionToken_ = consumes<edm::DetSetVector<PixelDigi>>
        (edm::InputTag("simSiPixelDigis"));

}
//
// destructor:
//
Resid::~Resid()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions:
// method called once each job just before starting event loop
//
void Resid::beginJob()
{
}

//----------------------------------------------------------------------
// method called for each event:

void Resid::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
    const int run = iRun.run();


}

// trajectory measurement: contains all info about measurement of a trajectory by a DetId(detector number)                                   
// - tracking rechit                                                                                                                         
// - predicted TrajectoryStateOnSurface (fitter and smoother)                                                                                
// - combination of TrajectoryStateOnSurfaces                                                                                                
// - compatibility estimate between tracking rechit and predicted state                                                                      

// find hits nearby                                                                                                                          
// given trajectory measurement                                                                                                              
// given clustercollection                                                                                                                   
// fiven trajtrackassociationcollection  
void Resid::checkAndSaveTrajMeasurementData
( const TrajectoryMeasurement& measurement,
  const edm::Handle<edmNew::DetSetVector<SiPixelCluster>>& clusterCollectionHandle,
  const edm::Handle<TrajTrackAssociationCollection>& trajTrackCollectionHandle,
  edm::OwnVector<TrackingRecHit> &layer1_hits
  ) {

    // Check if the measurement infos can be read                                                                                              
    if(!measurement.updatedState().isValid()) return;

    // given a measurement get tracking rechit                                                                                                 
    // for that detId that should be on pixel                                                                                                  
    TransientTrackingRecHit::ConstRecHitPointer recHit = measurement.recHit();
    DetId detId = recHit -> geographicalId();
    if(!ResidHelpers::detidIsOnPixel(detId)) return;

    // get estimated trajectory state at a surface close to where you want hits                                                                
    TrajectoryStateOnSurface trajStateOnSurface = ResidHelpers::getTrajectoryStateOnSurface(measurement);

    // skip hits with undeterminable positions                                                                                                 
    if(!(trajStateOnSurface.isValid())) return;

    // from estimated trajectory state get global and local position (x,y,z,validhit,roc)                                                      
    //GlobalPoint globalPosition     = trajStateOnSurface.globalPosition();
    LocalPoint  localPosition      = trajStateOnSurface.localPosition();
    //LocalError  localPositionError = trajStateOnSurface.localError().positionError();

    // get track local parameters: momentum,direction,angles                                                                                   
    LocalTrajectoryParameters trajectoryParameters = trajStateOnSurface.localParameters();
    //auto trajectoryMomentum = trajectoryParameters.momentum();
    //LocalVector localTrackDirection = trajectoryMomentum / trajectoryMomentum.mag();

    // from tracking rechit: check if is valid(what does that mean?) and a hit                                                                 
    // if yes: => save params                                                                                                                  
    // else: => loop over siPixelClusters on detector and get closest cluster to track local position                                          
    const SiPixelCluster* clust = nullptr;

    if(recHit -> isValid() && recHit -> hit() != 0) {
        const SiPixelRecHit *hit = static_cast<const SiPixelRecHit*>(recHit -> hit());
        hit_type = 1;
        clust = hit -> cluster().get();
        //float currentMinValueSquared = clusterPointDistanceSquared(detId, *clust_1D, localPosition, trajectoryParameters, false);
        //printf("recHit on trajectory measure (?), distance to traj %f and layer %i \n ",currentMinValueSquared,tTopo->pxbLayer(detId));
    } else {
        if(clusterCollectionHandle.isValid()) {
            const edmNew::DetSetVector<SiPixelCluster>::const_iterator clustersOnDet =
                clusterCollectionHandle -> find(detId);
            if(clustersOnDet != clusterCollectionHandle -> end()){
                hit_type = 2;
                clust = getClosestClusterOnDetSetToPoint(*clustersOnDet, localPosition, trajectoryParameters, false);
            }
        }
    }
    if(clust == nullptr){
        printf("Null ptrs \n");
        return;
    }
    clustSizeX = clust->sizeX();
    clustSizeY = clust->sizeY();
    clustxmin = clust->minPixelRow();
    clustxmax = clust->maxPixelRow();
    clustymin = clust->minPixelCol();
    clustymax = clust->maxPixelCol();

    onTrack = false;
    for(auto l1_rhit = layer1_hits.begin(); l1_rhit != layer1_hits.end(); l1_rhit++ ){
        const SiPixelRecHit *si_hit = static_cast<const SiPixelRecHit*>(l1_rhit -> hit());
        auto l1_clust = si_hit -> cluster().get();
        if((l1_clust->sizeX() == clustSizeX) && 
           (l1_clust->sizeY() == clustSizeY) && 
           (l1_clust->minPixelRow() == clustxmin) && 
           (l1_clust->minPixelCol() == clustymin)){
                onTrack = true;
                printf("onTrack = true \n");
        }


    }



    LocalPoint clust1DCoords, clust2DCoords;
    SiPixelRecHitQuality::QualWordType QW_1D, QW_2D;
    const GeomDetUnit* geomDetUnit = trackerGeometry_ -> idToDetUnit(detId);
    std::tie(clust1DCoords, std::ignore, QW_1D) = pixelClusterParameterEstimator1D_ -> getParameters(*clust, *geomDetUnit, trajectoryParameters);
    std::tie(clust2DCoords, std::ignore, QW_2D) = pixelClusterParameterEstimator2D_ -> getParameters(*clust, *geomDetUnit, trajectoryParameters);
    trackIntersectx = localPosition.x();                                                                                                   
    trackIntersecty = localPosition.y();                                                                                              
    // if cluster exists: get cluster parameters: charge,size, local and global position, and pixels charge and position                       

    resid1Dx = 1E4* (clust1DCoords.x() - trackIntersectx);
    resid1Dy = 1E4 *(clust1DCoords.y() - trackIntersecty);


    resid2Dx = 1E4* (clust2DCoords.x() - trackIntersectx);
    resid2Dy = 1E4 *(clust2DCoords.y() - trackIntersecty);

    probQ1D = recHitQuality.thePacking.probabilityQ(QW_1D);
    probXY1D = recHitQuality.thePacking.probabilityXY(QW_1D);
    qBin1D = recHitQuality.thePacking.qBin(QW_1D);
    hasBadPixels1D = recHitQuality.thePacking.hasBadPixels(QW_1D);
    isOnEdge1D = recHitQuality.thePacking.isOnEdge(QW_1D);
    hasFilledProb1D = recHitQuality.thePacking.isOnEdge(QW_1D);

    probQ2D = recHitQuality.thePacking.probabilityQ(QW_2D);
    probXY2D = recHitQuality.thePacking.probabilityXY(QW_2D);
    qBin2D = recHitQuality.thePacking.qBin(QW_2D);
    hasBadPixels2D = recHitQuality.thePacking.hasBadPixels(QW_2D);
    isOnEdge2D = recHitQuality.thePacking.isOnEdge(QW_2D);
    hasFilledProb2D = recHitQuality.thePacking.hasFilledProb(QW_2D);

    tree->Fill();
    printf("Valid Extrap hit \n");
    

    // Get closest other traj measurement                                                                                                      
    //getClosestOtherTrajMeasurementDistanceByLooping(measurement, trajTrackCollectionHandle, traj_.d_tr, traj_.dx_tr, traj_.dy_tr);
    //traj_.hit_near = (traj_.d_tr < 0.5); // 5 mm                                                                                               
    //traj_.clust_near = (traj_.d_cl != NOVAL_F && traj_.d_cl < HIT_CLUST_NEAR_CUT_VAL);
}

// get closest cluster to point given clusters on detector                                                                                   
    const SiPixelCluster* Resid::getClosestClusterOnDetSetToPoint
(const edmNew::DetSet<SiPixelCluster>& clustersOnDet, const LocalPoint& referencePoint, const LocalTrajectoryParameters & ltp, bool use2D)
{
    if(clustersOnDet.empty()) {
        printf("clusters not on det \n");
	return nullptr;
    }

    const DetId detId = clustersOnDet.id();
    const SiPixelCluster* minDistanceCluster = clustersOnDet.begin();
    float currentMinValueSquared = clusterPointDistanceSquared(detId, *minDistanceCluster, referencePoint, ltp, use2D);

    for(const auto& cluster: clustersOnDet) {
        float currentDistanceSquared = clusterPointDistanceSquared(detId, cluster, referencePoint, ltp, use2D);
        if(currentDistanceSquared < currentMinValueSquared) {
            currentMinValueSquared = std::move(currentDistanceSquared);
            minDistanceCluster = &cluster;
        }
    }
    printf("closest cluster on det: distance %f \n",currentMinValueSquared);
    return minDistanceCluster;
}

// given cluster and reference point calculate position from CPE and calculate distance to ref. point                                        
float Resid::clusterPointDistanceSquared(const DetId& detId, const SiPixelCluster& cluster, const LocalPoint& referencePoint,\
        const LocalTrajectoryParameters & ltp, bool use2D){
    const GeomDetUnit* geomDetUnit = trackerGeometry_ -> idToDetUnit(detId);
    LocalPoint clustLocalCoordinates;
    if (use2D) {
        std::tie(clustLocalCoordinates, std::ignore, std::ignore) = pixelClusterParameterEstimator2D_ -> getParameters(cluster, *geomDetUnit, ltp);
    }
    else{
        std::tie(clustLocalCoordinates, std::ignore, std::ignore) = pixelClusterParameterEstimator1D_ -> getParameters(cluster, *geomDetUnit, ltp);
    }

    float xDist = clustLocalCoordinates.x() - referencePoint.x();
    float yDist = clustLocalCoordinates.y() - referencePoint.y();
    float zDist = clustLocalCoordinates.z() - referencePoint.z();

    return xDist * xDist + yDist * yDist + zDist * zDist;
}



std::vector<TrajectoryMeasurement> Resid::getLayer1ExtrapolatedHitsFromMeas(const TrajectoryMeasurement& trajMeasurement)
{

    // Last layer 2 or disk 1 mesurement is to be propagated to layer 1 if possible
    // Only propagating valid measurements
    std::unique_ptr<LayerMeasurements> layerMeasurements
        (new LayerMeasurements(*measurementTracker_, *measurementTrackerEvent_));

    const DetLayer* pixelBarrelLayer1 = 
        measurementTracker_ -> geometricSearchTracker() -> pixelBarrelLayers().front();

    return layerMeasurements -> measurements(*pixelBarrelLayer1, trajMeasurement.updatedState(),
            *trackerPropagator_, *chi2MeasurementEstimator_);

}
void Resid::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup // , edm::ConsumesCollector && ic
        ){
    using namespace std;
    using namespace edm;
    using namespace reco;
    using namespace math;

    const double pi = 4*atan(1);
    const double wt = 180/pi;
    const double twopi = 2*pi;
    const double pihalf = 2*atan(1);
    //const double sqrtpihalf = sqrt(pihalf);

    //FOR NOW only do 1/3 events to get jobs to finish



    // Digis
    edm::Handle<edm::DetSetVector<PixelDigi>> digiCollectionHandle;
    if(saveDigiTree_ || npixFromDigiCollection_) iEvent.getByToken(pixelDigiCollectionToken_, digiCollectionHandle);

    // Get vertices
    edm::Handle<reco::VertexCollection>      vertexCollectionHandle;
    iEvent.getByToken(primaryVerticesToken_, vertexCollectionHandle);

    // Get trigger info
    edm::Handle<edm::TriggerResults> triggerResultsHandle;
    iEvent.getByToken(triggerResultsToken_, triggerResultsHandle);

    // Get pileup info
    edm::Handle<std::vector<PileupSummaryInfo>> puInfoCollectionHandle;
    iEvent.getByToken(pileupSummaryToken_,      puInfoCollectionHandle);

    // Get cluster collection
    edm::Handle<edmNew::DetSetVector<SiPixelCluster>> clusterCollectionHandle;
    iEvent.getByToken(clustersToken_,                 clusterCollectionHandle);

    // Get Traj-Track Collection
    edm::Handle<TrajTrackAssociationCollection>  trajTrackCollectionHandle;
    iEvent.getByToken(trajTrackCollectionToken_, trajTrackCollectionHandle);

    // TrackerTopology for module informations
    edm::ESHandle<TrackerTopology> trackerTopologyHandle;
    iSetup.get<TrackerTopologyRcd>().get(trackerTopologyHandle);
    trackerTopology_ = trackerTopologyHandle.product();

    // TrackerGeometry for module informations
    edm::ESHandle<TrackerGeometry> trackerGeometryHandle;
    iSetup.get<TrackerDigiGeometryRecord>().get(trackerGeometryHandle);
    trackerGeometry_ = trackerGeometryHandle.product();

    // Tracker propagator for propagating tracks to other layers
    edm::ESHandle<Propagator> propagatorHandle;
    iSetup.get<TrackingComponentsRecord>().get("PropagatorWithMaterial", propagatorHandle);
    std::unique_ptr<Propagator> propagatorUniquePtr(propagatorHandle.product() -> clone());
    trackerPropagator_ = propagatorUniquePtr.get();
    const_cast<Propagator*>(trackerPropagator_) -> setPropagationDirection(oppositeToMomentum);

    // Measurement Tracker Handle
    edm::ESHandle<MeasurementTracker> measurementTrackerHandle;
    iSetup.get<CkfComponentsRecord>().get(measurementTrackerHandle);
    measurementTracker_ = measurementTrackerHandle.product();

    // Measurement Tracker event
    edm::Handle<MeasurementTrackerEvent> measurementTrackerEventHandle;
    iEvent.getByToken(measurementTrackerEventToken_, measurementTrackerEventHandle);
    measurementTrackerEvent_ = measurementTrackerEventHandle.product();

    // Measurement estimator
    edm::ESHandle<Chi2MeasurementEstimatorBase> chi2MeasurementEstimatorHandle;
    iSetup.get<TrackingComponentsRecord>().get("Chi2", chi2MeasurementEstimatorHandle);
    chi2MeasurementEstimator_ = chi2MeasurementEstimatorHandle.product();

    // Pixel Parameter estimator 1D Reco                                                                                                     
    edm::ESHandle<PixelClusterParameterEstimator> pixelClusterParameterEstimatorHandle1D;
    iSetup.get<TkPixelCPERecord>().get("PixelCPETemplateReco", pixelClusterParameterEstimatorHandle1D);
    pixelClusterParameterEstimator1D_ = pixelClusterParameterEstimatorHandle1D.product();

    // Pixel Parameter estimator 2D Reco                                                                                                     
    edm::ESHandle<PixelClusterParameterEstimator> pixelClusterParameterEstimatorHandle2D;
    iSetup.get<TkPixelCPERecord>().get("PixelCPEClusterRepair", pixelClusterParameterEstimatorHandle2D);
    pixelClusterParameterEstimator2D_ = pixelClusterParameterEstimatorHandle2D.product();

    coord_.init(iSetup);


    int idbg = 0;  // printout for the first few events


    if( idbg ) {
        cout << endl;
        cout << "run " << iEvent.run();
        cout << ", LumiBlock " << iEvent.luminosityBlock();
        cout << ", event " << iEvent.eventAuxiliary().event();
        time_t unixZeit = iEvent.time().unixTime();
        cout << ", taken " << ctime(&unixZeit); // ctime has endline
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
    tTopo=tTopoH.product();

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
    if(maxSumPt < 1 ) return;

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
        if(pt < 10.) continue;

        //if( pt < 0.75 ) continue;// curls up
        //if( pt < 1.75 ) continue;// want sharper image

        //float tmp = abs(iTrack->dxy(vtxP))/iTrack->dxyError();
        //cout<<pt<<" "<<abs(iTrack->dxy(vtxP))<<" "<<iTrack->dxyError()<<" "<<tmp<<endl;


        if( (abs( iTrack->dxy(vtxP) ) > 5*iTrack->dxyError()) ) continue; // not prompt


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


        //--------------------------------------------------------------------------
        // loop over tracker detectors:


        //--------------------------------------------------------------------------
        // rec hits from track extra:

        if( iTrack->extra().isNull() ) continue;//next track
        if( ! iTrack->extra().isAvailable() ) continue;//next track


        //uint32_t outerDetId = 0;
        //float rmax = 0.;
        //
        uint32_t innerDetId = 0;
        uint32_t layer2_detid = 0;
        float rmin = 0.;


        // std::cout << " imod 0 "<< imod<< std::endl;

        edm::OwnVector<TrackingRecHit> recHitVector, recHitVectorNo1; // for seed
        bool found_layer2 = false;

        edm::OwnVector<TrackingRecHit> layer1_hits;

        TrackingRecHit *layer2_hit;


        Trajectory::RecHitContainer coTTRHvec, coTTRHvecNo1; // for fit, constant
        TrajectoryStateOnSurface initialTSOS = tTrack.innermostMeasurementState();

        // loop over recHits on this track:

        //TrajectoryStateOnSurface initialTSOS = tTrack.outermostMeasurementState();


        for( trackingRecHit_iterator irecHit = iTrack->recHitsBegin();
                irecHit != iTrack->recHitsEnd(); ++irecHit ) {

            DetId detId = (*irecHit)->geographicalId();
            uint32_t subDet = detId.subdetId();

            // enum Detector { Tracker=1, Muon=2, Ecal=3, Hcal=4, Calo=5 };
            if( detId.det() != 1 ){
                cout << "rec hit ID = " << detId.det() << " not in tracker!?!?\n";
                continue;
            }


            int ilay = tTopo->pxbLayer(detId); //PXBDetId(detId).layer();
            if (ilay == 1){
                printf("lay 1 hit \n");
                layer1_hits.push_back((*irecHit)->clone());
                continue;
            }
            recHitVector.push_back( (*irecHit)->clone() );


            // build transient hit:
            //#ifdef NEW_TRACKINGRECHITS
            auto tmprh =
                (*irecHit)->cloneForFit(*builder->geometry()->idToDet((**irecHit).geographicalId()));
            auto transRecHit =
                hitCloner.makeShared(tmprh, initialTSOS);
            //#else

            coTTRHvec.push_back( transRecHit );

            if( ! (*irecHit)->isValid() ) continue;

            double gX = transRecHit->globalPosition().x();
            double gY = transRecHit->globalPosition().y();
            double gZ = transRecHit->globalPosition().z();
            float gR = sqrt( gX*gX + gY*gY );
            if( gR < rmin ) {
                rmin = gR;
                innerDetId = detId.rawId();
            }




        }



        PTrajectoryStateOnDet PTraj = trajectoryStateTransform::persistentState( initialTSOS, innerDetId );
        const TrajectorySeed seed( PTraj, recHitVector, alongMomentum );

        //if( idbg ) cout << "  have seed\n";
        


        Trajectory refitTrajectory = theFitter->fitOne( seed, coTTRHvec, initialTSOS );
        if(refitTrajectory.isValid()){
            //if( refitTrajectoryCollection.size() > 0 ) { // should be either 0 or 1

                //const Trajectory& refitTrajectory = refitTrajectoryCollection.front();

                auto trajectoryMeasurements = refitTrajectory.measurements();
                // Trajectory.measurements:
                auto layer1TrajectoryMeasurements = getLayer1Intercept(trajectoryMeasurements);
                if (layer1TrajectoryMeasurements.size() == 0){
                    //dummy vec
                    continue;
                }
                if(layer1TrajectoryMeasurements.size() > 2) printf("More than 2 measures \n");


                //CRIS's FUNCTION GOES HERE!
                // should pass tree? find out where to end loop
                checkAndSaveTrajMeasurementData(layer1TrajectoryMeasurements.front(), clusterCollectionHandle,trajTrackCollectionHandle, layer1_hits);

            }//refitted trajectory
            else{
                printf("refit failed \n");
            }






    }// loop over tracks
    printf("loop over \n");

}
void Resid::getModuleData(ModuleData &mod, bool online, const DetId &detId)
{

    mod.init();

    mod.det  = detId.subdetId() - 1;
    mod.shl  = coord_.quadrant(detId);
    mod.side = coord_.side(detId);

    if(detId.subdetId() == PixelSubdetector::PixelBarrel) {

        mod.sec     = coord_.sector(detId);
        mod.half    = coord_.half(detId);
        mod.layer   = coord_.layer(detId);
        mod.flipped = coord_.flipped(detId); // opposite of outer
        if(online) {
            mod.ladder = coord_.signed_ladder(detId);
            mod.module = coord_.signed_module(detId);
        } else {
            mod.ladder = coord_.ladder(detId);
            mod.module = coord_.module(detId);
        }

    } else if(detId.subdetId() == PixelSubdetector::PixelEndcap) {

        mod.ring   = coord_.ring(detId);
        mod.panel  = coord_.panel(detId);
        mod.module = coord_.module(detId);
        if(online) {
            mod.disk  = coord_.signed_disk(detId);
            mod.blade = coord_.signed_blade(detId);
        } else {
            mod.disk  = coord_.disk(detId);
            mod.blade = coord_.blade(detId);
        }

    }

    mod.rawid = detId.rawId();
    mod.fedid = coord_.fedid(detId);

    // FED error
    std::map<uint32_t, int>::const_iterator federrors_it =
        federrors_.find(detId.rawId());
    mod.federr = (federrors_it != federrors_.end()) ? federrors_it->second : 0;

}

//event
std::vector<TrajectoryMeasurement> Resid::getLayer1Intercept (std::vector<TrajectoryMeasurement> &trajectoryMeasurements) 
{


    // std::cout << "Searching for the first layer 1 trajectory measurements... ";
    // Not for the faint hearted ...
    auto firstLayer2TrajMeasurementIt = 
        std::find_if(trajectoryMeasurements.begin(), 
                trajectoryMeasurements.end(),
                [&] (const TrajectoryMeasurement& measurement) {
                ModuleData mod;
                if (!(measurement.recHit() ->isValid())) return false;
                getModuleData(mod, 1, measurement.recHit() -> geographicalId());
                return mod.det == 0 && mod.layer == 2;
                });
    // Check if the last non-layer1 traj measurement is valid 
    //int idx = firstLayer2TrajMeasurementIt - trajectoryMeasurements.begin();
    //auto firstLayer2TrajMeasurementRecHit = trajectoryMeasurements.at(idx).recHit();
    
    std::vector<TrajectoryMeasurement> dummy_vec;
    if(firstLayer2TrajMeasurementIt == trajectoryMeasurements.end()){
            printf("can't find layer 2 hit \n");
            return dummy_vec;
    }
    auto firstLayer2TrajMeasurementRecHit = firstLayer2TrajMeasurementIt -> recHit();
    if(firstLayer2TrajMeasurementRecHit == nullptr){
        std::cout << "Invalid rechit pointer." << std::endl;
        return dummy_vec;
    }
    if(!(firstLayer2TrajMeasurementRecHit -> isValid())){
        printf("HIT NOT VALID!\n");
        return dummy_vec;
    }
    //std::cout << "***" << std::endl;
    //std::cout << "Number of traj. measurements in track: "
    //    << std::distance(trajectoryMeasurements.begin(), trajectoryMeasurements.end()) 
    //    << std::endl;
    //std::cout << "First layer 2 hit index:               " 
    //    << std::distance(trajectoryMeasurements.begin(), firstLayer2TrajMeasurementIt)
    //    << std::endl;
    //std::cout << "***" << std::endl;



    ModuleData mod;
    getModuleData(mod, 1, firstLayer2TrajMeasurementRecHit -> geographicalId());
    if(mod.layer != 2){
        printf("Couldn't find layer 2 hit \n");
        return dummy_vec;
    }


    std::vector<TrajectoryMeasurement> extrapolatedHitsOnLayer1
        (getLayer1ExtrapolatedHitsFromMeas(*firstLayer2TrajMeasurementIt));

    printf("Size of extrap hits is %i \n", (int) extrapolatedHitsOnLayer1.size());

    return extrapolatedHitsOnLayer1;
}

void Resid::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
void Resid::beginLuminosityBlock(edm::LuminosityBlock const& iLumi,
        edm::EventSetup const& iSetup) {}
void Resid::endLuminosityBlock(edm::LuminosityBlock const& iLumi,
        edm::EventSetup const& iSetup) {}
//----------------------------------------------------------------------
// method called just after ending the event loop:
//
void Resid::endJob() {


}

namespace ResidHelpers
{
    bool detidIsOnPixel(const DetId& detid) {
        if (detid.det()!=DetId::Tracker) return false;
        if (detid.subdetId() == PixelSubdetector::PixelBarrel) return true;
        if (detid.subdetId() == PixelSubdetector::PixelEndcap) return true;
        return false;
    }

    TrajectoryStateOnSurface getTrajectoryStateOnSurface(const TrajectoryMeasurement& measurement) {

        static TrajectoryStateCombiner trajStateCombiner;

        const auto& forwardPredictedState  = measurement.forwardPredictedState();
        const auto& backwardPredictedState = measurement.backwardPredictedState();

        if(forwardPredictedState.isValid() && backwardPredictedState.isValid())
            return trajStateCombiner(forwardPredictedState, backwardPredictedState);

        else if(backwardPredictedState.isValid())
            return backwardPredictedState;

        else if(forwardPredictedState.isValid())
            return forwardPredictedState;

        std::cout << "Error saving traj. measurement data." \
            " Trajectory state on surface cannot be determined." << std::endl;

        return TrajectoryStateOnSurface();

    }

    std::pair<float, float> getLocalXY(const TrajectoryMeasurement& measurement) {

        std::pair<float, float> returnValue;
        TrajectoryStateOnSurface trajStateOnSurface = getTrajectoryStateOnSurface(measurement);

        if(!(trajStateOnSurface.isValid())) return std::make_pair<float, float>(NOVAL_F, NOVAL_F);

        LocalPoint localPosition = trajStateOnSurface.localPosition();
        returnValue.first  = localPosition.x();
        returnValue.second = localPosition.y();

        return returnValue;

    }

    float trajMeasurementDistanceSquared
        (const TrajectoryMeasurement& lhs, const TrajectoryMeasurement& rhs) {

            std::pair<float, float> lhsLocalXY = getLocalXY(lhs);
            std::pair<float, float> rhsLocalXY = getLocalXY(rhs);

            float dxHit = lhsLocalXY.first  - rhsLocalXY.first;
            float dyHit = lhsLocalXY.second - rhsLocalXY.second;

            float distanceSquared = dxHit * dxHit + dyHit * dyHit;

            return distanceSquared;

        }

    void trajMeasurementDistanceSquared
        ( const TrajectoryMeasurement& lhs,
          const TrajectoryMeasurement& rhs,
          float& distanceSquared, float& dxSquared, float& dySquared) {

            std::pair<float, float> lhsLocalXY = getLocalXY(lhs);
            std::pair<float, float> rhsLocalXY = getLocalXY(rhs);

            float dxHit = lhsLocalXY.first  - rhsLocalXY.first;
            float dyHit = lhsLocalXY.second - rhsLocalXY.second;

            dxSquared = dxHit * dxHit;
            dySquared = dyHit * dyHit;

            distanceSquared = dxSquared + dySquared;

        }

    void trajMeasurementDistance
        ( const TrajectoryMeasurement& lhs,
          const TrajectoryMeasurement& rhs,
          float& distance, float& dx, float& dy) {

            trajMeasurementDistanceSquared(lhs, rhs, distance, dx, dy);

            distance = sqrt(distance);
            dx       = sqrt(dx);
            dy       = sqrt(dy);

            if((dx == NOVAL_F) || (dy == NOVAL_F)) distance = NOVAL_F;
        }

    void getClosestOtherTrajMeasurementDistanceByLooping(const TrajectoryMeasurement& measurement,
            const edm::Handle<TrajTrackAssociationCollection>& trajTrackCollectionHandle,
            float& distance, float& dx, float& dy) {
        dx = NOVAL_F;
        dy = NOVAL_F;

        DetId detId = measurement.recHit() -> geographicalId();

        std::vector<TrajectoryMeasurement>::const_iterator closestMeasurementIt =
            trajTrackCollectionHandle -> begin() -> key -> measurements().begin();

        if(&*closestMeasurementIt == &measurement) ++closestMeasurementIt;

        double closestTrajMeasurementDistanceSquared =
            trajMeasurementDistanceSquared(measurement, *closestMeasurementIt);

        for(const auto& otherTrackKeypair: *trajTrackCollectionHandle) 
        {
            const edm::Ref<std::vector<Trajectory>> otherTraj = otherTrackKeypair.key;
            for(auto otherTrajMeasurementIt = otherTraj -> measurements().begin();
                    otherTrajMeasurementIt != otherTraj -> measurements().end(); ++otherTrajMeasurementIt) {

                if(otherTrajMeasurementIt -> recHit() -> geographicalId() != detId) continue;

                if(&*otherTrajMeasurementIt == &measurement) continue;

                float distanceSquared = 
                    trajMeasurementDistanceSquared(measurement, *otherTrajMeasurementIt);

                if(distanceSquared < closestTrajMeasurementDistanceSquared) {
                    closestMeasurementIt = otherTrajMeasurementIt;
                    closestTrajMeasurementDistanceSquared = distanceSquared;
                }

            }

        }

        trajMeasurementDistance(measurement, *closestMeasurementIt, distance, dx, dy);

    }
} // ResidHelpers

//define this as a plug-in
DEFINE_FWK_MODULE(Resid);
*/
