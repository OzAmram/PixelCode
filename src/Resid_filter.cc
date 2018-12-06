// CMS and user include files:
#include "data_test/pixel/src/Resid_filter.h"

Resid_filter::Resid_filter(const edm::ParameterSet& iConfig// , edm::ConsumesCollector && ic
        )
{
    std::cout << "Resid_filter constructed\n";
    _triggerSrc = iConfig.getParameter<edm::InputTag>("triggerSource");
    _ttrhBuilder = iConfig.getParameter<std::string>("ttrhBuilder");
    std::cout<<_triggerSrc<<" "<<_triggerSrc.label()<<" "<<_triggerSrc.process()<<" "
        <<_triggerSrc.instance()<<" "<<std::endl;



    //Definition of parameters
    t_triggerSrc_ = consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggerSource"));
    t_offlineBeamSpot_ =    consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
    t_offlinePrimaryVertices_ =   consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));

    primaryVerticesToken_  = consumes<reco::VertexCollection> (edm::InputTag("offlinePrimaryVertices"));
    triggerResultsToken_ = consumes<edm::TriggerResults>(triggerTag_);
    pileupSummaryToken_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"));
    clustersToken_ = consumes<edmNew::DetSetVector<SiPixelCluster>>(edm::InputTag("siPixelClusters"));
    trackCollectionToken_ = consumes<reco::TrackCollection> (iConfig.getParameter<edm::InputTag>("trackInput"));
    trackCollectionGeneralToken_ = consumes<reco::TrackCollection> (iConfig.getParameter<edm::InputTag>("trackInputGeneral"));
    measurementTrackerEventToken_ = consumes<MeasurementTrackerEvent> (edm::InputTag("MeasurementTrackerEvent"));
    pixelDigiCollectionToken_ = consumes<edm::DetSetVector<PixelDigi>> (edm::InputTag("simSiPixelDigis"));
    dropLayer=iConfig.getParameter<int>("dropLayer");


    // tok_caloHH_ = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits", "HcalHits"));
    edm::Service<TFileService> fsT;
    tree = fsT->make<TTree>("tree", "tree");
    tree->Branch("trackIntersectx", &trackIntersectx);
    tree->Branch("trackIntersecty", &trackIntersecty);
    tree->Branch("trackIntersectErrorx", &trackIntersectErrorx);
    tree->Branch("trackIntersectErrory", &trackIntersectErrory);
    tree->Branch("cotalpha", &cotalpha);
    tree->Branch("cotbeta", &cotbeta);
    tree->Branch("resid1Dx", &resid1Dx);
    tree->Branch("resid1Dy", &resid1Dy);
    tree->Branch("resid2Dx", &resid2Dx);
    tree->Branch("resid2Dy", &resid2Dy);
    tree->Branch("err1Dx", &err1Dx);
    tree->Branch("err1Dy", &err1Dy);
    tree->Branch("err2Dx", &err2Dx);
    tree->Branch("err2Dy", &err2Dy);
    tree->Branch("trkEta", &trkEta);
    tree->Branch("trkPt", &trkPt);

    tree->Branch("clustSizeX", &clustSizeX);
    tree->Branch("clustSizeY", &clustSizeY);
    tree->Branch("clustxmax", &clustxmax);
    tree->Branch("clustxmin", &clustxmin);
    tree->Branch("clustymax", &clustymax);
    tree->Branch("clustymin", &clustymin);
    tree->Branch("clust_start_x", &clust_start_x);
    tree->Branch("clust_start_y", &clust_start_y);
    tree->Branch("clustMatrix", clustMatrix, "clustMatrix[13][21]/F");

    tree->Branch("clust2SizeX", &clust2SizeX);
    tree->Branch("clust2SizeY", &clust2SizeY);
    tree->Branch("clust2xmax", &clust2xmax);
    tree->Branch("clust2xmin", &clust2xmin);
    tree->Branch("clust2ymax", &clust2ymax);
    tree->Branch("clust2ymin", &clust2ymin);
    tree->Branch("clust2_start_x", &clust2_start_x);
    tree->Branch("clust2_start_y", &clust2_start_y);
    tree->Branch("clust2Matrix", clust2Matrix, "clust2Matrix[13][21]/F");

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

}
//
// destructor:
//
Resid_filter::~Resid_filter()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions:
// method called once each job just before starting event loop
//
void Resid_filter::beginJob()
{
}

//----------------------------------------------------------------------
// method called for each event:

void Resid_filter::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
    const int run = iRun.run();


}

bool Resid_filter::isClustInTrack(const SiPixelCluster *clust, reco::Track track){

    if(clust == nullptr){
        printf("Clust was null! \n");
        return false;
    }
    int minPixelRow = clust->minPixelRow();
    int maxPixelRow = clust->maxPixelRow();
    int minPixelCol = clust->minPixelCol();
    int maxPixelCol = clust->maxPixelCol();
    // rechits from track
    for( trackingRecHit_iterator irecHit = track.recHitsBegin();
            irecHit != track.recHitsEnd(); ++irecHit ) {

        DetId detId = (*irecHit)->geographicalId();
        if(!Resid_filterHelpers::detidIsOnPixel(detId)){
            //printf("rechit not on pixel \n");
            continue;
        }

        const SiPixelCluster* clustTrk = nullptr;
        const SiPixelRecHit *hitTrk = static_cast<const SiPixelRecHit*>((*irecHit) -> hit());
        clustTrk = hitTrk -> cluster().get();

        if(clust != nullptr && clustTrk!= nullptr){
            if(  clustTrk->maxPixelRow()  == maxPixelRow &&
                    clustTrk->minPixelRow()  == minPixelRow &&
                    clustTrk->maxPixelCol()  == maxPixelCol &&
                    clustTrk->maxPixelCol()  == minPixelCol ) {
                return true;
            }
        }
    }

    return false;
}

bool Resid_filter::isStripClustInTrack(const SiStripCluster *clust, reco::Track track){

    if(clust == nullptr){
        printf("Strip Clust was null! \n");
        return false;
    }
    unsigned int clustsize = clust->amplitudes().size();
    int first  = clust->firstStrip();     
    // rechits from track
    for( trackingRecHit_iterator irecHit = track.recHitsBegin();
            irecHit != track.recHitsEnd(); ++irecHit ) {

        DetId detId = (*irecHit)->geographicalId();
        if(!Resid_filterHelpers::detidIsOnPixel(detId)){
            //printf("rechit not on pixel \n");
            continue;
        }

        const SiStripCluster* clustTrk = nullptr;
        if(const SiStripRecHit2D * hitTrk = dynamic_cast<const SiStripRecHit2D *>((*irecHit) ->hit()))
            clustTrk = &(*hitTrk -> cluster());

        if(clust != nullptr && clustTrk!= nullptr){
            if(clustsize == clustTrk->amplitudes().size() &&
               first == clust->firstStrip())
                return true;
            }
        }
    

    return false;
}

const reco::Track* Resid_filter::associateInputTrack(const reco::Track iTrack,const edm::Handle<reco::TrackCollection>& tracksGeneral) {
    //  const reco::Track* pTrack = nullptr;

    // general tracks
    for( reco::TrackCollection::const_iterator iTrackGeneral = tracksGeneral->begin();
            iTrackGeneral != tracksGeneral->end(); ++iTrackGeneral ) {


        //we only use refit tracks with > 10 pt
        if(iTrackGeneral->pt() < 9.) continue;
        bool found_trk = true;
        printf("lets check new track \n");

        // check for which track all recthits from refitted track are on general track                                                  
        // rechits from refitted track
        for( trackingRecHit_iterator irecHit = iTrack.recHitsBegin();
                irecHit != iTrack.recHitsEnd(); ++irecHit ) {

            if(!(*irecHit)->isValid()) continue;
            DetId detId = (*irecHit)->geographicalId();
            if(Resid_filterHelpers::detidIsOnPixel(detId)) {
                const SiPixelRecHit *hit = static_cast<const SiPixelRecHit*>((*irecHit)->hit());
                const SiPixelCluster* clust = nullptr;
                clust = hit -> cluster().get();
                if(clust !=nullptr && !isClustInTrack(clust, *iTrackGeneral)){
                    printf("clust not in track, breaking \n");
                    found_trk = false;
                    break;

                }
                else{
                    printf("clust in track \n");
                }
            }
            else if(Resid_filterHelpers::detidIsOnStrips(detId)) {
                const SiStripCluster* clust = nullptr;

                if(const SiStripRecHit2D * rechit = 
                        dynamic_cast<const SiStripRecHit2D *>((*irecHit)->hit()))
                {   
                    printf("Strip type 1 \n");
                    clust = &(*rechit -> cluster());
                }
                //check if it is a matched SiStripMatchedRecHit2D
                else  if(const SiStripRecHit1D * rechit = 
                        dynamic_cast<const SiStripRecHit1D *>((*irecHit)->hit()))
                {   
                    printf("Strip type 2 \n");
                    clust = &(*rechit -> cluster());
                }
                //check if it is a matched SiStripMatchedRecHit2D
                else  if(const SiStripMatchedRecHit2D * rechit = 
                        dynamic_cast<const SiStripMatchedRecHit2D *>((*irecHit)->hit()))
                {   
                    printf("Strip type 3 \n");
                    continue;
                }
                //check if it is a  ProjectedSiStripRecHit2D
                else if(const ProjectedSiStripRecHit2D * rechit = 
                        dynamic_cast<const ProjectedSiStripRecHit2D *>((*irecHit)->hit())) {
                    printf("Strip type 4 \n");
                    continue;
                }

                if(clust !=nullptr && !isStripClustInTrack(clust, *iTrackGeneral)){
                    printf("strip clust not in track, breaking \n");
                    found_trk = false;
                    break;

                }
                else{
                    printf("strip clust in track \n");
                }
            }
            



        }// rechits from refitted track  

        if(found_trk) return &(*iTrackGeneral);

    } // general tracks
    printf("couldn't associate track!! \n");

    return nullptr;
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
void Resid_filter::checkAndSaveTrajMeasurementData
( const TrajectoryMeasurement& measurement,
  const edm::Handle<edmNew::DetSetVector<SiPixelCluster>>& clusterCollectionHandle,
  const reco::Track iTrack,
  const edm::Handle<reco::TrackCollection>& iGeneralTracks
  ) {

    // Check if the measurement infos can be read                                                                                              
    if(!measurement.updatedState().isValid()) return;

    // given a measurement get tracking rechit                                                                                                 
    // for that detId that should be on pixel                                                                                                  
    TransientTrackingRecHit::ConstRecHitPointer recHit = measurement.recHit();
    DetId detId = recHit -> geographicalId();
    if(!Resid_filterHelpers::detidIsOnPixel(detId)) return;

    // get estimated trajectory state at a surface close to where you want hits                                                                
    TrajectoryStateOnSurface trajStateOnSurface = Resid_filterHelpers::getTrajectoryStateOnSurface(measurement);

    // skip hits with undeterminable positions                                                                                                 
    if(!(trajStateOnSurface.isValid())) return;

    // from estimated trajectory state get global and local position (x,y,z,validhit,roc)                                                      
    //GlobalPoint globalPosition     = trajStateOnSurface.globalPosition();
    LocalPoint  localTrkPosition      = trajStateOnSurface.localPosition();
    LocalError  localTrkPositionError = trajStateOnSurface.localError().positionError();

    // get track local parameters: momentum,direction,angles                                                                                   
    LocalTrajectoryParameters trajectoryParameters = trajStateOnSurface.localParameters();
    auto trajectoryMomentum = trajectoryParameters.momentum();
    LocalVector localTrackDirection = trajectoryMomentum / trajectoryMomentum.mag();

    //float alpha = atan2(localDir.z(), localDir.x());
    //float beta = atan2(localDir.z(), localDir.y());
    cotalpha = trajectoryParameters.dxdz();
    cotbeta = trajectoryParameters.dydz();


    // from tracking rechit: check if is valid(what does that mean?) and a hit                                                                 
    // if yes: => save params                                                                                                                  
    // else: => loop over siPixelClusters on detector and get closest cluster to track local position                                          
    const SiPixelCluster* clust = nullptr;
    const SiPixelCluster* clust2 = nullptr;
    const edmNew::DetSetVector<SiPixelCluster>::const_iterator clustersOnDet =
        clusterCollectionHandle -> find(detId);

    if(recHit -> isValid() && recHit -> hit() != 0) {
        const SiPixelRecHit *hit = static_cast<const SiPixelRecHit*>(recHit -> hit());
        hit_type = 1;
        clust = hit -> cluster().get();
            if(clustersOnDet != clusterCollectionHandle -> end()){
                const SiPixelCluster *clust1a;
                std::tie(clust1a, clust2) = getClosestClustersOnDetSetToPoint(*clustersOnDet, localTrkPosition, trajectoryParameters, false);
                if(clust1a != nullptr){
                    float dist1= clusterPointDistanceSquared(detId, *clust, localTrkPosition, trajectoryParameters, false);
                    float dist2= clusterPointDistanceSquared(detId, *clust1a, localTrkPosition, trajectoryParameters, false);
                    printf("Clust 1 is clust 1a: %d \n", abs(dist1 -dist2) < 1E-6);
                }
            }
        //float currentMinValueSquared = clusterPointDistanceSquared(detId, *clust_1D, localPosition, trajectoryParameters, false);
        //printf("recHit on trajectory measure (?), distance to traj %f and layer %i \n ",currentMinValueSquared,tTopo->pxbLayer(detId));
    } else {
        if(clusterCollectionHandle.isValid()) {
            if(clustersOnDet != clusterCollectionHandle -> end()){
                hit_type = 2;
                std::tie(clust, clust2) = getClosestClustersOnDetSetToPoint(*clustersOnDet, localTrkPosition, trajectoryParameters, false);
            }
        }
    }
    if(clust == nullptr){
        printf("Null ptrs \n");
        return;
    }
    clust_start_x = (float) clust->minPixelRow();
    clust_start_y = (float) clust->minPixelCol();

    const PixelGeomDetUnit* pixelDet = dynamic_cast<const PixelGeomDetUnit*> (trackerGeometry_->idToDet(detId));
    auto topol = &(pixelDet->specificTopology());
    auto proxyT = dynamic_cast<const ProxyPixelTopology*>(topol);
    auto loc_trk_pred = Topology::LocalTrackPred(localTrkPosition.x(), localTrkPosition.y(), cotalpha, cotbeta);
    auto clust_start_lp = proxyT->localPosition(MeasurementPoint(clust_start_x, clust_start_y), loc_trk_pred );

    if(clust2 != nullptr){
        found2ndClust = true;
        clust2_start_x = (float) clust2->minPixelRow();
        clust2_start_y = (float) clust2->minPixelCol();
        auto clust2_start_lp = proxyT->localPosition(MeasurementPoint(clust2_start_x, clust2_start_y), loc_trk_pred );
        clust2_start_x = 1E4* clust2_start_lp.x();
        clust2_start_y = 1E4* clust2_start_lp.y();

        clust2SizeX = clust2->sizeX();
        clust2SizeY = clust2->sizeY();
        clust2xmin = clust2->minPixelRow();
        clust2xmax = clust2->maxPixelRow();
        clust2ymin = clust2->minPixelCol();
        clust2ymax = clust2->maxPixelCol();
        makeClustMatrix(clust2, clust2Matrix);
    }
    else{
        found2ndClust = false;
    }


    auto genAssocTrack = associateInputTrack(iTrack,iGeneralTracks);
    if(genAssocTrack != nullptr){
        onTrack = isClustInTrack(clust, *genAssocTrack);
    }
    else{
        onTrack = false;
    }
    printf("isOnTrack = %i \n", (int) onTrack);



    //printf("Row, col = %i, %i; localx, localy =  %.4f %.4f \n", clust->minPixelRow(), clust->minPixelCol(), clust_start_lp.x(), clust_start_lp.y());
    //printf("Track intersect is at %.4f +/- %.4f  %.4f +/- %.4f \n", localPosition.x(), sqrt(localPositionError.xx()), localPosition.y(), sqrt(localPositionError.yy()));

    trackIntersectx = 1E4 *localTrkPosition.x();                                                                                                   
    trackIntersecty = 1E4 *localTrkPosition.y();                                                                                              

    trackIntersectErrorx = 1E4 *sqrt(localTrkPositionError.xx());
    trackIntersectErrory = 1E4 *sqrt(localTrkPositionError.yy());

    clust_start_x = 1E4* clust_start_lp.x();
    clust_start_y = 1E4* clust_start_lp.y();

    clustSizeX = clust->sizeX();
    clustSizeY = clust->sizeY();
    clustxmin = clust->minPixelRow();
    clustxmax = clust->maxPixelRow();
    clustymin = clust->minPixelCol();
    clustymax = clust->maxPixelCol();
    makeClustMatrix(clust, clustMatrix);

    onTrack = false;


    LocalPoint clust1DCoords, clust2DCoords;
    LocalError clust1DError, clust2DError;
    SiPixelRecHitQuality::QualWordType QW_1D, QW_2D;
    const GeomDetUnit* geomDetUnit = trackerGeometry_ -> idToDetUnit(detId);
    std::tie(clust1DCoords, clust1DError, QW_1D) = pixelClusterParameterEstimator1D_ -> getParameters(*clust, *geomDetUnit, trajectoryParameters);
    std::tie(clust2DCoords, clust2DError, QW_2D) = pixelClusterParameterEstimator2D_ -> getParameters(*clust, *geomDetUnit, trajectoryParameters);
    // if cluster exists: get cluster parameters: charge,size, local and global position, and pixels charge and position                       

    resid1Dx = (1E4* clust1DCoords.x() - trackIntersectx);
    resid1Dy = (1E4* clust1DCoords.y() - trackIntersecty);

    err1Dx = 1E4 * sqrt(clust1DError.xx());
    err1Dy = 1E4 * sqrt(clust1DError.yy());


    resid2Dx = (1E4* clust2DCoords.x() - trackIntersectx);
    resid2Dy = (1E4* clust2DCoords.y() - trackIntersecty);

    err2Dx = 1E4 * sqrt(clust2DError.xx());
    err2Dy = 1E4 * sqrt(clust2DError.yy());

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
    printf("Valid Extrap hit, filling tree \n");
    

}
//
//make a 2d array of cluster charges
void Resid_filter::makeClustMatrix(const SiPixelCluster *clust, float clustMat[clust_mat_size_x][clust_mat_size_y]){

   int row_offset = clust->minPixelRow();
   int col_offset = clust->minPixelCol();

   //--- Copy clust's pixels (calibrated in electrons) into clustMatrix;
   memset( clustMat, 0, sizeof(float)*clust_mat_size_x*clust_mat_size_y );   // Wipe it clean.
   for (int i=0 ; i!=clust->size(); ++i )
   {
       auto pix = clust->pixel(i);
       int irow = int(pix.x) - row_offset;
       int icol = int(pix.y) - col_offset;
       // &&& Do we ever get pixels that are out of bounds ???  Need to check.
       if ( (irow<clust_mat_size_x) & (icol<clust_mat_size_y) ) clustMat[irow][icol] =  float(pix.adc);


   }
   return;
}

// get closest cluster to point given clusters on detector                                                                                   
    std::pair<const SiPixelCluster*, const SiPixelCluster*> Resid_filter::getClosestClustersOnDetSetToPoint
(const edmNew::DetSet<SiPixelCluster>& clustersOnDet, const LocalPoint& referencePoint, const LocalTrajectoryParameters & ltp, bool use2D)
{
    if(clustersOnDet.empty()) {
        printf("clusters not on det \n");
	return std::pair(nullptr, nullptr);
    }

    const DetId detId = clustersOnDet.id();
    const SiPixelCluster* minDistanceCluster = clustersOnDet.begin();
    float currentMinValueSquared = clusterPointDistanceSquared(detId, *minDistanceCluster, referencePoint, ltp, use2D);

    const SiPixelCluster* minDistance2Cluster = nullptr;
    float currentMin2ValueSquared = 99999.;
    

    for(const auto& cluster: clustersOnDet) {
        float currentDistanceSquared = clusterPointDistanceSquared(detId, cluster, referencePoint, ltp, use2D);
        if(currentDistanceSquared < currentMinValueSquared) {
            minDistance2Cluster = minDistanceCluster;
            currentMin2ValueSquared = currentMinValueSquared;
            currentMinValueSquared = std::move(currentDistanceSquared);
            minDistanceCluster = &cluster;
        }
        else if(currentDistanceSquared < currentMin2ValueSquared){
            currentMin2ValueSquared = std::move(currentDistanceSquared);
            minDistance2Cluster = &cluster;
        }
    }

    
    printf("closest cluster on det: distance %f \n",currentMinValueSquared);
    printf("2nd closest cluster on det: distance %f \n",currentMin2ValueSquared);
    return std::pair(minDistanceCluster, minDistance2Cluster);
}

// given cluster and reference point calculate position from CPE and calculate distance to ref. point                                        
float Resid_filter::clusterPointDistanceSquared(const DetId& detId, const SiPixelCluster& cluster, const LocalPoint& referencePoint,\
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



std::vector<TrajectoryMeasurement> Resid_filter::getLayerExtrapolatedHitsFromMeas(const TrajectoryMeasurement& trajMeasurement, int layer)
{

    // Last layer 2 or disk 1 mesurement is to be propagated to layer 1 if possible
    // Only propagating valid measurements
    std::unique_ptr<LayerMeasurements> layerMeasurements
        (new LayerMeasurements(*measurementTracker_, *measurementTrackerEvent_));

    const DetLayer* pixelBarrelLayer = 
        measurementTracker_ -> geometricSearchTracker() -> pixelBarrelLayers().at(layer-1);
    //printf("trying to get hit on layer %i \n", layer);
    //printf("seqNum is %i \n", pixelBarrelLayer->seqNum());
    //printf("size is %i \n", (int) measurementTracker_ -> geometricSearchTracker() -> pixelBarrelLayers().size());

    return layerMeasurements -> measurements(*pixelBarrelLayer, trajMeasurement.updatedState(),
            *trackerPropagator_, *chi2MeasurementEstimator_);

}
void Resid_filter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup // , edm::ConsumesCollector && ic
        ){
    using namespace std;
    using namespace edm;
    using namespace reco;
    using namespace math;
    printf("analyze start \n");

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



    //--------------------------------------------------------------------
    //Retrieve tracker topology from geometry
    edm::ESHandle<TrackerTopology> tTopoH;
    iSetup.get<TrackerTopologyRcd>().get(tTopoH);
    tTopo=tTopoH.product();




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
    iEvent.getByToken( trackCollectionToken_, tracks );

    if( tracks.failedToGet() || !tracks.isValid() ){
        printf("Couldn't get tracks!\n");
        return;
    }
    //
    // general tracks: (not dropping layer)
    Handle<reco::TrackCollection> tGeneralTracks;
    iEvent.getByToken( trackCollectionGeneralToken_, tGeneralTracks );
    if( tGeneralTracks.failedToGet() || !tGeneralTracks.isValid() ){
      printf("Couldn't get general tracks!\n");
      return;
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



        //--------------------------------------------------------------------------
        // rec hits from track extra:

        if( iTrack->extra().isNull() ) continue;//next track
        if( ! iTrack->extra().isAvailable() ) continue;//next track

        TransientTrack tTrack = theB->build(*iTrack);

        //uint32_t outerDetId = 0;
        //float rmax = 0.;
        //
        uint32_t innerDetId = 0;
        uint32_t layer2_detid = 0;
        float rmin = 0.;


        // std::cout << " imod 0 "<< imod<< std::endl;

        edm::OwnVector<TrackingRecHit> recHitVector, recHitVectorNo1; // for seed
        bool found_layer2 = false;


        TrackingRecHit *layer2_hit;


        Trajectory::RecHitContainer coTTRHvec, coTTRHvecNo1; // for fit, constant
        TrajectoryStateOnSurface initialTSOS = tTrack.innermostMeasurementState();

        // loop over recHits on this track:

        //TrajectoryStateOnSurface initialTSOS = tTrack.outermostMeasurementState();

        edm::OwnVector<TrackingRecHit> layeri_hits;
        bool skip_track = false;

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
            if (ilay == dropLayer){
                printf("lay %i hit. Should have been dropped?? \n", ilay);
                skip_track = true;
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
        if(skip_track) continue;



        PTrajectoryStateOnDet PTraj = trajectoryStateTransform::persistentState( initialTSOS, innerDetId );
        const TrajectorySeed seed( PTraj, recHitVector, alongMomentum );

        //if( idbg ) cout << "  have seed\n";
        


        Trajectory refitTrajectory = theFitter->fitOne( seed, coTTRHvec, initialTSOS );
        if(refitTrajectory.isValid()){

                auto trajectoryMeasurements = refitTrajectory.measurements();
                // Trajectory.measurements:
                auto layerTrajectoryMeasurements = getLayerIntercept(trajectoryMeasurements, dropLayer);
                if (layerTrajectoryMeasurements.size() == 0){
                    //dummy vec
                    continue;
                }
                if(layerTrajectoryMeasurements.size() > 2) printf("More than 2 measures \n");


                //CRIS's FUNCTION GOES HERE!
                // should pass tree? find out where to end loop
                checkAndSaveTrajMeasurementData(layerTrajectoryMeasurements.front(), clusterCollectionHandle, *iTrack, tGeneralTracks) ;

            }//refitted trajectory
            else{
                printf("refit failed \n");
            }






    }// loop over tracks
    printf("loop over \n");

}
void Resid_filter::getModuleData(ModuleData &mod, bool online, const DetId &detId)
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
std::vector<TrajectoryMeasurement> Resid_filter::getLayerIntercept (std::vector<TrajectoryMeasurement> &trajectoryMeasurements, int ilay) 
{


    // std::cout << "Searching for the first layer trajectory measurements... ";
    // Not for the faint hearted ...
    int layToFind;
    //if(ilay ==1) layToFind =2;
    //else layToFind = ilay-1;
    layToFind = ilay+1;
    auto firstLayerTrajMeasurementIt = 
        std::find_if(trajectoryMeasurements.begin(), 
                trajectoryMeasurements.end(),
                [&] (const TrajectoryMeasurement& measurement) {
                ModuleData mod;
                if (!(measurement.recHit() ->isValid())) return false;
                getModuleData(mod, 1, measurement.recHit() -> geographicalId());
                return mod.det == 0 && mod.layer == layToFind;
                });
    // Check if the last non-layer traj measurement is valid 
    //int idx = firstLayer2TrajMeasurementIt - trajectoryMeasurements.begin();
    //auto firstLayer2TrajMeasurementRecHit = trajectoryMeasurements.at(idx).recHit();
    
    std::vector<TrajectoryMeasurement> dummy_vec;
    if(firstLayerTrajMeasurementIt == trajectoryMeasurements.end()){
        printf("cant find hit on layer %i \n", layToFind);
        return dummy_vec;
    }
    auto firstLayerTrajMeasurementRecHit = firstLayerTrajMeasurementIt -> recHit();
    if(firstLayerTrajMeasurementRecHit == nullptr){
        std::cout << "Invalid rechit pointer." << std::endl;
        return dummy_vec;
    }
    if(!(firstLayerTrajMeasurementRecHit -> isValid())){
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




    std::vector<TrajectoryMeasurement> extrapolatedHitsOnLayer = getLayerExtrapolatedHitsFromMeas(*firstLayerTrajMeasurementIt,ilay);

    printf("Size of extrap hits is %i \n", (int) extrapolatedHitsOnLayer.size());

    return extrapolatedHitsOnLayer;
}

void Resid_filter::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
void Resid_filter::beginLuminosityBlock(edm::LuminosityBlock const& iLumi,
        edm::EventSetup const& iSetup) {}
void Resid_filter::endLuminosityBlock(edm::LuminosityBlock const& iLumi,
        edm::EventSetup const& iSetup) {}
//----------------------------------------------------------------------
// method called just after ending the event loop:
//
void Resid_filter::endJob() {


}


namespace Resid_filterHelpers
{




    bool detidIsOnPixel(const DetId& detid) {
        if (detid.det()!=DetId::Tracker) return false;
        if (detid.subdetId() == PixelSubdetector::PixelBarrel) return true;
        if (detid.subdetId() == PixelSubdetector::PixelEndcap) return true;
        return false;
    }
    bool detidIsOnStrips(const DetId &detid){
        if (detid.det()!=DetId::Tracker) return false;
        unsigned int subdetId = detid.subdetId();
        if (subdetId==SiStripDetId::TIB||subdetId==SiStripDetId::TOB
                || subdetId==SiStripDetId::TID||subdetId==SiStripDetId::TEC) return true;
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

} // Resid_filterHelpers

//define this as a plug-in
DEFINE_FWK_MODULE(Resid_filter);
