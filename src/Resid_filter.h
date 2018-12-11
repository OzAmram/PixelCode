
#ifndef PHASEIPIXELNTUPLIZER_H
#define PHASEIPIXELNTUPLIZER_H


// CMSSW code
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
//#include "DataFormats/TrackReco/interface/Track.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/Common/interface/ConditionsInEdm.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelRawData/interface/SiPixelRawDataError.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "TrackingTools/MeasurementDet/interface/LayerMeasurements.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitQuality.h"


// Datastructures - Keep all this in one file
// This has to be a versioned file
// It cannot go into separate files included from everywhere
//#include "../interface/DataStructures_v6.h" // 2017 July 25, CMSSW_9_2_7
#include "DPGAnalysis-SiPixelTools/PhaseIPixelNtuplizer/interface/DataStructures_v6.h"


// SiPixelCoordinates: new class for plotting Phase 0/1 Geometry
#include "DQM/SiPixelPhase1Common/interface/SiPixelCoordinates.h"

// Helpers to save canvases
//#include "DPGAnalysis-SiPixelTools/PhaseIPixelNtuplizer/interface/common_functions.h"

// ROOT Libraries
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
// #include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TEfficiency.h>

// C++
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>

// Compiler directives
#define EDM_ML_LOGDEBUG
#define ML_DEBUG


//Triplets 

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/Event.h"
#include <FWCore/Framework/interface/EventSetup.h>

#include <DataFormats/BeamSpot/interface/BeamSpot.h>


#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <DataFormats/TrackReco/interface/HitPattern.h>

#include <MagneticField/Engine/interface/MagneticField.h>

// To convert detId to subdet/layer number:
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h" //GeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"



#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"


#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryBuilder.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include <TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h>
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "DataFormats/TrackerRecHit2D/interface/TkCloner.h"

#include "TrackingTools/DetLayers/interface/DetLayer.h"


#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelLorentzAngle.h"
#include "CondFormats/DataRecord/interface/SiPixelLorentzAngleRcd.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"


class Resid_filter : public edm::EDAnalyzer
{
  using LumisectionCount = int;
    static constexpr int                  ZEROBIAS_TRIGGER_BIT           = 0;
    static constexpr int                  ZEROBIAS_BITMASK               = 1 << ZEROBIAS_TRIGGER_BIT;
    static constexpr int                  VERTEX_NUMTRACK_CUT_VAL        = 10;
    static constexpr int                  TRACK_QUALITY_HIGH_PURITY_BIT  = 2;
    static constexpr int                  TRACK_QUALITY_HIGH_PURITY_MASK = 1 << TRACK_QUALITY_HIGH_PURITY_BIT;
    static constexpr float                TRACK_PT_CUT_VAL               = 1.0f;
    static constexpr int                  TRACK_NSTRIP_CUT_VAL           = 10;
    static constexpr std::array<float, 4> TRACK_D0_CUT_BARREL_VAL        = {{0.01f, 0.02f, 0.02f, 0.02f}};
    static constexpr float                TRACK_D0_CUT_FORWARD_VAL       = 0.05f;
    static constexpr float                TRACK_DZ_CUT_BARREL_VAL        = 0.01f;
    static constexpr float                TRACK_DZ_CUT_FORWARD_VAL       = 0.5f;
    static constexpr float                MEAS_HITSEP_CUT_VAL            = 0.01f; //  100 um
    static constexpr float                HIT_CLUST_NEAR_CUT_VAL         = 0.10f; // 1000 um
    static constexpr float                BARREL_MODULE_EDGE_X_CUT       = 0.6f;
    static constexpr float                BARREL_MODULE_EDGE_Y_CUT       = 3.0f;

public:
  Resid_filter(edm::ParameterSet const& iConfig);
  virtual ~Resid_filter();
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

private:
  edm::ParameterSet iConfig_;
  std::string ntupleOutputFilename_;

  // States
  int isEventFromMc_;
  int nEvent_ = 0;
  LumisectionCount nLumisection_ = 0;

  // Options
  int isCosmicTracking_;
  int clusterSaveDownscaling_;
  int eventSaveDownscaling_;
  int saveDigiTree_;
  int npixFromDigiCollection_;
  int saveTrackTree_;
  int saveNonPropagatedExtraTrajTree_;
  int minVertexSize_;
  LumisectionCount efficiencyCalculationFrequency_;

  // Misc. data
  TFile*                                 ntupleOutputFile_;
  edm::Handle<edm::ConditionsInRunBlock> conditionsInRunBlock_;
  std::vector<std::string>               triggerNames_;
  edm::InputTag                          triggerTag_;
  std::map<uint32_t, int>                federrors_;
  // Trees
  TTree* eventTree_;
  TTree* lumiTree_;
  TTree* runTree_;
  TTree* digiTree_;
  TTree* clustTree_;
  TTree* trackTree_;
  TTree* trajTree_;
  TTree* nonPropagatedExtraTrajTree_;
  TTree* trajROCEfficiencyTree_;

  // Tree field definitions are in the interface directory
  EventData         evt_;
  LumiData          lumi_;
  RunData           run_;
  Digi              digi_;
  Cluster           clu_;
  TrackData         track_;
  TrajMeasurement   traj_;
  TrajROCEfficiency trajROCEff_;

  // Tokens
  edm::EDGetTokenT<edm::DetSetVector<SiPixelRawDataError>> rawDataErrorToken_;
  edm::EDGetTokenT<reco::VertexCollection>                 primaryVerticesToken_;
  edm::EDGetTokenT<edm::TriggerResults>                    triggerResultsToken_;
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>>   clustersToken_;
  edm::EDGetTokenT<reco::TrackCollection>                  trackCollectionToken_;
  edm::EDGetTokenT<reco::TrackCollection>                  trackCollectionGeneralToken_;
  edm::EDGetTokenT<MeasurementTrackerEvent>                measurementTrackerEventToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>>         pileupSummaryToken_;
  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>>           pixelDigiCollectionToken_;
#ifdef ADD_CHECK_PLOTS_TO_NTUPLE
  std::vector<edm::EDGetTokenT<std::vector<PSimHit>>>      simhitCollectionTokens_;
#endif
  edm::EDGetTokenT<edm::ConditionsInRunBlock>              conditionsInRunBlockToken_;

  int dropLayer = 1;

  // Tools
  SiPixelCoordinates coord_;
  const PixelClusterParameterEstimator* pixelClusterParameterEstimator1D_;
  const PixelClusterParameterEstimator* pixelClusterParameterEstimator2D_;
  const TrackerTopology*                trackerTopology_;
  const TrackerGeometry*                trackerGeometry_;
  const Propagator*                     trackerPropagator_;
  const MeasurementTracker*             measurementTracker_;
  const MeasurementTrackerEvent*        measurementTrackerEvent_;
  const MeasurementEstimator*           chi2MeasurementEstimator_;
  const TrackerTopology*                tTopo;
  //definitions from triplet code
    edm::InputTag _triggerSrc;
    std::string _ttrhBuilder;
    edm::EDGetTokenT<reco::BeamSpot>  t_offlineBeamSpot_;
    edm::EDGetTokenT<reco::VertexCollection> t_offlinePrimaryVertices_ ;
    edm::EDGetTokenT<edm::TriggerResults> t_triggerSrc_ ;
    edm::EDGetTokenT<reco::TrackCollection>  t_generalTracks_;
    // ----------member data:

    TTree *tree;
    SiPixelRecHitQuality recHitQuality;
    int hit_type;
    Double_t trackIntersectx, trackIntersecty, trackIntersectErrorx, trackIntersectErrory;
    Double_t resid1Dx, resid1Dy, resid2Dx, resid2Dy;
    Double_t err1Dx, err1Dy, err2Dx, err2Dy;
    Double_t trkEta, trkPt;
    float cotalpha, cotbeta;

    Int_t clustxmax, clustxmin, clustymin, clustymax, clustSizeY, clustSizeX;
    Int_t clust2xmax, clust2xmin, clust2ymin, clust2ymax, clust2SizeY, clust2SizeX;
    float probQ1D, probQ2D, probXY1D, probXY2D;
    int qBin1D, qBin2D;
    int numLayers;
    bool onTrack, assocTrack;
    bool hasBadPixels1D, isOnEdge1D, hasFilledProb1D;
    bool hasBadPixels2D, isOnEdge2D, hasFilledProb2D;
    bool found2ndClust;
    const static int clust_mat_size_x = 13;
    const static int clust_mat_size_y = 21;
    float clustMatrix[clust_mat_size_x][clust_mat_size_y];
    float clust_start_x, clust_start_y;

    float clust2Matrix[clust_mat_size_x][clust_mat_size_y];
    float clust2_start_x, clust2_start_y;

  // Private methods
  void setTriggerTable();

  void getEvtData(const edm::Event&, const edm::Handle<reco::VertexCollection>&,
		  const edm::Handle<edm::TriggerResults>&,
		  const edm::Handle<std::vector<PileupSummaryInfo>>&,
      const edm::Handle<edm::DetSetVector<PixelDigi>>&,
		  const edm::Handle<edmNew::DetSetVector<SiPixelCluster>>&,
		  const edm::Handle<TrajTrackAssociationCollection>&);

  int getTriggerInfo(const edm::Event&, const edm::Handle<edm::TriggerResults>&);

  float getPileupInfo(const edm::Handle<std::vector<PileupSummaryInfo>>&);

  void getDigiData(const edm::Handle<edm::DetSetVector<PixelDigi>>&);


  void getClustData(const edm::Handle<edmNew::DetSetVector<SiPixelCluster>>&);

  std::map<reco::TrackRef, TrackData> getTrackData(const edm::Handle<reco::VertexCollection>&,
						   const edm::Handle<TrajTrackAssociationCollection>&);

  void getTrajTrackData(const edm::Handle<reco::VertexCollection>&,
			const edm::Handle<edmNew::DetSetVector<SiPixelCluster>>&,
			const edm::Handle<TrajTrackAssociationCollection>&);

  bool isStripClustInTrack(const SiStripCluster *clust, reco::Track track,int clusttype);
  bool isClustInTrack(const SiPixelCluster *hit, const reco::Track track);
  const reco::Track* associateInputTrack(const reco::Track iTrack, const edm::Handle<reco::TrackCollection>& tracksGeneral);

  void getTrajTrackDataCosmics(const edm::Handle<reco::VertexCollection>&,
			       const edm::Handle<edmNew::DetSetVector<SiPixelCluster>>&,
			       const edm::Handle<TrajTrackAssociationCollection>&);

void checkAndSaveTrajMeasurementData
( const TrajectoryMeasurement& measurement,
  const edm::Handle<edmNew::DetSetVector<SiPixelCluster>>& clusterCollectionHandle,
  const reco::Track iTrack, const edm::Handle<reco::TrackCollection>& iGeneralTracks
  );

  std::vector<TrajectoryMeasurement> getLayerExtrapolatedHitsFromMeas(const TrajectoryMeasurement& trajMeasurement, int ilay);
  std::vector<TrajectoryMeasurement> getLayerIntercept (std::vector<TrajectoryMeasurement> &trajectoryMeasurements, int ilay);

  enum TrajectoryMeasurementEfficiencyQualification
  {
    EXCLUDED,
    VALIDHIT,
    MISSING
  };

  TrajectoryMeasurementEfficiencyQualification getTrajMeasurementEfficiencyQualification(const TrajectoryMeasurement& t_measurement);

  void getDisk1PropagationData(const edm::Handle<TrajTrackAssociationCollection>&);

  std::vector<TEfficiency> getDetectorPartEfficienciesInTrajTreeEntryRange(const TrajMeasurement& t_trajField, const Long64_t& t_minEntry, const Long64_t& t_maxEntry);
  void generateROCEfficiencyTree();


  void handleDefaultError(const std::string&, const std::string&, std::string);

  void handleDefaultError(const std::string&, const std::string&, std::vector<std::string>);

  void printEvtInfo(const std::string&);

  void getModuleData(ModuleData&, bool, const DetId&);

  void getRocData(ModuleData&, bool, const DetId&, const PixelDigi*);

  void getRocData(ModuleData&, bool, const DetId&, const SiPixelCluster*);

  void getRocData(ModuleData&, bool, const SiPixelRecHit*);

  void propagateTrackToLayer1(const edm::Ref<std::vector<Trajectory>>&, const reco::TrackRef);

  std::tuple<std::vector<TrajectoryMeasurement>::const_iterator, float>
  findMatchingTrajMeasurement(const GlobalPoint&, const ModuleData&,
			      const std::vector<TrajectoryMeasurement>&);

  void makeClustMatrix(const SiPixelCluster *clust, float clustMat[clust_mat_size_x][clust_mat_size_y]);

  std::pair<const SiPixelCluster*, const SiPixelCluster*> getClosestClustersOnDetSetToPoint(const edmNew::DetSet<SiPixelCluster>&,
							 const LocalPoint&, 
							 const LocalTrajectoryParameters&,
							 bool);

  float trajMeasGlobalPointDistanceSquared(const TrajectoryMeasurement&, const GlobalPoint&);

  float clusterPointDistanceSquared(const DetId&, const SiPixelCluster&, const LocalPoint&, const LocalTrajectoryParameters&, bool);

  LocalPoint clusterPointDistanceVector(const DetId&, const SiPixelCluster&, const LocalPoint&);

  float clusterPointDistance(const DetId&, const SiPixelCluster&, const LocalPoint&);

  void printTrackCompositionInfo(const edm::Ref<std::vector<Trajectory>>&,
				 const reco::TrackRef&,
				 const edm::Handle<edmNew::DetSetVector<SiPixelCluster>>,
				 const edm::Handle<reco::VertexCollection>&);
};

namespace Resid_filterHelpers 
{
  std::map<uint32_t, int>
  getFedErrors(const edm::Event&,
	       const edm::EDGetTokenT<edm::DetSetVector<SiPixelRawDataError>>&);


  bool detidIsOnPixel(const DetId&);
  bool detidIsOnStrips(const DetId&);

  bool areIdenticalModules(const ModuleData&, const ModuleData&);

  int trajectoryHasPixelHit(const edm::Ref<std::vector<Trajectory>>&);

  reco::VertexCollection::const_iterator
  findClosestVertexToTrack(const reco::TrackRef&,
			   const edm::Handle<reco::VertexCollection>&, const unsigned int&);

  TrajectoryStateOnSurface getTrajectoryStateOnSurface(const TrajectoryMeasurement&);

  std::pair<float, float> getLocalXY(const TrajectoryMeasurement&);

  float trajMeasurementDistanceSquared(const TrajectoryMeasurement&, const TrajectoryMeasurement&);

  void trajMeasurementDistanceSquared(const TrajectoryMeasurement&, const TrajectoryMeasurement&,
				      float&, float&, float&);

  void trajMeasurementDistance(const TrajectoryMeasurement&, const TrajectoryMeasurement&,
			       float&, float&, float&);


  // int getTrackParentVtxNumTracks(const edm::Handle<reco::VertexCollection>&, const reco::TrackRef);

} // NtuplizerHelpers

#endif
