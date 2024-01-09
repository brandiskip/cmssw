// Package:    SiOuterTrackerV
// Class:      SiOuterTrackerV

// Original Author:  Emily MacDonald

// system include files
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

// user include files
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "OuterTrackerMonitorTrackingParticles.h"

OuterTrackerMonitorTrackingParticles::OuterTrackerMonitorTrackingParticles(const edm::ParameterSet &iConfig)
    : m_topoToken(esConsumes()), conf_(iConfig) {
  topFolderName_ = conf_.getParameter<std::string>("TopFolderName");
  trackingParticleToken_ =
      consumes<std::vector<TrackingParticle>>(conf_.getParameter<edm::InputTag>("trackingParticleToken"));
  ttStubMCTruthToken_ =
      consumes<TTStubAssociationMap<Ref_Phase2TrackerDigi_>>(conf_.getParameter<edm::InputTag>("MCTruthStubInputTag"));
  ttClusterMCTruthToken_ = consumes<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>>(
      conf_.getParameter<edm::InputTag>("MCTruthClusterInputTag"));
  ttTrackMCTruthToken_ = consumes<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>>(
      conf_.getParameter<edm::InputTag>("MCTruthTrackInputTag"));
  magneticFieldToken_ = esConsumes<MagneticField, IdealMagneticFieldRecord>(magneticFieldInputTag_);
  getTokenTrackerGeom_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();
  edm::InputTag L1StubInputTag = conf_.getParameter<edm::InputTag>("L1StubInputTag");
  ttStubToken_ = consumes<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > >(L1StubInputTag);
  L1Tk_minNStub = conf_.getParameter<int>("L1Tk_minNStub");         // min number of stubs in the track
  L1Tk_maxChi2dof = conf_.getParameter<double>("L1Tk_maxChi2dof");  // maximum chi2/dof of the track
  TP_minNStub = conf_.getParameter<int>("TP_minNStub");             // min number of stubs in the tracking particle to
  TP_minNLayersStub = conf_.getParameter<int>("TP_minNLayersStub"); //min number of layers with stubs in the tracking particle to consider matching
  TP_minPt = conf_.getParameter<double>("TP_minPt");      // min pT to consider matching
  TP_maxEta = conf_.getParameter<double>("TP_maxEta");    // max eta to consider matching
  TP_maxVtxZ = conf_.getParameter<double>("TP_maxVtxZ");  // max vertZ (or z0) to consider matching
}

OuterTrackerMonitorTrackingParticles::~OuterTrackerMonitorTrackingParticles() = default;

// member functions

float OuterTrackerMonitorTrackingParticles::phiOverBendCorrection(bool isBarrel, float stub_z, float stub_r, const TrackerTopology* tTopo, uint32_t detid, const GeomDetUnit* det0, const GeomDetUnit* det1) {
    // Get R0, R1, Z0, Z1 values
    // split resolution values between positive and negative z
    float R0 = det0->position().perp();
    float R1 = det1->position().perp();
    float Z0 = det0->position().z();
    float Z1 = det1->position().z();

    bool isTiltedBarrel = (isBarrel && tTopo->tobSide(detid) != 3);
    
    float tiltAngle = 0; // Initialize to 0 (meaning no tilt, in the endcaps)
    if (isTiltedBarrel) {
        float deltaR = std::abs(R1 - R0);
        float deltaZ = (R1 - R0 > 0) ? (Z1 - Z0) : -(Z1 - Z0); // if module parallel, tilt angle should be π/2 and deltaZ would approach zero
        hist_deltaZ->Fill(deltaZ);
        hist_deltaR->Fill(deltaR);
        tiltAngle = atan(deltaR / std::abs(deltaZ));
        hist_tiltAngle->Fill(tiltAngle);
        hist_tiltAngle_vs_Z0->Fill(Z0, tiltAngle);
        hist_deltaR_vs_deltaZ->Fill(Z0, R0);
        hist_Z0_vs_deltaZ->Fill(deltaZ, Z0);
        hist_R0_vs_deltaR->Fill(deltaR, R0);
      }

    float correction;
    if (isBarrel && tTopo->tobSide(detid) != 3) {  // Assuming this condition represents tiltedBarrel
        correction = cos(tiltAngle) * std::abs(stub_z) / stub_r + sin(tiltAngle);
        hist_cosTiltAngle->Fill(cos(tiltAngle));
        hist_sinTiltAngle->Fill(sin(tiltAngle));
    } else if (isBarrel) {
        correction = 1;
    } else {
        correction = std::abs(stub_z) / stub_r; // if tiltAngle = 0, stub (not module) is parallel to the beam line, if tiltAngle = 90, stub is perpendicular to beamline
    }

    return correction;
  }

  //std::vector<double> OuterTrackerMonitorTrackingParticles::getTPDerivedCoords(edm::Ptr<TrackingParticle> my_tp, bool isBarrel, double modMaxZ, double modMinZ, float modMaxR, float modMinR) const {
  std::vector<double> OuterTrackerMonitorTrackingParticles::getTPDerivedCoords(edm::Ptr<TrackingParticle> my_tp, bool isBarrel, double modMaxZ, double modMinZ, float stub_r) const {


  double tp_phi = -99;
  double tp_r = -99;
  double tp_z = -99;
  double bfield_ = 3.8;  //B-field in T
  double c_ = 2.99792458E10;  //speed of light cm/s

  // Get values from the tracking particle my_tp
  double myTP_pt = my_tp->pt();
  double myTP_charge = my_tp->charge();
  float myTP_z0 = my_tp->vertex().z();
  double myTP_t = my_tp->tanl();
  double myTP_rinv = (myTP_charge * bfield_) / (myTP_pt);

  if (isBarrel) { 
      //tp_r = (modMaxR + modMinR) / 2;
      tp_r = stub_r;
      tp_phi = my_tp->p4().phi() - std::asin(tp_r * myTP_rinv * c_ / 2.0E13);
      tp_phi = reco::reduceRange(tp_phi);  
      tp_z = myTP_z0 + (2.0E13 / c_) * myTP_t * (1 / myTP_rinv) * std::asin(tp_r * myTP_rinv * c_ / 2.0E13);
  } else {
      tp_z = (modMaxZ + modMinZ) / 2;
      tp_phi = my_tp->p4().phi() - (tp_z - myTP_z0) * myTP_rinv * c_ / 2.0E13 / myTP_t; 
      tp_phi = reco::reduceRange(tp_phi);
      tp_r = 2 / myTP_rinv * std::sin((tp_z - myTP_z0) * myTP_rinv * c_ / 2.0E13 / myTP_t);
  }

  hist_tp_phi->Fill(tp_phi);

  std::vector<double> tpDerived_coords{tp_r, tp_z, tp_phi};
  return tpDerived_coords;
}

  
// ------------ method called for each event  ------------
void OuterTrackerMonitorTrackingParticles::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // Tracking Particles
  edm::Handle<std::vector<TrackingParticle>> trackingParticleHandle;
  iEvent.getByToken(trackingParticleToken_, trackingParticleHandle);

  // Truth Association Maps
  edm::Handle<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>> MCTruthTTTrackHandle;
  iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);
  edm::Handle<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>> MCTruthTTClusterHandle;
  iEvent.getByToken(ttClusterMCTruthToken_, MCTruthTTClusterHandle);
  edm::Handle<TTStubAssociationMap<Ref_Phase2TrackerDigi_>> MCTruthTTStubHandle;
  iEvent.getByToken(ttStubMCTruthToken_, MCTruthTTStubHandle);
  edm::Handle<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > > TTStubHandle;
  iEvent.getByToken(ttStubToken_, TTStubHandle);
  edm::ESHandle<TrackerGeometry> tGeomHandle = iSetup.getHandle(getTokenTrackerGeom_);

  // Geometries
  const TrackerGeometry *const theTrackerGeom = tGeomHandle.product();
  const TrackerTopology *const tTopo = &iSetup.getData(m_topoToken);

  edm::InputTag L1StubInputTag("TTStubsFromPhase2TrackerDigis","StubAccepted");

  // Loop over tracking particles
  int this_tp = 0;
  for (const auto &iterTP : *trackingParticleHandle) {
    edm::Ptr<TrackingParticle> tp_ptr(trackingParticleHandle, this_tp);
    this_tp++;

    // int tmp_eventid = iterTP.eventId().event();
    float tmp_tp_pt = iterTP.pt();
    float tmp_tp_phi = iterTP.phi();
    float tmp_tp_eta = iterTP.eta();

    //Calculate nLayers variable
    std::vector<edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>, TTStub<Ref_Phase2TrackerDigi_>>>
        theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);

    int hasStubInLayer[11] = {0};
    for (unsigned int is = 0; is < theStubRefs.size(); is++) {
      DetId detid(theStubRefs.at(is)->getDetId());
      int layer = -1;
      if (detid.subdetId() == StripSubdetector::TOB)
        layer = static_cast<int>(tTopo->layer(detid)) - 1;  //fill in array as entries 0-5
      else if (detid.subdetId() == StripSubdetector::TID)
        layer = static_cast<int>(tTopo->layer(detid)) + 5;  //fill in array as entries 6-10

      //treat genuine stubs separately (==2 is genuine, ==1 is not)
      if (MCTruthTTStubHandle->findTrackingParticlePtr(theStubRefs.at(is)).isNull() && hasStubInLayer[layer] < 2)
        hasStubInLayer[layer] = 1;
      else
        hasStubInLayer[layer] = 2;
    }

    int nStubLayerTP = 0;
    for (int isum = 0; isum < 11; isum++) {
      if (hasStubInLayer[isum] >= 1)
        nStubLayerTP += 1;
    }

    if (std::fabs(tmp_tp_eta) > TP_maxEta)
      continue;
    // Fill the 1D distribution plots for tracking particles, to monitor change in stub definition
    if (tmp_tp_pt > TP_minPt && nStubLayerTP >= TP_minNLayersStub) {
      trackParts_Pt->Fill(tmp_tp_pt);
      trackParts_Eta->Fill(tmp_tp_eta);
      trackParts_Phi->Fill(tmp_tp_phi);
    }

    // if (TP_select_eventid == 0 && tmp_eventid != 0)
    //   continue;  //only care about tracking particles from the primary interaction for efficiency/resolution
    int nStubTP = -1;
    if (MCTruthTTStubHandle.isValid()) {
      std::vector<edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>, TTStub<Ref_Phase2TrackerDigi_>>>
          theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);
      nStubTP = (int)theStubRefs.size();
      }
    if (MCTruthTTClusterHandle.isValid() && MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).empty())
      continue;

    float tmp_tp_vz = iterTP.vz();
    float tmp_tp_vx = iterTP.vx();
    float tmp_tp_vy = iterTP.vy();
    float tmp_tp_charge = tp_ptr->charge();
    int tmp_tp_pdgid = iterTP.pdgId();

    // ----------------------------------------------------------------------------------------------
    // calculate d0 and VtxZ propagated back to the IP, pass if greater than max
    // VtxZ
    float tmp_tp_t = tan(2.0 * atan(1.0) - 2.0 * atan(exp(-tmp_tp_eta)));
    float delx = -tmp_tp_vx;
    float dely = -tmp_tp_vy;
    float K = 0.01 * 0.5696 / tmp_tp_pt * tmp_tp_charge;  // curvature correction
    float A = 1. / (2. * K);
    float tmp_tp_x0p = delx - A * sin(tmp_tp_phi);
    float tmp_tp_y0p = dely + A * cos(tmp_tp_phi);
    float tmp_tp_rp = sqrt(tmp_tp_x0p * tmp_tp_x0p + tmp_tp_y0p * tmp_tp_y0p);
    static double pi = 4.0 * atan(1.0);
    float delphi = tmp_tp_phi - atan2(-K * tmp_tp_x0p, K * tmp_tp_y0p);
    if (delphi < -pi)
      delphi += 2.0 * pi;
    if (delphi > pi)
      delphi -= 2.0 * pi;

    float tmp_tp_VtxZ = tmp_tp_vz + tmp_tp_t * delphi / (2.0 * K);
    float tmp_tp_VtxR = sqrt(tmp_tp_vx * tmp_tp_vx + tmp_tp_vy * tmp_tp_vy);
    float tmp_tp_d0 = tmp_tp_charge * tmp_tp_rp - (1. / (2. * K));

    // simpler formula for d0, in cases where the charge is zero:
    // https://github.com/cms-sw/cmssw/blob/master/DataFormats/TrackReco/interface/TrackBase.h
    float other_d0 = -tmp_tp_vx * sin(tmp_tp_phi) + tmp_tp_vy * cos(tmp_tp_phi);
    tmp_tp_d0 = tmp_tp_d0 * (-1);  // fix d0 sign
    if (K == 0) {
      tmp_tp_d0 = other_d0;
      tmp_tp_VtxZ = tmp_tp_vz;
    }
    if (std::fabs(tmp_tp_VtxZ) > TP_maxVtxZ)
      continue;

    // To make efficiency plots where the denominator has NO stub cuts
    if (tmp_tp_VtxR < 1.0) {
      tp_pt->Fill(tmp_tp_pt);  //pT effic, no cut on pT, but VtxR cut
      if (tmp_tp_pt <= 10)
        tp_pt_zoom->Fill(tmp_tp_pt);  //pT effic, no cut on pT, but VtxR cut
    }
    if (tmp_tp_pt < TP_minPt)
      continue;
    tp_VtxR->Fill(tmp_tp_VtxR);  // VtxR efficiency has no cut on VtxR
    if (tmp_tp_VtxR > 1.0)
      continue;
    tp_eta->Fill(tmp_tp_eta);
    tp_d0->Fill(tmp_tp_d0);
    tp_VtxZ->Fill(tmp_tp_VtxZ);

    if (nStubTP < TP_minNStub || nStubLayerTP < TP_minNLayersStub)
      continue;  //nStub cut not included in denominator of efficiency plots

    // Find clusters related to tracking particle
    std::vector<edm::Ref<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_> >, TTCluster<Ref_Phase2TrackerDigi_> >> associatedClusters = MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr);

    std::map<DetId, int> clusterToStubCountMap;

    // Loop through each cluster and check if it's genuine
    for (std::size_t k = 0; k < associatedClusters.size(); ++k) {

      DetId clusdetid = associatedClusters[k]->getDetId();
      clusterToStubCountMap[clusdetid] = 0;
      DetId stackDetid = tTopo->stack(clusdetid);
      const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit(clusdetid);
      const PixelGeomDetUnit* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(det0);
      const PixelTopology* topol = dynamic_cast<const PixelTopology*>(&(theGeomDet->specificTopology())); 
      GlobalPoint coords = theGeomDet->surface().toGlobal(topol->localPosition(associatedClusters[k]->findAverageLocalCoordinatesCentered()));

      bool isGenuine = MCTruthTTClusterHandle->isGenuine(associatedClusters[k]);
      if (!isGenuine) 
        continue;

      gen_clusters->Fill(tmp_tp_pt);
      if (tmp_tp_pt > 0 && tmp_tp_pt <= 10)
          gen_clusters_zoom->Fill(tmp_tp_pt);
      
      // if there are stubs on the same detid, loop on those stubs
      if (TTStubHandle->find(stackDetid) != TTStubHandle->end()) {
        edmNew::DetSet< TTStub<Ref_Phase2TrackerDigi_> > stubs = (*TTStubHandle)[stackDetid];
        for (auto stubIter = stubs.begin(); stubIter != stubs.end(); ++stubIter) {
          auto stubRef = edmNew::makeRefTo(TTStubHandle, stubIter);
          if (!MCTruthTTStubHandle->isGenuine(stubRef)) {
              continue; // Skip to the next iteration if the stub is not genuine
          }

          GlobalPoint coords0 = theGeomDet->surface().toGlobal(topol->localPosition(stubIter->clusterRef(0)->findAverageLocalCoordinatesCentered()));
          GlobalPoint coords1 = theGeomDet->surface().toGlobal(topol->localPosition(stubIter->clusterRef(1)->findAverageLocalCoordinatesCentered()));

          if (coords.x() == coords0.x() || coords.x() == coords1.x()) {
            edm::Ptr<TrackingParticle> stubTP = MCTruthTTStubHandle->findTrackingParticlePtr(edmNew::makeRefTo(TTStubHandle, stubIter));
            if (stubTP.isNull()) 
              continue;

            float stub_tp_pt = stubTP->pt();
            if (stub_tp_pt == tmp_tp_pt){
              std::cout << "stub_tp_pt: " << stub_tp_pt << std::endl;
              std::cout << "tmp_tp_pt: " << tmp_tp_pt << std::endl;
              clusterToStubCountMap[clusdetid]++;
              gen_clusters_if_stub->Fill(tmp_tp_pt);
              if (tmp_tp_pt > 0 && tmp_tp_pt <= 10)
                gen_clusters_if_stub_zoom->Fill(tmp_tp_pt);
              }
            }
          }
        }
      }

      // Check for clusters matched with multiple stubs
      for (const auto& entry : clusterToStubCountMap) {
        if (entry.second > 1){
          std::cout << "Cluster ID " << entry.first << " has " << entry.second << " matched stub(s)." << std::endl;
        }  
      }


    // ----------------------------------------------------------------------------------------------
    // look for L1 tracks matched to the tracking particle
    int tp_nMatch = 0;
    int i_track = -1;
    float i_chi2dof = 99999;
    if (MCTruthTTTrackHandle.isValid()) {
      std::vector<edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>>> matchedTracks =
          MCTruthTTTrackHandle->findTTTrackPtrs(tp_ptr);

      // ----------------------------------------------------------------------------------------------
      // loop over matched L1 tracks
      // here, "match" means tracks that can be associated to a TrackingParticle
      // with at least one hit of at least one of its clusters
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SLHCTrackerTriggerSWTools#MC_truth_for_TTTrack
      int trkCounter = 0;
      for (const auto &thisTrack : matchedTracks) {
        if (!MCTruthTTTrackHandle->isGenuine(thisTrack))
          continue;
        // ----------------------------------------------------------------------------------------------
        // further require L1 track to be (loosely) genuine, that there is only
        // one TP matched to the track
        // + have >= L1Tk_minNStub stubs for it to be a valid match
        int tmp_trk_nstub = thisTrack->getStubRefs().size();
        if (tmp_trk_nstub < L1Tk_minNStub)
          continue;
        float dmatch_pt = 999;
        float dmatch_eta = 999;
        float dmatch_phi = 999;
        int match_id = 999;

        edm::Ptr<TrackingParticle> my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(thisTrack);
        dmatch_pt = std::fabs(my_tp->p4().pt() - tmp_tp_pt);
        dmatch_eta = std::fabs(my_tp->p4().eta() - tmp_tp_eta);
        dmatch_phi = std::fabs(my_tp->p4().phi() - tmp_tp_phi);
        match_id = my_tp->pdgId();
        float tmp_trk_chi2dof = thisTrack->chi2Red();

        // ensure that track is uniquely matched to the TP we are looking at!
        if (dmatch_pt < 0.1 && dmatch_eta < 0.1 && dmatch_phi < 0.1 && tmp_tp_pdgid == match_id) {
          tp_nMatch++;
          if (i_track < 0 || tmp_trk_chi2dof < i_chi2dof) {
            i_track = trkCounter;
            i_chi2dof = tmp_trk_chi2dof;
          }
        }
        trkCounter++;
      }  // end loop over matched L1 tracks

      if (tp_nMatch < 1)
        continue;
      // Get information on the matched tracks
      float tmp_matchtrk_pt = -999;
      float tmp_matchtrk_eta = -999;
      float tmp_matchtrk_phi = -999;
      float tmp_matchtrk_VtxZ = -999;
      float tmp_matchtrk_chi2dof = -999;
      int tmp_matchTrk_nStub = -999;
      float tmp_matchtrk_d0 = -999;

      tmp_matchtrk_pt = matchedTracks[i_track]->momentum().perp();
      tmp_matchtrk_eta = matchedTracks[i_track]->momentum().eta();
      tmp_matchtrk_phi = matchedTracks[i_track]->momentum().phi();
      tmp_matchtrk_VtxZ = matchedTracks[i_track]->z0();
      tmp_matchtrk_chi2dof = matchedTracks[i_track]->chi2Red();
      tmp_matchTrk_nStub = (int)matchedTracks[i_track]->getStubRefs().size();

      //for d0
      float tmp_matchtrk_x0 = matchedTracks[i_track]->POCA().x();
      float tmp_matchtrk_y0 = matchedTracks[i_track]->POCA().y();
      tmp_matchtrk_d0 = -tmp_matchtrk_x0 * sin(tmp_matchtrk_phi) + tmp_matchtrk_y0 * cos(tmp_matchtrk_phi);

      //Add cuts for the matched tracks, numerator
      if (tmp_matchTrk_nStub < L1Tk_minNStub || tmp_matchtrk_chi2dof > L1Tk_maxChi2dof)
        continue;

      // fill matched track histograms (if passes all criteria)
      match_tp_pt->Fill(tmp_tp_pt);
      if (tmp_tp_pt > 0 && tmp_tp_pt <= 10)
        match_tp_pt_zoom->Fill(tmp_tp_pt);
      match_tp_eta->Fill(tmp_tp_eta);
      match_tp_d0->Fill(tmp_tp_d0);
      match_tp_VtxR->Fill(tmp_tp_VtxR);
      match_tp_VtxZ->Fill(tmp_tp_VtxZ);

      // Eta and pT histograms for resolution
      float pt_diff = tmp_matchtrk_pt - tmp_tp_pt;
      float pt_res = pt_diff / tmp_tp_pt;
      float eta_res = tmp_matchtrk_eta - tmp_tp_eta;
      float phi_res = tmp_matchtrk_phi - tmp_tp_phi;
      float VtxZ_res = tmp_matchtrk_VtxZ - tmp_tp_VtxZ;
      float d0_res = tmp_matchtrk_d0 - tmp_tp_d0;

      // fill total resolution histograms
      res_pt->Fill(pt_diff);
      res_ptRel->Fill(pt_res);
      res_eta->Fill(eta_res);

      // Fill resolution plots for different abs(eta) bins:
      // (0, 0.7), (0.7, 1.0), (1.0, 1.2), (1.2, 1.6), (1.6, 2.0), (2.0, 2.4)
      if (std::fabs(tmp_tp_eta) >= 0 && std::fabs(tmp_tp_eta) < 0.7) {
        reseta_eta0to0p7->Fill(eta_res);
        resphi_eta0to0p7->Fill(phi_res);
        resVtxZ_eta0to0p7->Fill(VtxZ_res);
        resd0_eta0to0p7->Fill(d0_res);
        if (tmp_tp_pt >= 2 && tmp_tp_pt < 3)
          respt_eta0to0p7_pt2to3->Fill(pt_res);
        else if (tmp_tp_pt >= 3 && tmp_tp_pt < 8)
          respt_eta0to0p7_pt3to8->Fill(pt_res);
        else if (tmp_tp_pt >= 8)
          respt_eta0to0p7_pt8toInf->Fill(pt_res);
      } else if (std::fabs(tmp_tp_eta) >= 0.7 && std::fabs(tmp_tp_eta) < 1.0) {
        reseta_eta0p7to1->Fill(eta_res);
        resphi_eta0p7to1->Fill(phi_res);
        resVtxZ_eta0p7to1->Fill(VtxZ_res);
        resd0_eta0p7to1->Fill(d0_res);
        if (tmp_tp_pt >= 2 && tmp_tp_pt < 3)
          respt_eta0p7to1_pt2to3->Fill(pt_res);
        else if (tmp_tp_pt >= 3 && tmp_tp_pt < 8)
          respt_eta0p7to1_pt3to8->Fill(pt_res);
        else if (tmp_tp_pt >= 8)
          respt_eta0p7to1_pt8toInf->Fill(pt_res);
      } else if (std::fabs(tmp_tp_eta) >= 1.0 && std::fabs(tmp_tp_eta) < 1.2) {
        reseta_eta1to1p2->Fill(eta_res);
        resphi_eta1to1p2->Fill(phi_res);
        resVtxZ_eta1to1p2->Fill(VtxZ_res);
        resd0_eta1to1p2->Fill(d0_res);
        if (tmp_tp_pt >= 2 && tmp_tp_pt < 3)
          respt_eta1to1p2_pt2to3->Fill(pt_res);
        else if (tmp_tp_pt >= 3 && tmp_tp_pt < 8)
          respt_eta1to1p2_pt3to8->Fill(pt_res);
        else if (tmp_tp_pt >= 8)
          respt_eta1to1p2_pt8toInf->Fill(pt_res);
      } else if (std::fabs(tmp_tp_eta) >= 1.2 && std::fabs(tmp_tp_eta) < 1.6) {
        reseta_eta1p2to1p6->Fill(eta_res);
        resphi_eta1p2to1p6->Fill(phi_res);
        resVtxZ_eta1p2to1p6->Fill(VtxZ_res);
        resd0_eta1p2to1p6->Fill(d0_res);
        if (tmp_tp_pt >= 2 && tmp_tp_pt < 3)
          respt_eta1p2to1p6_pt2to3->Fill(pt_res);
        else if (tmp_tp_pt >= 3 && tmp_tp_pt < 8)
          respt_eta1p2to1p6_pt3to8->Fill(pt_res);
        else if (tmp_tp_pt >= 8)
          respt_eta1p2to1p6_pt8toInf->Fill(pt_res);
      } else if (std::fabs(tmp_tp_eta) >= 1.6 && std::fabs(tmp_tp_eta) < 2.0) {
        reseta_eta1p6to2->Fill(eta_res);
        resphi_eta1p6to2->Fill(phi_res);
        resVtxZ_eta1p6to2->Fill(VtxZ_res);
        resd0_eta1p6to2->Fill(d0_res);
        if (tmp_tp_pt >= 2 && tmp_tp_pt < 3)
          respt_eta1p6to2_pt2to3->Fill(pt_res);
        else if (tmp_tp_pt >= 3 && tmp_tp_pt < 8)
          respt_eta1p6to2_pt3to8->Fill(pt_res);
        else if (tmp_tp_pt >= 8)
          respt_eta1p6to2_pt8toInf->Fill(pt_res);
      } else if (std::fabs(tmp_tp_eta) >= 2.0 && std::fabs(tmp_tp_eta) <= 2.4) {
        reseta_eta2to2p4->Fill(eta_res);
        resphi_eta2to2p4->Fill(phi_res);
        resVtxZ_eta2to2p4->Fill(VtxZ_res);
        resd0_eta2to2p4->Fill(d0_res);
        if (tmp_tp_pt >= 2 && tmp_tp_pt < 3)
          respt_eta2to2p4_pt2to3->Fill(pt_res);
        else if (tmp_tp_pt >= 3 && tmp_tp_pt < 8)
          respt_eta2to2p4_pt3to8->Fill(pt_res);
        else if (tmp_tp_pt >= 8)
          respt_eta2to2p4_pt8toInf->Fill(pt_res);
      }
    }  //if MC TTTrack handle is valid 
  }    //end loop over tracking particles

  // ----------------------------------------------------------------------------------------------
  // float myTP_phi = -999;
  int myTP_charge = -999;
  float myTP_pt = -999;
  float myTP_eta = -999;
  float myTP_dxy = -999;
  float stub_r = -999;
  float stub_phi = -999;
  double bfield_{3.8};  //B-field in T
  double c_{2.99792458E10};  //speed of light cm/s

  // Loop over L1 stubs
  for (auto gd = theTrackerGeom->dets().begin(); gd != theTrackerGeom->dets().end(); gd++) {
    DetId detid = (*gd)->geographicalId();
    if (detid.subdetId() != StripSubdetector::TOB && detid.subdetId() != StripSubdetector::TID)
        continue;

    if (!tTopo->isLower(detid))
        continue; // Only process lower part of the stack

    DetId stackDetid = tTopo->stack(detid); // Get Stub module DetId

    if (TTStubHandle->find(stackDetid) == TTStubHandle->end())
        continue;

    // Get the DetSets of the Clusters
    edmNew::DetSet<TTStub<Ref_Phase2TrackerDigi_> > stubs = (*TTStubHandle)[stackDetid];
    numOfStubs->Fill(stubs.size());

    const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit(detid);
    const GeomDetUnit* det1 = theTrackerGeom->idToDetUnit(tTopo->partnerDetId(detid));

    // Calculate detector module min and max positions
    float modMinR = std::min(det0->position().perp(), det1->position().perp());
    float modMaxR = std::max(det0->position().perp(), det1->position().perp());
    float modMinZ = std::min(det0->position().z(), det1->position().z());
    float modMaxZ = std::max(det0->position().z(), det1->position().z());
    float sensorSpacing = sqrt((modMaxR - modMinR) * (modMaxR - modMinR) + (modMaxZ - modMinZ) * (modMaxZ - modMinZ));

    // Calculate the strip pitch from the detector's topology
    const PixelGeomDetUnit* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(det0);
    const PixelTopology& topo = theGeomDet->specificTopology();
    float stripPitch = topo.pitch().first;

    // Loop over all individual stubs in a specific detector unit
    for (auto stubIter = stubs.begin(); stubIter != stubs.end(); ++stubIter) {
      edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> >, TTStub<Ref_Phase2TrackerDigi_> > tempStubPtr = edmNew::makeRefTo(TTStubHandle, stubIter);

      if (!MCTruthTTStubHandle->isGenuine(tempStubPtr)) {
          continue; // Skip non-genuine stubs
      }

      // Process tracking particle information
      edm::Ptr<TrackingParticle> my_tp = MCTruthTTStubHandle->findTrackingParticlePtr(tempStubPtr);
      if (my_tp.isNull())
          continue;

      int isBarrel = 0;
        int layer = -999999;
        if (detid.subdetId() == StripSubdetector::TOB) {
          isBarrel = 1;
          layer = static_cast<int>(tTopo->layer(detid));
        } else if (detid.subdetId() == StripSubdetector::TID) {
          isBarrel = 0;
          layer = static_cast<int>(tTopo->layer(detid));
        } else {
          edm::LogVerbatim("Tracklet") << "WARNING -- neither TOB or TID stub, shouldn't happen...";
          layer = -1;
        }

      int isPSmodule = 0;
      if (topo.nrows() == 960)
        isPSmodule = 1;

      // Calculate local coordinates of clusters
      MeasurementPoint innerClusterCoords = tempStubPtr->clusterRef(0)->findAverageLocalCoordinatesCentered();
      MeasurementPoint outerClusterCoords = tempStubPtr->clusterRef(1)->findAverageLocalCoordinatesCentered();

      // Convert local coordinates to global positions
      Global3DPoint innerClusterGlobalPos = theGeomDet->surface().toGlobal(theGeomDet->topology().localPosition(innerClusterCoords));
      Global3DPoint outerClusterGlobalPos = theGeomDet->surface().toGlobal(theGeomDet->topology().localPosition(outerClusterCoords));

      // Determine maximum and minimum Z positions of the stubs
      float stub_maxZ = std::max(innerClusterGlobalPos.z(), outerClusterGlobalPos.z());
      float stub_minZ = std::min(innerClusterGlobalPos.z(), outerClusterGlobalPos.z());

      // stub parameters
      stub_phi = innerClusterGlobalPos.phi();
      stub_r = innerClusterGlobalPos.perp();
      float stub_z = innerClusterGlobalPos.z();

      myTP_charge = my_tp->charge();
      myTP_pt = my_tp->p4().pt();
      myTP_eta = my_tp->p4().eta();

      float myTP_x0 = my_tp->vertex().x();
      float myTP_y0 = my_tp->vertex().y();
      myTP_dxy = sqrt(myTP_x0 * myTP_x0 + myTP_y0 * myTP_y0);
      
      if (myTP_charge == 0) continue;
      if (myTP_pt < TP_minPt) continue;
      if (std::abs(myTP_eta) > TP_maxEta) continue;

      std::vector<double> tpDerivedCoords = getTPDerivedCoords(my_tp, isBarrel, modMaxZ, modMinZ, stub_r);
      float tp_z = tpDerivedCoords[1];
      float tp_phi = tpDerivedCoords[2];

      float trigDisplace = tempStubPtr->rawBend();
      float trigOffset = tempStubPtr->bendOffset();
      float trigPos = tempStubPtr->innerClusterPosition();
      float trigBend = tempStubPtr->bendFE();

      if (!isBarrel && stub_maxZ < 0.0){
        trigBend = -trigBend; 
      }

      bool isTiltedBarrel = (isBarrel && tTopo->tobSide(detid) != 3);
      
      float correctionValue = phiOverBendCorrection(isBarrel, stub_z, stub_r, tTopo, detid, det0, det1);
      float trackBend = -(sensorSpacing * stub_r * bfield_ * c_ * myTP_charge) /
                      (stripPitch * 2.0E13 * myTP_pt * correctionValue);
      
      float bendRes = trackBend - trigBend;
      float zRes = tp_z - stub_z;
      float phiRes = tp_phi - stub_phi;

      // Fill is PS module histograms
      if (isPSmodule) {
        z_res_isPS->Fill(zRes);
        stub_phi_res_isPS->Fill(phiRes);
        if (isBarrel) {
          z_res_isPS_barrel->Fill(zRes);
          stub_phi_res_isPS_barrel->Fill(phiRes);
        } else {
          z_res_isPS_endcap->Fill(zRes);
          stub_phi_res_isPS_endcap->Fill(phiRes);
        }
      }
      // Fill is 2S module histograms
      else {
        z_res_is2S->Fill(zRes);
        stub_phi_res_is2S->Fill(phiRes);
        if (isBarrel) {
          z_res_is2S_barrel->Fill(zRes);
          stub_phi_res_is2S_barrel->Fill(phiRes);
        } else {
          z_res_is2S_endcap->Fill(zRes);
          stub_phi_res_is2S_endcap->Fill(phiRes);
        }
      }

      if (stub_maxZ > 0){
        stub_z_res_g0->Fill(zRes);
      } else {
        stub_z_res_l0->Fill(zRes);
      }

      
      if (std::abs(bendRes) > 1.5){
        TP_pT_bendres_g1p5->Fill(myTP_pt);
        TP_eta_bendres_g1p5->Fill(myTP_eta);
        TP_dxy_bendres_g1p5->Fill(myTP_dxy);
    } else {
        TP_pT_bendres_0_to_1p5->Fill(myTP_pt);
        TP_eta_bendres_0_to_1p5->Fill(myTP_eta);
        TP_dxy_bendres_0_to_1p5->Fill(myTP_dxy);
    }

      // fill histograms for associated tracking particle from genuine stub
      TP_pT->Fill(myTP_pt);
      track_bend->Fill(trackBend);

      // associated stub is genuine and has an associated tracking particle
      stub_rawBend->Fill(trigDisplace);
      stub_bendOffset->Fill(trigOffset);
      stub_inClusPos->Fill(trigPos);
      stub_bendFE->Fill(trigBend);
      hist_stub_phi->Fill(stub_phi);
      hist_tp_z->Fill(tp_z);
      hist_stub_z->Fill(stub_z);
      bend_res->Fill(bendRes);
      stub_z_res->Fill(zRes);
      stub_phi_res->Fill(phiRes);
      
      if (isBarrel){
        stub_phi_res_barrel->Fill(phiRes);
        z_res_barrel->Fill(zRes);
        trackPhi_vs_stubPhi_barrel->Fill(tp_phi, stub_phi);
        hist_tp_z_barrel->Fill(tp_z);
        hist_stub_z_barrel->Fill(stub_z);
      } else{
        stub_phi_res_endcap->Fill(phiRes);
        z_res_endcap->Fill(zRes);
        trackPhi_vs_stubPhi_endcap->Fill(tp_phi, stub_phi);
        hist_tp_z_endcap->Fill(tp_z);
        hist_stub_z_endcap->Fill(stub_z);
      }

      // Fill 2D histogram
      trackBend_vs_stubBend->Fill(trackBend, trigBend);
      trackPhi_vs_stubPhi->Fill(tp_phi, stub_phi);
      stub_maxZ_vs_minZ->Fill(stub_maxZ, stub_minZ);
      modMaxZ_vs_modMinZ->Fill(modMaxZ, modMinZ);
      stub_Z_vs_tpZ->Fill(stub_z, tp_z);
      modMaxR_vs_modMinR->Fill(modMaxR, modMinR);

      if (myTP_pt > 10){
        stub_maxZ_vs_minZ_highPt->Fill(stub_maxZ, stub_minZ);
        modMaxZ_vs_modMinZ_highPt->Fill(modMaxZ, modMinZ);
        stub_Z_vs_tpZ_highPt->Fill(stub_z, tp_z);
      }

      // Fill the histogram for barrel stubs
      if (isBarrel == 1) {
        barrelHistogram_genuine->Fill(layer); // layer is the variable determined from your provided code
        barrel_trackBend_vs_stubBend->Fill(trackBend, trigBend);
        bend_res_barrel->Fill(bendRes);
        if (layer == 1){
          bend_res_barrel_L1->Fill(bendRes);
          barrel_trackBend_vs_stubBend_L1->Fill(trackBend, trigBend);
          z_res_barrel_L1->Fill(zRes);
          stub_Z_vs_tpZ_L1->Fill(stub_z, tp_z);
          if (isTiltedBarrel){
            stub_Z_vs_tpZ_L1_tilted->Fill(stub_z, tp_z);
            if (myTP_pt >=2 && myTP_pt < 3){
              stub_Z_vs_tpZ_L1_tilted_pT_2to3->Fill(stub_z, tp_z);
            }
            else if (myTP_pt >= 3 && myTP_pt < 5){
              stub_Z_vs_tpZ_L1_tilted_pT_3to5->Fill(stub_z, tp_z);
            }
            else if (myTP_pt >= 5 && myTP_pt <= 10){
              stub_Z_vs_tpZ_L1_tilted_pT_5to10->Fill(stub_z, tp_z);
            }
          } 
      } else if (layer == 2){
          bend_res_barrel_L2->Fill(bendRes);
          barrel_trackBend_vs_stubBend_L2->Fill(trackBend, trigBend);
          z_res_barrel_L2->Fill(zRes);
          stub_Z_vs_tpZ_L2->Fill(stub_z, tp_z);
          if (isTiltedBarrel){
            stub_Z_vs_tpZ_L2_tilted->Fill(stub_z, tp_z);
            if (myTP_pt >=2 && myTP_pt < 3){
              stub_Z_vs_tpZ_L2_tilted_pT_2to3->Fill(stub_z, tp_z);
            }
            else if (myTP_pt >= 3 && myTP_pt < 5){
              stub_Z_vs_tpZ_L2_tilted_pT_3to5->Fill(stub_z, tp_z);
            }
            else if (myTP_pt >= 5 && myTP_pt <= 10){
              stub_Z_vs_tpZ_L2_tilted_pT_5to10->Fill(stub_z, tp_z);
            }
          }
      } else if (layer == 3){
          bend_res_barrel_L3->Fill(bendRes);
          barrel_trackBend_vs_stubBend_L3->Fill(trackBend, trigBend);
          z_res_barrel_L3->Fill(zRes);
          stub_Z_vs_tpZ_L3->Fill(stub_z, tp_z);
          if (isTiltedBarrel){
            stub_Z_vs_tpZ_L3_tilted->Fill(stub_z, tp_z);
            if (myTP_pt >=2 && myTP_pt < 3){
              stub_Z_vs_tpZ_L3_tilted_pT_2to3->Fill(stub_z, tp_z);
            }
            else if (myTP_pt >= 3 && myTP_pt < 5){
              stub_Z_vs_tpZ_L3_tilted_pT_3to5->Fill(stub_z, tp_z);
            }
            else if (myTP_pt >= 5 && myTP_pt <= 10){
              stub_Z_vs_tpZ_L3_tilted_pT_5to10->Fill(stub_z, tp_z);
            }
          }
      } else if (layer == 4){
          bend_res_barrel_L4->Fill(bendRes);
          barrel_trackBend_vs_stubBend_L4->Fill(trackBend, trigBend);
          z_res_barrel_L4->Fill(zRes);
          stub_Z_vs_tpZ_L4->Fill(stub_z, tp_z);
      } else if (layer == 5){
          bend_res_barrel_L5->Fill(bendRes);
          barrel_trackBend_vs_stubBend_L5->Fill(trackBend, trigBend);
          z_res_barrel_L5->Fill(zRes);
          stub_Z_vs_tpZ_L5->Fill(stub_z, tp_z);
      } else if (layer == 6){
          bend_res_barrel_L6->Fill(bendRes);
          barrel_trackBend_vs_stubBend_L6->Fill(trackBend, trigBend);
          z_res_barrel_L6->Fill(zRes);
          stub_Z_vs_tpZ_L6->Fill(stub_z, tp_z);
      }
        
    } else if (isBarrel == 0) {
        endcapHistogram_genuine->Fill(layer);
        endcap_trackBend_vs_stubBend->Fill(trackBend, trigBend);
        bend_res_endcap->Fill(bendRes);
        if (stub_maxZ > 0){
          endcap_disc_Fw_genuine->Fill(layer);
          endcap_fw_trackBend_vs_stubBend->Fill(trackBend, trigBend);
          bend_res_fw_endcap->Fill(bendRes);
          if (layer == 1){
            bend_res_fw_endcap_D1->Fill(bendRes);
        } else if (layer == 2) {
            bend_res_fw_endcap_D2->Fill(bendRes);
        } else if (layer == 3) {
            bend_res_fw_endcap_D3->Fill(bendRes);
        } else if (layer == 4) {
            bend_res_fw_endcap_D4->Fill(bendRes);
        } else if (layer == 5) {
            bend_res_fw_endcap_D5->Fill(bendRes);
        }
      } else {
          endcap_disc_Bw_genuine->Fill(layer);
          endcap_bw_trackBend_vs_stubBend->Fill(trackBend, trigBend);
          bend_res_bw_endcap->Fill(bendRes);
          if (layer == 1){
            bend_res_bw_endcap_D1->Fill(bendRes);
        } else if (layer == 2) {
            bend_res_bw_endcap_D2->Fill(bendRes);
        } else if (layer == 3) {
            bend_res_bw_endcap_D3->Fill(bendRes);
        } else if (layer == 4) {
            bend_res_bw_endcap_D4->Fill(bendRes);
        } else if (layer == 5) {
            bend_res_bw_endcap_D5->Fill(bendRes);
          } 
        }
      }
    }
  }
} // end of method

// ------------ method called once each job just before starting event loop
// ------------
void OuterTrackerMonitorTrackingParticles::bookHistograms(DQMStore::IBooker &iBooker,
                                                          edm::Run const &run,
                                                          edm::EventSetup const &es) {
  // Histogram setup and definitions
  std::string HistoName;
  iBooker.setCurrentFolder(topFolderName_ + "/trackParticles");

  // 1D: pT
  edm::ParameterSet psTrackParts_Pt = conf_.getParameter<edm::ParameterSet>("TH1TrackParts_Pt");
  HistoName = "trackParts_Pt";
  trackParts_Pt = iBooker.book1D(HistoName,
                                 HistoName,
                                 psTrackParts_Pt.getParameter<int32_t>("Nbinsx"),
                                 psTrackParts_Pt.getParameter<double>("xmin"),
                                 psTrackParts_Pt.getParameter<double>("xmax"));
  trackParts_Pt->setAxisTitle("p_{T} [GeV]", 1);
  trackParts_Pt->setAxisTitle("# tracking particles", 2);

  // 1D: eta
  edm::ParameterSet psTrackParts_Eta = conf_.getParameter<edm::ParameterSet>("TH1TrackParts_Eta");
  HistoName = "trackParts_Eta";
  trackParts_Eta = iBooker.book1D(HistoName,
                                  HistoName,
                                  psTrackParts_Eta.getParameter<int32_t>("Nbinsx"),
                                  psTrackParts_Eta.getParameter<double>("xmin"),
                                  psTrackParts_Eta.getParameter<double>("xmax"));
  trackParts_Eta->setAxisTitle("#eta", 1);
  trackParts_Eta->setAxisTitle("# tracking particles", 2);

  // 1D: phi
  edm::ParameterSet psTrackParts_Phi = conf_.getParameter<edm::ParameterSet>("TH1TrackParts_Phi");
  HistoName = "trackParts_Phi";
  trackParts_Phi = iBooker.book1D(HistoName,
                                  HistoName,
                                  psTrackParts_Phi.getParameter<int32_t>("Nbinsx"),
                                  psTrackParts_Phi.getParameter<double>("xmin"),
                                  psTrackParts_Phi.getParameter<double>("xmax"));
  trackParts_Phi->setAxisTitle("#phi", 1);
  trackParts_Phi->setAxisTitle("# tracking particles", 2);

  // 1D plots for efficiency
  iBooker.setCurrentFolder(topFolderName_ + "/Tracks/Efficiency");
  // pT
  edm::ParameterSet psEffic_pt = conf_.getParameter<edm::ParameterSet>("TH1Effic_pt");
  HistoName = "tp_pt";
  tp_pt = iBooker.book1D(HistoName,
                         HistoName,
                         psEffic_pt.getParameter<int32_t>("Nbinsx"),
                         psEffic_pt.getParameter<double>("xmin"),
                         psEffic_pt.getParameter<double>("xmax"));
  tp_pt->setAxisTitle("p_{T} [GeV]", 1);
  tp_pt->setAxisTitle("# tracking particles", 2);

  // Matched TP's pT
  HistoName = "match_tp_pt";
  match_tp_pt = iBooker.book1D(HistoName,
                               HistoName,
                               psEffic_pt.getParameter<int32_t>("Nbinsx"),
                               psEffic_pt.getParameter<double>("xmin"),
                               psEffic_pt.getParameter<double>("xmax"));
  match_tp_pt->setAxisTitle("p_{T} [GeV]", 1);
  match_tp_pt->setAxisTitle("# matched tracking particles", 2);

  HistoName = "gen_clusters";
  gen_clusters = iBooker.book1D(HistoName,
                         HistoName,
                         psEffic_pt.getParameter<int32_t>("Nbinsx"),
                         psEffic_pt.getParameter<double>("xmin"),
                         psEffic_pt.getParameter<double>("xmax"));
  gen_clusters->setAxisTitle("p_{T} [GeV]", 1);
  gen_clusters->setAxisTitle("# tracking particles", 2);

  HistoName = "gen_clusters_if_stub";
  gen_clusters_if_stub = iBooker.book1D(HistoName,
                         HistoName,
                         psEffic_pt.getParameter<int32_t>("Nbinsx"),
                         psEffic_pt.getParameter<double>("xmin"),
                         psEffic_pt.getParameter<double>("xmax"));
  gen_clusters_if_stub->setAxisTitle("p_{T} [GeV]", 1);
  gen_clusters_if_stub->setAxisTitle("# tracking particles", 2);

  // pT zoom (0-10 GeV)
  edm::ParameterSet psEffic_pt_zoom = conf_.getParameter<edm::ParameterSet>("TH1Effic_pt_zoom");
  HistoName = "tp_pt_zoom";
  tp_pt_zoom = iBooker.book1D(HistoName,
                              HistoName,
                              psEffic_pt_zoom.getParameter<int32_t>("Nbinsx"),
                              psEffic_pt_zoom.getParameter<double>("xmin"),
                              psEffic_pt_zoom.getParameter<double>("xmax"));
  tp_pt_zoom->setAxisTitle("p_{T} [GeV]", 1);
  tp_pt_zoom->setAxisTitle("# tracking particles", 2);

  // Matched pT zoom (0-10 GeV)
  HistoName = "match_tp_pt_zoom";
  match_tp_pt_zoom = iBooker.book1D(HistoName,
                                    HistoName,
                                    psEffic_pt_zoom.getParameter<int32_t>("Nbinsx"),
                                    psEffic_pt_zoom.getParameter<double>("xmin"),
                                    psEffic_pt_zoom.getParameter<double>("xmax"));
  match_tp_pt_zoom->setAxisTitle("p_{T} [GeV]", 1);
  match_tp_pt_zoom->setAxisTitle("# matched tracking particles", 2);

  HistoName = "gen_clusters_zoom";
  gen_clusters_zoom = iBooker.book1D(HistoName,
                         HistoName,
                         psEffic_pt_zoom.getParameter<int32_t>("Nbinsx"),
                         psEffic_pt_zoom.getParameter<double>("xmin"),
                         psEffic_pt_zoom.getParameter<double>("xmax"));
  gen_clusters_zoom->setAxisTitle("p_{T} [GeV]", 1);
  gen_clusters_zoom->setAxisTitle("# tracking particles", 2);

  HistoName = "gen_clusters_if_stub_zoom";
  gen_clusters_if_stub_zoom = iBooker.book1D(HistoName,
                         HistoName,
                         psEffic_pt_zoom.getParameter<int32_t>("Nbinsx"),
                         psEffic_pt_zoom.getParameter<double>("xmin"),
                         psEffic_pt_zoom.getParameter<double>("xmax"));
  gen_clusters_if_stub_zoom->setAxisTitle("p_{T} [GeV]", 1);
  gen_clusters_if_stub_zoom->setAxisTitle("# tracking particles", 2);

  // eta
  edm::ParameterSet psEffic_eta = conf_.getParameter<edm::ParameterSet>("TH1Effic_eta");
  HistoName = "tp_eta";
  tp_eta = iBooker.book1D(HistoName,
                          HistoName,
                          psEffic_eta.getParameter<int32_t>("Nbinsx"),
                          psEffic_eta.getParameter<double>("xmin"),
                          psEffic_eta.getParameter<double>("xmax"));
  tp_eta->setAxisTitle("#eta", 1);
  tp_eta->setAxisTitle("# tracking particles", 2);

  // Matched eta
  HistoName = "match_tp_eta";
  match_tp_eta = iBooker.book1D(HistoName,
                                HistoName,
                                psEffic_eta.getParameter<int32_t>("Nbinsx"),
                                psEffic_eta.getParameter<double>("xmin"),
                                psEffic_eta.getParameter<double>("xmax"));
  match_tp_eta->setAxisTitle("#eta", 1);
  match_tp_eta->setAxisTitle("# matched tracking particles", 2);

  // d0
  edm::ParameterSet psEffic_d0 = conf_.getParameter<edm::ParameterSet>("TH1Effic_d0");
  HistoName = "tp_d0";
  tp_d0 = iBooker.book1D(HistoName,
                         HistoName,
                         psEffic_d0.getParameter<int32_t>("Nbinsx"),
                         psEffic_d0.getParameter<double>("xmin"),
                         psEffic_d0.getParameter<double>("xmax"));
  tp_d0->setAxisTitle("d_{0} [cm]", 1);
  tp_d0->setAxisTitle("# tracking particles", 2);

  // Matched d0
  HistoName = "match_tp_d0";
  match_tp_d0 = iBooker.book1D(HistoName,
                               HistoName,
                               psEffic_d0.getParameter<int32_t>("Nbinsx"),
                               psEffic_d0.getParameter<double>("xmin"),
                               psEffic_d0.getParameter<double>("xmax"));
  match_tp_d0->setAxisTitle("d_{0} [cm]", 1);
  match_tp_d0->setAxisTitle("# matched tracking particles", 2);

  // VtxR (also known as vxy)
  edm::ParameterSet psEffic_VtxR = conf_.getParameter<edm::ParameterSet>("TH1Effic_VtxR");
  HistoName = "tp_VtxR";
  tp_VtxR = iBooker.book1D(HistoName,
                           HistoName,
                           psEffic_VtxR.getParameter<int32_t>("Nbinsx"),
                           psEffic_VtxR.getParameter<double>("xmin"),
                           psEffic_VtxR.getParameter<double>("xmax"));
  tp_VtxR->setAxisTitle("d_{xy} [cm]", 1);
  tp_VtxR->setAxisTitle("# tracking particles", 2);

  // Matched VtxR
  HistoName = "match_tp_VtxR";
  match_tp_VtxR = iBooker.book1D(HistoName,
                                 HistoName,
                                 psEffic_VtxR.getParameter<int32_t>("Nbinsx"),
                                 psEffic_VtxR.getParameter<double>("xmin"),
                                 psEffic_VtxR.getParameter<double>("xmax"));
  match_tp_VtxR->setAxisTitle("d_{xy} [cm]", 1);
  match_tp_VtxR->setAxisTitle("# matched tracking particles", 2);

  // VtxZ
  edm::ParameterSet psEffic_VtxZ = conf_.getParameter<edm::ParameterSet>("TH1Effic_VtxZ");
  HistoName = "tp_VtxZ";
  tp_VtxZ = iBooker.book1D(HistoName,
                           HistoName,
                           psEffic_VtxZ.getParameter<int32_t>("Nbinsx"),
                           psEffic_VtxZ.getParameter<double>("xmin"),
                           psEffic_VtxZ.getParameter<double>("xmax"));
  tp_VtxZ->setAxisTitle("z_{0} [cm]", 1);
  tp_VtxZ->setAxisTitle("# tracking particles", 2);

  // Matched d0
  HistoName = "match_tp_VtxZ";
  match_tp_VtxZ = iBooker.book1D(HistoName,
                                 HistoName,
                                 psEffic_VtxZ.getParameter<int32_t>("Nbinsx"),
                                 psEffic_VtxZ.getParameter<double>("xmin"),
                                 psEffic_VtxZ.getParameter<double>("xmax"));
  match_tp_VtxZ->setAxisTitle("z_{0} [cm]", 1);
  match_tp_VtxZ->setAxisTitle("# matched tracking particles", 2);

  // 1D plots for resolution
  iBooker.setCurrentFolder(topFolderName_ + "/Tracks/Resolution");

  
  edm::ParameterSet psDelta_Z = conf_.getParameter<edm::ParameterSet>("TH1delta_Z");
  HistoName = "delta_Z";
  hist_deltaZ = iBooker.book1D(HistoName,
                          HistoName,
                          psDelta_Z.getParameter<int32_t>("Nbinsx"),
                          psDelta_Z.getParameter<double>("xmin"),
                          psDelta_Z.getParameter<double>("xmax"));
  hist_deltaZ->setAxisTitle("delta_Z [cm]", 1);
  hist_deltaZ->setAxisTitle("count", 2);

  edm::ParameterSet psDelta_R = conf_.getParameter<edm::ParameterSet>("TH1delta_R");
  HistoName = "delta_R";
  hist_deltaR = iBooker.book1D(HistoName,
                          HistoName,
                          psDelta_R.getParameter<int32_t>("Nbinsx"),
                          psDelta_R.getParameter<double>("xmin"),
                          psDelta_R.getParameter<double>("xmax"));
  hist_deltaR->setAxisTitle("delta_R [cm]", 1);
  hist_deltaR->setAxisTitle("count", 2);

  edm::ParameterSet pstiltAngle = conf_.getParameter<edm::ParameterSet>("TH1tilt_Angle");
  HistoName = "tiltAngle";
  hist_tiltAngle = iBooker.book1D(HistoName,
                          HistoName,
                          pstiltAngle.getParameter<int32_t>("Nbinsx"),
                          pstiltAngle.getParameter<double>("xmin"),
                          pstiltAngle.getParameter<double>("xmax"));
  hist_tiltAngle->setAxisTitle("tiltAngle [radians]", 1);
  hist_tiltAngle->setAxisTitle("count", 2);

  edm::ParameterSet pstp_phi = conf_.getParameter<edm::ParameterSet>("TH1tp_phi");
  HistoName = "tp_phi";
  hist_tp_phi = iBooker.book1D(HistoName,
                          HistoName,
                          pstp_phi.getParameter<int32_t>("Nbinsx"),
                          pstp_phi.getParameter<double>("xmin"),
                          pstp_phi.getParameter<double>("xmax"));
  hist_tp_phi->setAxisTitle("tp_phi [radians]", 1);
  hist_tp_phi->setAxisTitle("count", 2);

  HistoName = "stub_phi";
  hist_stub_phi = iBooker.book1D(HistoName,
                          HistoName,
                          pstp_phi.getParameter<int32_t>("Nbinsx"),
                          pstp_phi.getParameter<double>("xmin"),
                          pstp_phi.getParameter<double>("xmax"));
  hist_stub_phi->setAxisTitle("stub_phi [radians]", 1);
  hist_stub_phi->setAxisTitle("count", 2);

  // full pT
  edm::ParameterSet psRes_pt = conf_.getParameter<edm::ParameterSet>("TH1Res_pt");
  HistoName = "res_pt";
  res_pt = iBooker.book1D(HistoName,
                          HistoName,
                          psRes_pt.getParameter<int32_t>("Nbinsx"),
                          psRes_pt.getParameter<double>("xmin"),
                          psRes_pt.getParameter<double>("xmax"));
  res_pt->setAxisTitle("p_{T} [GeV]", 1);
  res_pt->setAxisTitle("# tracking particles", 2);

  // Full eta
  edm::ParameterSet psRes_eta = conf_.getParameter<edm::ParameterSet>("TH1Res_eta");
  HistoName = "res_eta";
  res_eta = iBooker.book1D(HistoName,
                           HistoName,
                           psRes_eta.getParameter<int32_t>("Nbinsx"),
                           psRes_eta.getParameter<double>("xmin"),
                           psRes_eta.getParameter<double>("xmax"));
  res_eta->setAxisTitle("#eta", 1);
  res_eta->setAxisTitle("# tracking particles", 2);

  edm::ParameterSet psTP_eta = conf_.getParameter<edm::ParameterSet>("TH1TP_eta");
  HistoName = "TP_eta_bendres_g1p5";
  TP_eta_bendres_g1p5 = iBooker.book1D(HistoName,
                          HistoName,
                          psTP_eta.getParameter<int32_t>("Nbinsx"),
                          psTP_eta.getParameter<double>("xmin"),
                          psTP_eta.getParameter<double>("xmax"));
  TP_eta_bendres_g1p5->setAxisTitle("#eta_{trk}", 1);
  TP_eta_bendres_g1p5->setAxisTitle("# associated tracking particles", 2);

  HistoName = "TP_eta_bendres_0_to_1p5";
  TP_eta_bendres_0_to_1p5= iBooker.book1D(HistoName,
                          HistoName,
                          psTP_eta.getParameter<int32_t>("Nbinsx"),
                          psTP_eta.getParameter<double>("xmin"),
                          psTP_eta.getParameter<double>("xmax"));
  TP_eta_bendres_0_to_1p5->setAxisTitle("#eta_{trk}", 1);
  TP_eta_bendres_0_to_1p5->setAxisTitle("# associated tracking particles", 2);

  // Stub associated TP dxy
  edm::ParameterSet psTP_dxy = conf_.getParameter<edm::ParameterSet>("TH1TP_dxy");
  HistoName = "TP_dxy_bendres_g1p5";
  TP_dxy_bendres_g1p5 = iBooker.book1D(HistoName,
                          HistoName,
                          psTP_dxy.getParameter<int32_t>("Nbinsx"),
                          psTP_dxy.getParameter<double>("xmin"),
                          psTP_dxy.getParameter<double>("xmax"));
  TP_dxy_bendres_g1p5->setAxisTitle("dxy", 1);
  TP_dxy_bendres_g1p5->setAxisTitle("# associated tracking particles", 2);

  HistoName = "TP_dxy_bendres_0_to_1p5";
  TP_dxy_bendres_0_to_1p5 = iBooker.book1D(HistoName,
                          HistoName,
                          psTP_dxy.getParameter<int32_t>("Nbinsx"),
                          psTP_dxy.getParameter<double>("xmin"),
                          psTP_dxy.getParameter<double>("xmax"));
  TP_dxy_bendres_0_to_1p5->setAxisTitle("dxy", 1);
  TP_dxy_bendres_0_to_1p5->setAxisTitle("# associated tracking particles", 2);

  // Stub associated tp pT
  edm::ParameterSet psTP_pt = conf_.getParameter<edm::ParameterSet>("TH1TP_pt");
  HistoName = "all stub associated tp pT";
  TP_pT = iBooker.book1D(HistoName,
                          HistoName,
                          psTP_pt.getParameter<int32_t>("Nbinsx"),
                          psTP_pt.getParameter<double>("xmin"),
                          psTP_pt.getParameter<double>("xmax"));
  TP_pT->setAxisTitle("p_{T} [GeV]", 1);
  TP_pT->setAxisTitle("# associated tracking particles", 2);

  HistoName = "TP_pT_bendres_g1p5";
  TP_pT_bendres_g1p5 = iBooker.book1D(HistoName,
                          HistoName,
                          psTP_pt.getParameter<int32_t>("Nbinsx"),
                          psTP_pt.getParameter<double>("xmin"),
                          psTP_pt.getParameter<double>("xmax"));
  TP_pT_bendres_g1p5->setAxisTitle("p_{T} [GeV]", 1);
  TP_pT_bendres_g1p5->setAxisTitle("# associated tracking particles", 2);

  HistoName = "TP_pT_bendres_0_to_1p5";
  TP_pT_bendres_0_to_1p5 = iBooker.book1D(HistoName,
                          HistoName,
                          psTP_pt.getParameter<int32_t>("Nbinsx"),
                          psTP_pt.getParameter<double>("xmin"),
                          psTP_pt.getParameter<double>("xmax"));
  TP_pT_bendres_0_to_1p5->setAxisTitle("p_{T} [GeV]", 1);
  TP_pT_bendres_0_to_1p5->setAxisTitle("# associated tracking particles", 2);

  // Stub in z
  edm::ParameterSet posStubz = conf_.getParameter<edm::ParameterSet>("TH1Stub_z");
  HistoName = "stub_z";
  hist_stub_z = iBooker.book1D(HistoName,
                           HistoName,
                           posStubz.getParameter<int32_t>("Nbinsx"),
                           posStubz.getParameter<double>("xmin"),
                           posStubz.getParameter<double>("xmax"));
  hist_stub_z->setAxisTitle("radius [cm]", 1);
  hist_stub_z->setAxisTitle("counts ", 2);

  HistoName = "stub_z_barrel";
  hist_stub_z_barrel = iBooker.book1D(HistoName,
                           HistoName,
                           posStubz.getParameter<int32_t>("Nbinsx"),
                           posStubz.getParameter<double>("xmin"),
                           posStubz.getParameter<double>("xmax"));
  hist_stub_z_barrel->setAxisTitle("radius [cm]", 1);
  hist_stub_z_barrel->setAxisTitle("counts ", 2);

  HistoName = "stub_z_endcap";
  hist_stub_z_endcap = iBooker.book1D(HistoName,
                           HistoName,
                           posStubz.getParameter<int32_t>("Nbinsx"),
                           posStubz.getParameter<double>("xmin"),
                           posStubz.getParameter<double>("xmax"));
  hist_stub_z_endcap->setAxisTitle("radius [cm]", 1);
  hist_stub_z_endcap->setAxisTitle("counts ", 2);

  HistoName = "tp_z";
  hist_tp_z = iBooker.book1D(HistoName,
                           HistoName,
                           posStubz.getParameter<int32_t>("Nbinsx"),
                           posStubz.getParameter<double>("xmin"),
                           posStubz.getParameter<double>("xmax"));
  hist_tp_z->setAxisTitle("radius [cm]", 1);
  hist_tp_z->setAxisTitle("counts ", 2);

  HistoName = "tp_z_barrel";
  hist_tp_z_barrel = iBooker.book1D(HistoName,
                           HistoName,
                           posStubz.getParameter<int32_t>("Nbinsx"),
                           posStubz.getParameter<double>("xmin"),
                           posStubz.getParameter<double>("xmax"));
  hist_tp_z_barrel->setAxisTitle("radius [cm]", 1);
  hist_tp_z_barrel->setAxisTitle("counts ", 2);

  HistoName = "tp_z_endcap";
  hist_tp_z_endcap = iBooker.book1D(HistoName,
                           HistoName,
                           posStubz.getParameter<int32_t>("Nbinsx"),
                           posStubz.getParameter<double>("xmin"),
                           posStubz.getParameter<double>("xmax"));
  hist_tp_z_endcap->setAxisTitle("radius [cm]", 1);
  hist_tp_z_endcap->setAxisTitle("counts ", 2);

  HistoName = "TP_z0";
  TP_z0 = iBooker.book1D(HistoName,
                           HistoName,
                           posStubz.getParameter<int32_t>("Nbinsx"),
                           posStubz.getParameter<double>("xmin"),
                           posStubz.getParameter<double>("xmax"));
  TP_z0->setAxisTitle("radius [cm]", 1);
  TP_z0->setAxisTitle("counts ", 2);

  // Stub raw bend
  edm::ParameterSet stub_raw_bend = conf_.getParameter<edm::ParameterSet>("TH1Stub_rawBend");
  HistoName = "stub_rawBend";
  stub_rawBend = iBooker.book1D(HistoName,
                           HistoName,
                           stub_raw_bend.getParameter<int32_t>("Nbinsx"),
                           stub_raw_bend.getParameter<double>("xmin"),
                           stub_raw_bend.getParameter<double>("xmax"));
  stub_rawBend->setAxisTitle("rawBend", 1);
  stub_rawBend->setAxisTitle("counts ", 2);

  // Stub bend offset
  edm::ParameterSet stub_bend_offset = conf_.getParameter<edm::ParameterSet>("TH1Stub_bendOffset");
  HistoName = "stub_bendOffset";
  stub_bendOffset = iBooker.book1D(HistoName,
                           HistoName,
                           stub_bend_offset.getParameter<int32_t>("Nbinsx"),
                           stub_bend_offset.getParameter<double>("xmin"),
                           stub_bend_offset.getParameter<double>("xmax"));
  stub_bendOffset->setAxisTitle("rawBend", 1);
  stub_bendOffset->setAxisTitle("counts", 2);

  // Stub trigPos
  edm::ParameterSet stub_inClus_Pos = conf_.getParameter<edm::ParameterSet>("TH1Stub_inClusPos");
  HistoName = "stub_inClusPos";
  stub_inClusPos = iBooker.book1D(HistoName,
                           HistoName,
                           stub_inClus_Pos.getParameter<int32_t>("Nbinsx"),
                           stub_inClus_Pos.getParameter<double>("xmin"),
                           stub_inClus_Pos.getParameter<double>("xmax"));
  stub_inClusPos->setAxisTitle("position", 1);
  stub_inClusPos->setAxisTitle("counts ", 2);

  // Stub bendFE
  edm::ParameterSet stub_bend_FE = conf_.getParameter<edm::ParameterSet>("TH1Stub_bendFE");
  HistoName = "stub_bendFE";
  stub_bendFE = iBooker.book1D(HistoName,
                           HistoName,
                           stub_bend_FE.getParameter<int32_t>("Nbinsx"),
                           stub_bend_FE.getParameter<double>("xmin"),
                           stub_bend_FE.getParameter<double>("xmax"));
  stub_bendFE->setAxisTitle("position", 1);
  stub_bendFE->setAxisTitle("counts ", 2);

  // track bend
  edm::ParameterSet psTrackBend = conf_.getParameter<edm::ParameterSet>("TH1Track_Bend");
  HistoName = "track_bend";
  track_bend = iBooker.book1D(HistoName,
                           HistoName,
                           psTrackBend.getParameter<int32_t>("Nbinsx"),
                           psTrackBend.getParameter<double>("xmin"),
                           psTrackBend.getParameter<double>("xmax"));
  track_bend->setAxisTitle("position", 1);
  track_bend->setAxisTitle("counts ", 2);

  // Stubs in event
  edm::ParameterSet stubs_In_Event = conf_.getParameter<edm::ParameterSet>("TH1StubInEvent");
  HistoName = "# of stubs";
  numOfStubs = iBooker.book1D(HistoName,
                            HistoName,
                            stubs_In_Event.getParameter<int32_t>("Nbinsx"),
                            stubs_In_Event.getParameter<double>("xmin"),
                            stubs_In_Event.getParameter<double>("xmax"));
  numOfStubs->setAxisTitle("det layer", 1);
  numOfStubs->setAxisTitle("# of stubs ", 2);

  // stub vs tp z-coord resolution
  edm::ParameterSet psZ_Res = conf_.getParameter<edm::ParameterSet>("TH1Z_Res");
  HistoName = "z-coordinate resolution";
  stub_z_res = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  stub_z_res->setAxisTitle("tp_z - stub_z", 1);
  stub_z_res->setAxisTitle("events ", 2);

  HistoName = "z-coordinate resolution stub_z G 0";
  stub_z_res_g0 = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  stub_z_res_g0->setAxisTitle("tp_z - stub_z", 1);
  stub_z_res_g0->setAxisTitle("events ", 2);

  HistoName = "z-coordinate resolution stub_z L 0";
  stub_z_res_l0 = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  stub_z_res_l0->setAxisTitle("tp_z - stub_z", 1);
  stub_z_res_l0->setAxisTitle("events ", 2);

  // z-resoution for barrel region
  HistoName = "#Delta z barrel";
  z_res_barrel = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  z_res_barrel->setAxisTitle("tp_z - stub_z", 1);
  z_res_barrel->setAxisTitle("events ", 2);

  HistoName = "z-coordinate resolution barrel L1";
  z_res_barrel_L1 = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  z_res_barrel_L1->setAxisTitle("tp_z - stub_z", 1);
  z_res_barrel_L1->setAxisTitle("events ", 2);

  HistoName = "z-coordinate resolution barrel L2";
  z_res_barrel_L2 = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  z_res_barrel_L2->setAxisTitle("tp_z - stub_z", 1);
  z_res_barrel_L2->setAxisTitle("events ", 2);

  HistoName = "z-coordinate resolution barrel L3";
  z_res_barrel_L3 = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  z_res_barrel_L3->setAxisTitle("tp_z - stub_z", 1);
  z_res_barrel_L3->setAxisTitle("events ", 2);

  HistoName = "z-coordinate resolution barrel L4";
  z_res_barrel_L4 = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  z_res_barrel_L4->setAxisTitle("tp_z - stub_z", 1);
  z_res_barrel_L4->setAxisTitle("events ", 2);

  HistoName = "z-coordinate resolution barrel L5";
  z_res_barrel_L5 = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  z_res_barrel_L5->setAxisTitle("tp_z - stub_z", 1);
  z_res_barrel_L5->setAxisTitle("events ", 2);

  HistoName = "z-coordinate resolution barrel L6";
  z_res_barrel_L6 = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  z_res_barrel_L6->setAxisTitle("tp_z - stub_z", 1);
  z_res_barrel_L6->setAxisTitle("events ", 2);

  // z-resolution for PS modules
  HistoName = "z-coordinate resolution PS modules";
  z_res_isPS = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  z_res_isPS->setAxisTitle("tp_z - stub_z", 1);
  z_res_isPS->setAxisTitle("events ", 2);

  HistoName = "#Delta z barrel PS modules";
  z_res_isPS_barrel = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  z_res_isPS_barrel->setAxisTitle("tp_z - stub_z [cm]", 1);
  z_res_isPS_barrel->setAxisTitle("events ", 2);

  // z-resolution for 2S modules
  HistoName = "z-coordinate resolution 2S modules";
  z_res_is2S = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  z_res_is2S->setAxisTitle("tp_z - stub_z [cm]", 1);
  z_res_is2S->setAxisTitle("events ", 2);

   HistoName = "#Delta z barrel 2S modules";
  z_res_is2S_barrel = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res.getParameter<int32_t>("Nbinsx"),
                            psZ_Res.getParameter<double>("xmin"),
                            psZ_Res.getParameter<double>("xmax"));
  z_res_is2S_barrel->setAxisTitle("tp_z - stub_z [cm]", 1);
  z_res_is2S_barrel->setAxisTitle("events ", 2);

  // z-resolution endcaps only
  edm::ParameterSet psZ_Res_Endcap = conf_.getParameter<edm::ParameterSet>("TH1Z_Res_Endcap");
  HistoName = "z-coordinate resolution endcaps PS modules";
  z_res_isPS_endcap = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res_Endcap.getParameter<int32_t>("Nbinsx"),
                            psZ_Res_Endcap.getParameter<double>("xmin"),
                            psZ_Res_Endcap.getParameter<double>("xmax"));
  z_res_isPS_endcap->setAxisTitle("tp_z - stub_z", 1);
  z_res_isPS_endcap->setAxisTitle("events ", 2);

  HistoName = "z-coordinate resolution endcaps 2S modules";
  z_res_is2S_endcap = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res_Endcap.getParameter<int32_t>("Nbinsx"),
                            psZ_Res_Endcap.getParameter<double>("xmin"),
                            psZ_Res_Endcap.getParameter<double>("xmax"));
  z_res_is2S_endcap->setAxisTitle("tp_z - stub_z", 1);
  z_res_is2S_endcap->setAxisTitle("events ", 2);

  // z_resolution for endcaps
  HistoName = "z-coordinate resolution endcaps";
  z_res_endcap = iBooker.book1D(HistoName,
                            HistoName,
                            psZ_Res_Endcap.getParameter<int32_t>("Nbinsx"),
                            psZ_Res_Endcap.getParameter<double>("xmin"),
                            psZ_Res_Endcap.getParameter<double>("xmax"));
  z_res_endcap->setAxisTitle("tp_z - stub_z", 1);
  z_res_endcap->setAxisTitle("events ", 2);

  // stub vs tp phi resolution
  edm::ParameterSet psPhi_Res = conf_.getParameter<edm::ParameterSet>("TH1Phi_Res");
  HistoName = "#Delta #phi";
  stub_phi_res = iBooker.book1D(HistoName,
                                    HistoName,
                                    psPhi_Res.getParameter<int32_t>("Nbinsx"),
                                    psPhi_Res.getParameter<double>("xmin"),
                                    psPhi_Res.getParameter<double>("xmax"));
  stub_phi_res->setAxisTitle("#phi_{tp} - #phi_{stub}", 1);
  stub_phi_res->setAxisTitle("# counts", 2);

  HistoName = "stub phi resolution in barrel";
  stub_phi_res_barrel = iBooker.book1D(HistoName,
                                    HistoName,
                                    psPhi_Res.getParameter<int32_t>("Nbinsx"),
                                    psPhi_Res.getParameter<double>("xmin"),
                                    psPhi_Res.getParameter<double>("xmax"));
  stub_phi_res_barrel->setAxisTitle("#phi_{tp} - #phi_{stub}", 1);
  stub_phi_res_barrel->setAxisTitle("# counts", 2);

  HistoName = "stub phi resolution in endcap";
  stub_phi_res_endcap = iBooker.book1D(HistoName,
                                    HistoName,
                                    psPhi_Res.getParameter<int32_t>("Nbinsx"),
                                    psPhi_Res.getParameter<double>("xmin"),
                                    psPhi_Res.getParameter<double>("xmax"));
  stub_phi_res_endcap->setAxisTitle("#phi_{tp} - #phi_{stub}", 1);
  stub_phi_res_endcap->setAxisTitle("# counts", 2);

  HistoName = "stub phi resolution PS modules";
  stub_phi_res_isPS = iBooker.book1D(HistoName,
                                    HistoName,
                                    psPhi_Res.getParameter<int32_t>("Nbinsx"),
                                    psPhi_Res.getParameter<double>("xmin"),
                                    psPhi_Res.getParameter<double>("xmax"));
  stub_phi_res_isPS->setAxisTitle("#phi_{tp} - #phi_{stub}", 1);
  stub_phi_res_isPS->setAxisTitle("# counts", 2);

  HistoName = "stub phi resolution barrel PS module";
  stub_phi_res_isPS_barrel = iBooker.book1D(HistoName,
                                    HistoName,
                                    psPhi_Res.getParameter<int32_t>("Nbinsx"),
                                    psPhi_Res.getParameter<double>("xmin"),
                                    psPhi_Res.getParameter<double>("xmax"));
  stub_phi_res_isPS_barrel->setAxisTitle("#phi_{tp} - #phi_{stub}", 1);
  stub_phi_res_isPS_barrel->setAxisTitle("# counts", 2);

  HistoName = "stub phi resolution endcap PS module";
  stub_phi_res_isPS_endcap = iBooker.book1D(HistoName,
                                    HistoName,
                                    psPhi_Res.getParameter<int32_t>("Nbinsx"),
                                    psPhi_Res.getParameter<double>("xmin"),
                                    psPhi_Res.getParameter<double>("xmax"));
  stub_phi_res_isPS_endcap->setAxisTitle("#phi_{tp} - #phi_{stub}", 1);
  stub_phi_res_isPS_endcap->setAxisTitle("# counts", 2);

  HistoName = "stub phi resolution 2S modules";
  stub_phi_res_is2S = iBooker.book1D(HistoName,
                                    HistoName,
                                    psPhi_Res.getParameter<int32_t>("Nbinsx"),
                                    psPhi_Res.getParameter<double>("xmin"),
                                    psPhi_Res.getParameter<double>("xmax"));
  stub_phi_res_is2S->setAxisTitle("#phi_{tp} - #phi_{stub}", 1);
  stub_phi_res_is2S->setAxisTitle("# counts", 2);

  HistoName = "stub phi resolution barrel 2S module";
  stub_phi_res_is2S_barrel = iBooker.book1D(HistoName,
                                    HistoName,
                                    psPhi_Res.getParameter<int32_t>("Nbinsx"),
                                    psPhi_Res.getParameter<double>("xmin"),
                                    psPhi_Res.getParameter<double>("xmax"));
  stub_phi_res_is2S_barrel->setAxisTitle("#phi_{tp} - #phi_{stub}", 1);
  stub_phi_res_is2S_barrel->setAxisTitle("# counts", 2);

  HistoName = "stub phi resolution endcap 2S module";
  stub_phi_res_is2S_endcap = iBooker.book1D(HistoName,
                                    HistoName,
                                    psPhi_Res.getParameter<int32_t>("Nbinsx"),
                                    psPhi_Res.getParameter<double>("xmin"),
                                    psPhi_Res.getParameter<double>("xmax"));
  stub_phi_res_is2S_endcap->setAxisTitle("#phi_{tp} - #phi_{stub}", 1);
  stub_phi_res_is2S_endcap->setAxisTitle("# counts", 2);

  edm::ParameterSet pscos_tiltAngle = conf_.getParameter<edm::ParameterSet>("TH1cosTiltAngle");
  HistoName = "cos(tiltAngle)";
  hist_cosTiltAngle = iBooker.book1D(HistoName,
                            HistoName,
                            pscos_tiltAngle.getParameter<int32_t>("Nbinsx"),
                            pscos_tiltAngle.getParameter<double>("xmin"),
                            pscos_tiltAngle.getParameter<double>("xmax"));
  hist_cosTiltAngle->setAxisTitle("cos(tiltAngle)", 1);
  hist_cosTiltAngle->setAxisTitle("events ", 2);

  edm::ParameterSet pssin_tiltAngle = conf_.getParameter<edm::ParameterSet>("TH1sinTiltAngle");
  HistoName = "sin(tiltAngle)";
  hist_sinTiltAngle = iBooker.book1D(HistoName,
                            HistoName,
                            pssin_tiltAngle.getParameter<int32_t>("Nbinsx"),
                            pssin_tiltAngle.getParameter<double>("xmin"),
                            pssin_tiltAngle.getParameter<double>("xmax"));
  hist_sinTiltAngle->setAxisTitle("sin(tiltAngle)", 1);
  hist_sinTiltAngle->setAxisTitle("events ", 2);

  // bend resolution
  edm::ParameterSet psBend_Res = conf_.getParameter<edm::ParameterSet>("TH1Bend_Res");
  HistoName = "#Delta bend";
  bend_res = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res->setAxisTitle("stub bend - tp bend", 1);
  bend_res->setAxisTitle("events ", 2);

  HistoName = "bend resolution barrel";
  bend_res_barrel = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_barrel->setAxisTitle("stub bend - tp bend", 1);
  bend_res_barrel->setAxisTitle("counts", 2);

  HistoName = "bend resolution L1";
  bend_res_barrel_L1 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_barrel_L1->setAxisTitle("stub bend - tp bend", 1);
  bend_res_barrel_L1->setAxisTitle("counts", 2);

  HistoName = "bend resolution L2";
  bend_res_barrel_L2 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_barrel_L2->setAxisTitle("stub bend - tp bend", 1);
  bend_res_barrel_L2->setAxisTitle("counts", 2);

  HistoName = "bend resolution L3";
  bend_res_barrel_L3 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_barrel_L3->setAxisTitle("stub bend - tp bend", 1);
  bend_res_barrel_L3->setAxisTitle("counts", 2);

  HistoName = "bend resolution L4";
  bend_res_barrel_L4 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_barrel_L4->setAxisTitle("stub bend - tp bend", 1);
  bend_res_barrel_L4->setAxisTitle("counts", 2);

  HistoName = "bend resolution L5";
  bend_res_barrel_L5 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_barrel_L5->setAxisTitle("stub bend - tp bend", 1);
  bend_res_barrel_L5->setAxisTitle("counts", 2);

  HistoName = "bend resolution L6";
  bend_res_barrel_L6 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_barrel_L6->setAxisTitle("stub bend - tp bend", 1);
  bend_res_barrel_L6->setAxisTitle("counts", 2);

  HistoName = "bend resolution endcaps";
  bend_res_endcap = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_endcap->setAxisTitle("stub bend - tp bend", 1);
  bend_res_endcap->setAxisTitle("counts", 2);

  HistoName = "bend resolution endcaps fw";
  bend_res_fw_endcap = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_fw_endcap->setAxisTitle("stub bend - tp bend", 1);
  bend_res_fw_endcap->setAxisTitle("counts", 2);

  HistoName = "bend resolution +D1";
  bend_res_fw_endcap_D1 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_fw_endcap_D1->setAxisTitle("stub bend - tp bend", 1);
  bend_res_fw_endcap_D1->setAxisTitle("counts", 2);

  HistoName = "bend resolution +D2";
  bend_res_fw_endcap_D2 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_fw_endcap_D2->setAxisTitle("stub bend - tp bend", 1);
  bend_res_fw_endcap_D1->setAxisTitle("counts", 2);

  HistoName = "bend resolution +D3";
  bend_res_fw_endcap_D3 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_fw_endcap_D3->setAxisTitle("stub bend - tp bend", 1);
  bend_res_fw_endcap_D3->setAxisTitle("counts", 2);

  HistoName = "bend resolution +D4";
  bend_res_fw_endcap_D4 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_fw_endcap_D4->setAxisTitle("stub bend - tp bend", 1);
  bend_res_fw_endcap_D4->setAxisTitle("counts", 2);

  HistoName = "bend resolution +D5";
  bend_res_fw_endcap_D5 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_fw_endcap_D5->setAxisTitle("stub bend - tp bend", 1);
  bend_res_fw_endcap_D5->setAxisTitle("counts", 2);

  HistoName = "bend resolution endcaps bw";
  bend_res_bw_endcap = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_bw_endcap->setAxisTitle("stub bend - tp bend", 1);
  bend_res_bw_endcap->setAxisTitle("counts", 2);

  HistoName = "bend resolution -D1";
  bend_res_bw_endcap_D1 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_bw_endcap_D1->setAxisTitle("stub bend - tp bend", 1);
  bend_res_bw_endcap_D1->setAxisTitle("counts", 2);

  HistoName = "bend resolution -D2";
  bend_res_bw_endcap_D2 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_bw_endcap_D2->setAxisTitle("stub bend - tp bend", 1);
  bend_res_bw_endcap_D1->setAxisTitle("counts", 2);

  HistoName = "bend resolution -D3";
  bend_res_bw_endcap_D3 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_bw_endcap_D3->setAxisTitle("stub bend - tp bend", 1);
  bend_res_bw_endcap_D3->setAxisTitle("counts", 2);

  HistoName = "bend resolution -D4";
  bend_res_bw_endcap_D4 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_bw_endcap_D4->setAxisTitle("stub bend - tp bend", 1);
  bend_res_bw_endcap_D4->setAxisTitle("counts", 2);

  HistoName = "bend resolution -D5";
  bend_res_bw_endcap_D5 = iBooker.book1D(HistoName,
                            HistoName,
                            psBend_Res.getParameter<int32_t>("Nbinsx"),
                            psBend_Res.getParameter<double>("xmin"),
                            psBend_Res.getParameter<double>("xmax"));
  bend_res_bw_endcap_D5->setAxisTitle("stub bend - tp bend", 1);
  bend_res_bw_endcap_D5->setAxisTitle("counts", 2);

  // stub in barrel layers
  edm::ParameterSet barrel_layers = conf_.getParameter<edm::ParameterSet>("TH1Barrel_Layers");
  HistoName = "# genuine stubs in barrel layers";
  barrelHistogram_genuine = iBooker.book1D(HistoName,
                                  HistoName,
                                  barrel_layers.getParameter<int32_t>("Nbinsx"),
                                  barrel_layers.getParameter<double>("xmin"),
                                  barrel_layers.getParameter<double>("xmax"));
  barrelHistogram_genuine->setAxisTitle("barrel layer", 1);
  barrelHistogram_genuine->setAxisTitle("# of genuine stubs", 2);

  // stub in endcap disc
  edm::ParameterSet endcap_layers = conf_.getParameter<edm::ParameterSet>("TH1Endcap_Layers");
  HistoName = "# genuine stubs in endcap disc";
  endcapHistogram_genuine = iBooker.book1D(HistoName,
                                  HistoName,
                                  endcap_layers.getParameter<int32_t>("Nbinsx"),
                                  endcap_layers.getParameter<double>("xmin"),
                                  endcap_layers.getParameter<double>("xmax"));
  endcapHistogram_genuine->setAxisTitle("endcap disc", 1);
  endcapHistogram_genuine->setAxisTitle("# of genuine stubs", 2);

  // stub in front end endcap disc
  HistoName = "# genuine stubs in Forward Endcap Disc";
  endcap_disc_Fw_genuine = iBooker.book1D(HistoName,
                                  HistoName,
                                  endcap_layers.getParameter<int32_t>("Nbinsx"),
                                  endcap_layers.getParameter<double>("xmin"),
                                  endcap_layers.getParameter<double>("xmax"));
  endcap_disc_Fw_genuine->setAxisTitle("forward endcap disc", 1);
  endcap_disc_Fw_genuine->setAxisTitle("# of genuine stubs", 2);

  // stub in backend endcap disc
  HistoName = "# genuine stubs in Backward Endcap Disc";
  endcap_disc_Bw_genuine = iBooker.book1D(HistoName,
                                  HistoName,
                                  endcap_layers.getParameter<int32_t>("Nbinsx"),
                                  endcap_layers.getParameter<double>("xmin"),
                                  endcap_layers.getParameter<double>("xmax"));
  endcap_disc_Bw_genuine->setAxisTitle("backward endcap disc", 1);
  endcap_disc_Bw_genuine->setAxisTitle("# of genuine stubs", 2);

  // Relative pT
  edm::ParameterSet psRes_ptRel = conf_.getParameter<edm::ParameterSet>("TH1Res_ptRel");
  HistoName = "res_ptRel";
  res_ptRel = iBooker.book1D(HistoName,
                             HistoName,
                             psRes_ptRel.getParameter<int32_t>("Nbinsx"),
                             psRes_ptRel.getParameter<double>("xmin"),
                             psRes_ptRel.getParameter<double>("xmax"));
  res_ptRel->setAxisTitle("Relative p_{T} [GeV]", 1);
  res_ptRel->setAxisTitle("# tracking particles", 2);

  // Eta parts (for resolution)
  // Eta 1 (0 to 0.7)
  HistoName = "reseta_eta0to0p7";
  reseta_eta0to0p7 = iBooker.book1D(HistoName,
                                    HistoName,
                                    psRes_eta.getParameter<int32_t>("Nbinsx"),
                                    psRes_eta.getParameter<double>("xmin"),
                                    psRes_eta.getParameter<double>("xmax"));
  reseta_eta0to0p7->setAxisTitle("#eta_{trk} - #eta_{tp}", 1);
  reseta_eta0to0p7->setAxisTitle("# tracking particles", 2);

  // Eta 2 (0.7 to 1.0)
  HistoName = "reseta_eta0p7to1";
  reseta_eta0p7to1 = iBooker.book1D(HistoName,
                                    HistoName,
                                    psRes_eta.getParameter<int32_t>("Nbinsx"),
                                    psRes_eta.getParameter<double>("xmin"),
                                    psRes_eta.getParameter<double>("xmax"));
  reseta_eta0p7to1->setAxisTitle("#eta_{trk} - #eta_{tp}", 1);
  reseta_eta0p7to1->setAxisTitle("# tracking particles", 2);

  // Eta 3 (1.0 to 1.2)
  HistoName = "reseta_eta1to1p2";
  reseta_eta1to1p2 = iBooker.book1D(HistoName,
                                    HistoName,
                                    psRes_eta.getParameter<int32_t>("Nbinsx"),
                                    psRes_eta.getParameter<double>("xmin"),
                                    psRes_eta.getParameter<double>("xmax"));
  reseta_eta1to1p2->setAxisTitle("#eta_{trk} - #eta_{tp}", 1);
  reseta_eta1to1p2->setAxisTitle("# tracking particles", 2);

  // Eta 4 (1.2 to 1.6)
  HistoName = "reseta_eta1p2to1p6";
  reseta_eta1p2to1p6 = iBooker.book1D(HistoName,
                                      HistoName,
                                      psRes_eta.getParameter<int32_t>("Nbinsx"),
                                      psRes_eta.getParameter<double>("xmin"),
                                      psRes_eta.getParameter<double>("xmax"));
  reseta_eta1p2to1p6->setAxisTitle("#eta_{trk} - #eta_{tp}", 1);
  reseta_eta1p2to1p6->setAxisTitle("# tracking particles", 2);

  // Eta 5 (1.6 to 2.0)
  HistoName = "reseta_eta1p6to2";
  reseta_eta1p6to2 = iBooker.book1D(HistoName,
                                    HistoName,
                                    psRes_eta.getParameter<int32_t>("Nbinsx"),
                                    psRes_eta.getParameter<double>("xmin"),
                                    psRes_eta.getParameter<double>("xmax"));
  reseta_eta1p6to2->setAxisTitle("#eta_{trk} - #eta_{tp}", 1);
  reseta_eta1p6to2->setAxisTitle("# tracking particles", 2);

  // Eta 6 (2.0 to 2.4)
  HistoName = "reseta_eta2to2p4";
  reseta_eta2to2p4 = iBooker.book1D(HistoName,
                                    HistoName,
                                    psRes_eta.getParameter<int32_t>("Nbinsx"),
                                    psRes_eta.getParameter<double>("xmin"),
                                    psRes_eta.getParameter<double>("xmax"));
  reseta_eta2to2p4->setAxisTitle("#eta_{trk} - #eta_{tp}", 1);
  reseta_eta2to2p4->setAxisTitle("# tracking particles", 2);

  // pT parts for resolution (pT res vs eta)
  // pT a (2 to 3 GeV)
  // Eta 1 (0 to 0.7)
  HistoName = "respt_eta0to0p7_pt2to3";
  respt_eta0to0p7_pt2to3 = iBooker.book1D(HistoName,
                                          HistoName,
                                          psRes_pt.getParameter<int32_t>("Nbinsx"),
                                          psRes_pt.getParameter<double>("xmin"),
                                          psRes_pt.getParameter<double>("xmax"));
  respt_eta0to0p7_pt2to3->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta0to0p7_pt2to3->setAxisTitle("# tracking particles", 2);

  // Eta 2 (0.7 to 1.0)
  HistoName = "respt_eta0p7to1_pt2to3";
  respt_eta0p7to1_pt2to3 = iBooker.book1D(HistoName,
                                          HistoName,
                                          psRes_pt.getParameter<int32_t>("Nbinsx"),
                                          psRes_pt.getParameter<double>("xmin"),
                                          psRes_pt.getParameter<double>("xmax"));
  respt_eta0p7to1_pt2to3->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta0p7to1_pt2to3->setAxisTitle("# tracking particles", 2);

  // Eta 3 (1.0 to 1.2)
  HistoName = "respt_eta1to1p2_pt2to3";
  respt_eta1to1p2_pt2to3 = iBooker.book1D(HistoName,
                                          HistoName,
                                          psRes_pt.getParameter<int32_t>("Nbinsx"),
                                          psRes_pt.getParameter<double>("xmin"),
                                          psRes_pt.getParameter<double>("xmax"));
  respt_eta1to1p2_pt2to3->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta1to1p2_pt2to3->setAxisTitle("# tracking particles", 2);

  // Eta 4 (1.2 to 1.6)
  HistoName = "respt_eta1p2to1p6_pt2to3";
  respt_eta1p2to1p6_pt2to3 = iBooker.book1D(HistoName,
                                            HistoName,
                                            psRes_pt.getParameter<int32_t>("Nbinsx"),
                                            psRes_pt.getParameter<double>("xmin"),
                                            psRes_pt.getParameter<double>("xmax"));
  respt_eta1p2to1p6_pt2to3->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta1p2to1p6_pt2to3->setAxisTitle("# tracking particles", 2);

  // Eta 5 (1.6 to 2.0)
  HistoName = "respt_eta1p6to2_pt2to3";
  respt_eta1p6to2_pt2to3 = iBooker.book1D(HistoName,
                                          HistoName,
                                          psRes_pt.getParameter<int32_t>("Nbinsx"),
                                          psRes_pt.getParameter<double>("xmin"),
                                          psRes_pt.getParameter<double>("xmax"));
  respt_eta1p6to2_pt2to3->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta1p6to2_pt2to3->setAxisTitle("# tracking particles", 2);

  // Eta 6 (2.0 to 2.4)
  HistoName = "respt_eta2to2p4_pt2to3";
  respt_eta2to2p4_pt2to3 = iBooker.book1D(HistoName,
                                          HistoName,
                                          psRes_pt.getParameter<int32_t>("Nbinsx"),
                                          psRes_pt.getParameter<double>("xmin"),
                                          psRes_pt.getParameter<double>("xmax"));
  respt_eta2to2p4_pt2to3->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta2to2p4_pt2to3->setAxisTitle("# tracking particles", 2);

  // pT b (3 to 8 GeV)
  // Eta 1 (0 to 0.7)
  HistoName = "respt_eta0to0p7_pt3to8";
  respt_eta0to0p7_pt3to8 = iBooker.book1D(HistoName,
                                          HistoName,
                                          psRes_pt.getParameter<int32_t>("Nbinsx"),
                                          psRes_pt.getParameter<double>("xmin"),
                                          psRes_pt.getParameter<double>("xmax"));
  respt_eta0to0p7_pt3to8->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta0to0p7_pt3to8->setAxisTitle("# tracking particles", 2);

  // Eta 2 (0.7 to 1.0)
  HistoName = "respt_eta0p7to1_pt3to8";
  respt_eta0p7to1_pt3to8 = iBooker.book1D(HistoName,
                                          HistoName,
                                          psRes_pt.getParameter<int32_t>("Nbinsx"),
                                          psRes_pt.getParameter<double>("xmin"),
                                          psRes_pt.getParameter<double>("xmax"));
  respt_eta0p7to1_pt3to8->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta0p7to1_pt3to8->setAxisTitle("# tracking particles", 2);

  // Eta 3 (1.0 to 1.2)
  HistoName = "respt_eta1to1p2_pt3to8";
  respt_eta1to1p2_pt3to8 = iBooker.book1D(HistoName,
                                          HistoName,
                                          psRes_pt.getParameter<int32_t>("Nbinsx"),
                                          psRes_pt.getParameter<double>("xmin"),
                                          psRes_pt.getParameter<double>("xmax"));
  respt_eta1to1p2_pt3to8->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta1to1p2_pt3to8->setAxisTitle("# tracking particles", 2);

  // Eta 4 (1.2 to 1.6)
  HistoName = "respt_eta1p2to1p6_pt3to8";
  respt_eta1p2to1p6_pt3to8 = iBooker.book1D(HistoName,
                                            HistoName,
                                            psRes_pt.getParameter<int32_t>("Nbinsx"),
                                            psRes_pt.getParameter<double>("xmin"),
                                            psRes_pt.getParameter<double>("xmax"));
  respt_eta1p2to1p6_pt3to8->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta1p2to1p6_pt3to8->setAxisTitle("# tracking particles", 2);

  // Eta 5 (1.6 to 2.0)
  HistoName = "respt_eta1p6to2_pt3to8";
  respt_eta1p6to2_pt3to8 = iBooker.book1D(HistoName,
                                          HistoName,
                                          psRes_pt.getParameter<int32_t>("Nbinsx"),
                                          psRes_pt.getParameter<double>("xmin"),
                                          psRes_pt.getParameter<double>("xmax"));
  respt_eta1p6to2_pt3to8->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta1p6to2_pt3to8->setAxisTitle("# tracking particles", 2);

  // Eta 6 (2.0 to 2.4)
  HistoName = "respt_eta2to2p4_pt3to8";
  respt_eta2to2p4_pt3to8 = iBooker.book1D(HistoName,
                                          HistoName,
                                          psRes_pt.getParameter<int32_t>("Nbinsx"),
                                          psRes_pt.getParameter<double>("xmin"),
                                          psRes_pt.getParameter<double>("xmax"));
  respt_eta2to2p4_pt3to8->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta2to2p4_pt3to8->setAxisTitle("# tracking particles", 2);

  // pT c (>8 GeV)
  // Eta 1 (0 to 0.7)
  HistoName = "respt_eta0to0p7_pt8toInf";
  respt_eta0to0p7_pt8toInf = iBooker.book1D(HistoName,
                                            HistoName,
                                            psRes_pt.getParameter<int32_t>("Nbinsx"),
                                            psRes_pt.getParameter<double>("xmin"),
                                            psRes_pt.getParameter<double>("xmax"));
  respt_eta0to0p7_pt8toInf->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta0to0p7_pt8toInf->setAxisTitle("# tracking particles", 2);

  // Eta 2 (0.7 to 1.0)
  HistoName = "respt_eta0p7to1_pt8toInf";
  respt_eta0p7to1_pt8toInf = iBooker.book1D(HistoName,
                                            HistoName,
                                            psRes_pt.getParameter<int32_t>("Nbinsx"),
                                            psRes_pt.getParameter<double>("xmin"),
                                            psRes_pt.getParameter<double>("xmax"));
  respt_eta0p7to1_pt8toInf->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta0p7to1_pt8toInf->setAxisTitle("# tracking particles", 2);

  // Eta 3 (1.0 to 1.2)
  HistoName = "respt_eta1to1p2_pt8toInf";
  respt_eta1to1p2_pt8toInf = iBooker.book1D(HistoName,
                                            HistoName,
                                            psRes_pt.getParameter<int32_t>("Nbinsx"),
                                            psRes_pt.getParameter<double>("xmin"),
                                            psRes_pt.getParameter<double>("xmax"));
  respt_eta1to1p2_pt8toInf->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta1to1p2_pt8toInf->setAxisTitle("# tracking particles", 2);

  // Eta 4 (1.2 to 1.6)
  HistoName = "respt_eta1p2to1p6_pt8toInf";
  respt_eta1p2to1p6_pt8toInf = iBooker.book1D(HistoName,
                                              HistoName,
                                              psRes_pt.getParameter<int32_t>("Nbinsx"),
                                              psRes_pt.getParameter<double>("xmin"),
                                              psRes_pt.getParameter<double>("xmax"));
  respt_eta1p2to1p6_pt8toInf->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta1p2to1p6_pt8toInf->setAxisTitle("# tracking particles", 2);

  // Eta 5 (1.6 to 2.0)
  HistoName = "respt_eta1p6to2_pt8toInf";
  respt_eta1p6to2_pt8toInf = iBooker.book1D(HistoName,
                                            HistoName,
                                            psRes_pt.getParameter<int32_t>("Nbinsx"),
                                            psRes_pt.getParameter<double>("xmin"),
                                            psRes_pt.getParameter<double>("xmax"));
  respt_eta1p6to2_pt8toInf->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta1p6to2_pt8toInf->setAxisTitle("# tracking particles", 2);

  // Eta 6 (2.0 to 2.4)
  HistoName = "respt_eta2to2p4_pt8toInf";
  respt_eta2to2p4_pt8toInf = iBooker.book1D(HistoName,
                                            HistoName,
                                            psRes_pt.getParameter<int32_t>("Nbinsx"),
                                            psRes_pt.getParameter<double>("xmin"),
                                            psRes_pt.getParameter<double>("xmax"));
  respt_eta2to2p4_pt8toInf->setAxisTitle("(p_{T}(trk) - p_{T}(tp))/p_{T}(tp)", 1);
  respt_eta2to2p4_pt8toInf->setAxisTitle("# tracking particles", 2);

  // Phi parts (for resolution)
  // Eta 1 (0 to 0.7)
  edm::ParameterSet psRes_phi = conf_.getParameter<edm::ParameterSet>("TH1Res_phi");
  HistoName = "resphi_eta0to0p7";
  resphi_eta0to0p7 = iBooker.book1D(HistoName,
                                    HistoName,
                                    psRes_phi.getParameter<int32_t>("Nbinsx"),
                                    psRes_phi.getParameter<double>("xmin"),
                                    psRes_phi.getParameter<double>("xmax"));
  resphi_eta0to0p7->setAxisTitle("#phi_{trk} - #phi_{tp}", 1);
  resphi_eta0to0p7->setAxisTitle("# tracking particles", 2);

  // Eta 2 (0.7 to 1.0)
  HistoName = "resphi_eta0p7to1";
  resphi_eta0p7to1 = iBooker.book1D(HistoName,
                                    HistoName,
                                    psRes_phi.getParameter<int32_t>("Nbinsx"),
                                    psRes_phi.getParameter<double>("xmin"),
                                    psRes_phi.getParameter<double>("xmax"));
  resphi_eta0p7to1->setAxisTitle("#phi_{trk} - #phi_{tp}", 1);
  resphi_eta0p7to1->setAxisTitle("# tracking particles", 2);

  // Eta 3 (1.0 to 1.2)
  HistoName = "resphi_eta1to1p2";
  resphi_eta1to1p2 = iBooker.book1D(HistoName,
                                    HistoName,
                                    psRes_phi.getParameter<int32_t>("Nbinsx"),
                                    psRes_phi.getParameter<double>("xmin"),
                                    psRes_phi.getParameter<double>("xmax"));
  resphi_eta1to1p2->setAxisTitle("#phi_{trk} - #phi_{tp}", 1);
  resphi_eta1to1p2->setAxisTitle("# tracking particles", 2);

  // Eta 4 (1.2 to 1.6)
  HistoName = "resphi_eta1p2to1p6";
  resphi_eta1p2to1p6 = iBooker.book1D(HistoName,
                                      HistoName,
                                      psRes_phi.getParameter<int32_t>("Nbinsx"),
                                      psRes_phi.getParameter<double>("xmin"),
                                      psRes_phi.getParameter<double>("xmax"));
  resphi_eta1p2to1p6->setAxisTitle("#phi_{trk} - #phi_{tp}", 1);
  resphi_eta1p2to1p6->setAxisTitle("# tracking particles", 2);

  // Eta 5 (1.6 to 2.0)
  HistoName = "resphi_eta1p6to2";
  resphi_eta1p6to2 = iBooker.book1D(HistoName,
                                    HistoName,
                                    psRes_phi.getParameter<int32_t>("Nbinsx"),
                                    psRes_phi.getParameter<double>("xmin"),
                                    psRes_phi.getParameter<double>("xmax"));
  resphi_eta1p6to2->setAxisTitle("#phi_{trk} - #phi_{tp}", 1);
  resphi_eta1p6to2->setAxisTitle("# tracking particles", 2);

  // Eta 6 (2.0 to 2.4)
  HistoName = "resphi_eta2to2p4";
  resphi_eta2to2p4 = iBooker.book1D(HistoName,
                                    HistoName,
                                    psRes_phi.getParameter<int32_t>("Nbinsx"),
                                    psRes_phi.getParameter<double>("xmin"),
                                    psRes_phi.getParameter<double>("xmax"));
  resphi_eta2to2p4->setAxisTitle("#phi_{trk} - #phi_{tp}", 1);
  resphi_eta2to2p4->setAxisTitle("# tracking particles", 2);

  // VtxZ parts (for resolution)
  // Eta 1 (0 to 0.7)
  edm::ParameterSet psRes_VtxZ = conf_.getParameter<edm::ParameterSet>("TH1Res_VtxZ");
  HistoName = "resVtxZ_eta0to0p7";
  resVtxZ_eta0to0p7 = iBooker.book1D(HistoName,
                                     HistoName,
                                     psRes_VtxZ.getParameter<int32_t>("Nbinsx"),
                                     psRes_VtxZ.getParameter<double>("xmin"),
                                     psRes_VtxZ.getParameter<double>("xmax"));
  resVtxZ_eta0to0p7->setAxisTitle("VtxZ_{trk} - VtxZ_{tp} [cm]", 1);
  resVtxZ_eta0to0p7->setAxisTitle("# tracking particles", 2);

  // Eta 2 (0.7 to 1.0)
  HistoName = "resVtxZ_eta0p7to1";
  resVtxZ_eta0p7to1 = iBooker.book1D(HistoName,
                                     HistoName,
                                     psRes_VtxZ.getParameter<int32_t>("Nbinsx"),
                                     psRes_VtxZ.getParameter<double>("xmin"),
                                     psRes_VtxZ.getParameter<double>("xmax"));
  resVtxZ_eta0p7to1->setAxisTitle("VtxZ_{trk} - VtxZ_{tp} [cm]", 1);
  resVtxZ_eta0p7to1->setAxisTitle("# tracking particles", 2);

  // Eta 3 (1.0 to 1.2)
  HistoName = "resVtxZ_eta1to1p2";
  resVtxZ_eta1to1p2 = iBooker.book1D(HistoName,
                                     HistoName,
                                     psRes_VtxZ.getParameter<int32_t>("Nbinsx"),
                                     psRes_VtxZ.getParameter<double>("xmin"),
                                     psRes_VtxZ.getParameter<double>("xmax"));
  resVtxZ_eta1to1p2->setAxisTitle("VtxZ_{trk} - VtxZ_{tp} [cm]", 1);
  resVtxZ_eta1to1p2->setAxisTitle("# tracking particles", 2);

  // Eta 4 (1.2 to 1.6)
  HistoName = "resVtxZ_eta1p2to1p6";
  resVtxZ_eta1p2to1p6 = iBooker.book1D(HistoName,
                                       HistoName,
                                       psRes_VtxZ.getParameter<int32_t>("Nbinsx"),
                                       psRes_VtxZ.getParameter<double>("xmin"),
                                       psRes_VtxZ.getParameter<double>("xmax"));
  resVtxZ_eta1p2to1p6->setAxisTitle("VtxZ_{trk} - VtxZ_{tp} [cm]", 1);
  resVtxZ_eta1p2to1p6->setAxisTitle("# tracking particles", 2);

  // Eta 5 (1.6 to 2.0)
  HistoName = "resVtxZ_eta1p6to2";
  resVtxZ_eta1p6to2 = iBooker.book1D(HistoName,
                                     HistoName,
                                     psRes_VtxZ.getParameter<int32_t>("Nbinsx"),
                                     psRes_VtxZ.getParameter<double>("xmin"),
                                     psRes_VtxZ.getParameter<double>("xmax"));
  resVtxZ_eta1p6to2->setAxisTitle("VtxZ_{trk} - VtxZ_{tp} [cm]", 1);
  resVtxZ_eta1p6to2->setAxisTitle("# tracking particles", 2);

  // Eta 6 (2.0 to 2.4)
  HistoName = "resVtxZ_eta2to2p4";
  resVtxZ_eta2to2p4 = iBooker.book1D(HistoName,
                                     HistoName,
                                     psRes_VtxZ.getParameter<int32_t>("Nbinsx"),
                                     psRes_VtxZ.getParameter<double>("xmin"),
                                     psRes_VtxZ.getParameter<double>("xmax"));
  resVtxZ_eta2to2p4->setAxisTitle("VtxZ_{trk} - VtxZ_{tp} [cm]", 1);
  resVtxZ_eta2to2p4->setAxisTitle("# tracking particles", 2);

  // d0 parts (for resolution)
  // Eta 1 (0 to 0.7)
  edm::ParameterSet psRes_d0 = conf_.getParameter<edm::ParameterSet>("TH1Res_d0");
  HistoName = "resd0_eta0to0p7";
  resd0_eta0to0p7 = iBooker.book1D(HistoName,
                                   HistoName,
                                   psRes_d0.getParameter<int32_t>("Nbinsx"),
                                   psRes_d0.getParameter<double>("xmin"),
                                   psRes_d0.getParameter<double>("xmax"));
  resd0_eta0to0p7->setAxisTitle("d0_{trk} - d0_{tp} [cm]", 1);
  resd0_eta0to0p7->setAxisTitle("# tracking particles", 2);

  // Eta 2 (0.7 to 1.0)
  HistoName = "resd0_eta0p7to1";
  resd0_eta0p7to1 = iBooker.book1D(HistoName,
                                   HistoName,
                                   psRes_d0.getParameter<int32_t>("Nbinsx"),
                                   psRes_d0.getParameter<double>("xmin"),
                                   psRes_d0.getParameter<double>("xmax"));
  resd0_eta0p7to1->setAxisTitle("d0_{trk} - d0_{tp} [cm]", 1);
  resd0_eta0p7to1->setAxisTitle("# tracking particles", 2);

  // Eta 3 (1.0 to 1.2)
  HistoName = "resd0_eta1to1p2";
  resd0_eta1to1p2 = iBooker.book1D(HistoName,
                                   HistoName,
                                   psRes_d0.getParameter<int32_t>("Nbinsx"),
                                   psRes_d0.getParameter<double>("xmin"),
                                   psRes_d0.getParameter<double>("xmax"));
  resd0_eta1to1p2->setAxisTitle("d0_{trk} - d0_{tp} [cm]", 1);
  resd0_eta1to1p2->setAxisTitle("# tracking particles", 2);

  // Eta 4 (1.2 to 1.6)
  HistoName = "resd0_eta1p2to1p6";
  resd0_eta1p2to1p6 = iBooker.book1D(HistoName,
                                     HistoName,
                                     psRes_d0.getParameter<int32_t>("Nbinsx"),
                                     psRes_d0.getParameter<double>("xmin"),
                                     psRes_d0.getParameter<double>("xmax"));
  resd0_eta1p2to1p6->setAxisTitle("d0_{trk} - d0_{tp} [cm]", 1);
  resd0_eta1p2to1p6->setAxisTitle("# tracking particles", 2);

  // Eta 5 (1.6 to 2.0)
  HistoName = "resd0_eta1p6to2";
  resd0_eta1p6to2 = iBooker.book1D(HistoName,
                                   HistoName,
                                   psRes_d0.getParameter<int32_t>("Nbinsx"),
                                   psRes_d0.getParameter<double>("xmin"),
                                   psRes_d0.getParameter<double>("xmax"));
  resd0_eta1p6to2->setAxisTitle("d0_{trk} - d0_{tp} [cm]", 1);
  resd0_eta1p6to2->setAxisTitle("# tracking particles", 2);

  // Eta 6 (2.0 to 2.4)
  HistoName = "resd0_eta2to2p4";
  resd0_eta2to2p4 = iBooker.book1D(HistoName,
                                   HistoName,
                                   psRes_d0.getParameter<int32_t>("Nbinsx"),
                                   psRes_d0.getParameter<double>("xmin"),
                                   psRes_d0.getParameter<double>("xmax"));
  resd0_eta2to2p4->setAxisTitle("d0_{trk} - d0_{tp} [cm]", 1);
  resd0_eta2to2p4->setAxisTitle("# tracking particles", 2);

  // 2D: trackBend vs. stubBend
  edm::ParameterSet pstiltAngleVsZ0 = conf_.getParameter<edm::ParameterSet>("TH2tiltAngleVsZ0");
  HistoName = "tiltAngle_vs_Z0";
  hist_tiltAngle_vs_Z0= iBooker.book2D(
      HistoName, 
      HistoName,
      pstiltAngleVsZ0.getParameter<int32_t>("Nbinsx"),
      pstiltAngleVsZ0.getParameter<double>("xmin"),
      pstiltAngleVsZ0.getParameter<double>("xmax"),
      pstiltAngleVsZ0.getParameter<int32_t>("Nbinsy"),
      pstiltAngleVsZ0.getParameter<double>("ymin"),
      pstiltAngleVsZ0.getParameter<double>("ymax"));
  hist_tiltAngle_vs_Z0->setAxisTitle("Z0 [cm]", 1);
  hist_tiltAngle_vs_Z0->setAxisTitle("tiltAngle [radians]", 2);

  edm::ParameterSet psDeltaRVsDeltaZ = conf_.getParameter<edm::ParameterSet>("TH2DeltaRVsDeltaZ");
  HistoName = "R0_vs_Z0";
  hist_deltaR_vs_deltaZ= iBooker.book2D(
      HistoName, 
      HistoName,
      psDeltaRVsDeltaZ.getParameter<int32_t>("Nbinsx"),
      psDeltaRVsDeltaZ.getParameter<double>("xmin"),
      psDeltaRVsDeltaZ.getParameter<double>("xmax"),
      psDeltaRVsDeltaZ.getParameter<int32_t>("Nbinsy"),
      psDeltaRVsDeltaZ.getParameter<double>("ymin"),
      psDeltaRVsDeltaZ.getParameter<double>("ymax"));
  hist_deltaR_vs_deltaZ->setAxisTitle("Z0 [cm]", 1);
  hist_deltaR_vs_deltaZ->setAxisTitle("R0 [cm]", 2);

  edm::ParameterSet psZ0VsDeltaZ = conf_.getParameter<edm::ParameterSet>("TH2Z0VsDeltaZ");
  HistoName = "Z0_vs_deltaZ";
  hist_Z0_vs_deltaZ = iBooker.book2D(
      HistoName, 
      HistoName,
      psZ0VsDeltaZ.getParameter<int32_t>("Nbinsx"),
      psZ0VsDeltaZ.getParameter<double>("xmin"),
      psZ0VsDeltaZ.getParameter<double>("xmax"),
      psZ0VsDeltaZ.getParameter<int32_t>("Nbinsy"),
      psZ0VsDeltaZ.getParameter<double>("ymin"),
      psZ0VsDeltaZ.getParameter<double>("ymax"));
  hist_Z0_vs_deltaZ->setAxisTitle("delta_Z [cm]", 1);
  hist_Z0_vs_deltaZ->setAxisTitle("Z0 [cm]", 2);

  edm::ParameterSet psR0VsDeltaR = conf_.getParameter<edm::ParameterSet>("TH2R0VsDeltaR");
  HistoName = "R0_vs_deltaR";
  hist_R0_vs_deltaR = iBooker.book2D(
      HistoName, 
      HistoName,
      psR0VsDeltaR.getParameter<int32_t>("Nbinsx"),
      psR0VsDeltaR.getParameter<double>("xmin"),
      psR0VsDeltaR.getParameter<double>("xmax"),
      psR0VsDeltaR.getParameter<int32_t>("Nbinsy"),
      psR0VsDeltaR.getParameter<double>("ymin"),
      psR0VsDeltaR.getParameter<double>("ymax"));
  hist_R0_vs_deltaR->setAxisTitle("delta_R [cm]", 1);
  hist_R0_vs_deltaR->setAxisTitle("R0 [cm]", 2);

  edm::ParameterSet psTpPhiVsStubPhi = conf_.getParameter<edm::ParameterSet>("TH2TpPhiVsStubPhi");
  HistoName = "trackPhi_vs_stubPhi";
  trackPhi_vs_stubPhi = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpPhiVsStubPhi.getParameter<int32_t>("Nbinsx"),
      psTpPhiVsStubPhi.getParameter<double>("xmin"),
      psTpPhiVsStubPhi.getParameter<double>("xmax"),
      psTpPhiVsStubPhi.getParameter<int32_t>("Nbinsy"),
      psTpPhiVsStubPhi.getParameter<double>("ymin"),
      psTpPhiVsStubPhi.getParameter<double>("ymax"));
  trackPhi_vs_stubPhi->setAxisTitle("track phi [radians]", 1);
  trackPhi_vs_stubPhi->setAxisTitle("stub phi [radians]", 2);

  HistoName = "trackPhi_vs_stubPhi_barrel";
  trackPhi_vs_stubPhi_barrel = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpPhiVsStubPhi.getParameter<int32_t>("Nbinsx"),
      psTpPhiVsStubPhi.getParameter<double>("xmin"),
      psTpPhiVsStubPhi.getParameter<double>("xmax"),
      psTpPhiVsStubPhi.getParameter<int32_t>("Nbinsy"),
      psTpPhiVsStubPhi.getParameter<double>("ymin"),
      psTpPhiVsStubPhi.getParameter<double>("ymax"));
  trackPhi_vs_stubPhi_barrel->setAxisTitle("track phi [radians]", 1);
  trackPhi_vs_stubPhi_barrel->setAxisTitle("stub phi [radians]", 2);

  HistoName = "trackPhi_vs_stubPhi_endcap";
  trackPhi_vs_stubPhi_endcap = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpPhiVsStubPhi.getParameter<int32_t>("Nbinsx"),
      psTpPhiVsStubPhi.getParameter<double>("xmin"),
      psTpPhiVsStubPhi.getParameter<double>("xmax"),
      psTpPhiVsStubPhi.getParameter<int32_t>("Nbinsy"),
      psTpPhiVsStubPhi.getParameter<double>("ymin"),
      psTpPhiVsStubPhi.getParameter<double>("ymax"));
  trackPhi_vs_stubPhi_endcap->setAxisTitle("track phi [radians]", 1);
  trackPhi_vs_stubPhi_endcap->setAxisTitle("stub phi [radians]", 2);

  edm::ParameterSet psTpZVsStubZ = conf_.getParameter<edm::ParameterSet>("TH2TpZVsStubZ");
  HistoName = "stub_maxZ_vs_minZ";
  stub_maxZ_vs_minZ = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_maxZ_vs_minZ->setAxisTitle("track z [cm]", 1);
  stub_maxZ_vs_minZ->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "modMaxZ_vs_modMinZ";
  modMaxZ_vs_modMinZ = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  modMaxZ_vs_modMinZ->setAxisTitle("track z [cm]", 1);
  modMaxZ_vs_modMinZ->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ";
  stub_Z_vs_tpZ = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_maxZ_vs_minZ_highPt";
  stub_maxZ_vs_minZ_highPt = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_maxZ_vs_minZ_highPt->setAxisTitle("track z [cm]", 1);
  stub_maxZ_vs_minZ_highPt->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L1";
  stub_Z_vs_tpZ_L1 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L1->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L1->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L1_tilted";
  stub_Z_vs_tpZ_L1_tilted = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L1_tilted->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L1_tilted->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L1_tilted_pT_2to3";
  stub_Z_vs_tpZ_L1_tilted_pT_2to3 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L1_tilted_pT_2to3->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L1_tilted_pT_2to3->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L1_tilted_pT_3to5";
  stub_Z_vs_tpZ_L1_tilted_pT_3to5 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L1_tilted_pT_3to5->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L1_tilted_pT_3to5->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L1_tilted_pT_5to10";
  stub_Z_vs_tpZ_L1_tilted_pT_5to10 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L1_tilted_pT_5to10->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L1_tilted_pT_5to10->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L2";
  stub_Z_vs_tpZ_L2 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L2->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L2->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L2_tilted";
  stub_Z_vs_tpZ_L2_tilted = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L2_tilted->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L2_tilted->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L2_tilted_pT_2to3";
  stub_Z_vs_tpZ_L2_tilted_pT_2to3 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L2_tilted_pT_2to3->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L2_tilted_pT_2to3->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L2_tilted_pT_3to5";
  stub_Z_vs_tpZ_L2_tilted_pT_3to5 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L2_tilted_pT_3to5->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L2_tilted_pT_3to5->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L2_tilted_pT_5to10";
  stub_Z_vs_tpZ_L2_tilted_pT_5to10 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L2_tilted_pT_5to10->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L2_tilted_pT_5to10->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L3";
  stub_Z_vs_tpZ_L3 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L3->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L3->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L3_tilted";
  stub_Z_vs_tpZ_L3_tilted = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L3_tilted->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L3_tilted->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L3_tilted_pT_2to3";
  stub_Z_vs_tpZ_L3_tilted_pT_2to3 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L3_tilted_pT_2to3->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L3_tilted_pT_2to3->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L3_tilted_pT_3to5";
  stub_Z_vs_tpZ_L3_tilted_pT_3to5 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L3_tilted_pT_3to5->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L3_tilted_pT_3to5->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L3_tilted_pT_5to10";
  stub_Z_vs_tpZ_L3_tilted_pT_5to10 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L3_tilted_pT_5to10->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L3_tilted_pT_5to10->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L4";
  stub_Z_vs_tpZ_L4 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L4->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L4->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L5";
  stub_Z_vs_tpZ_L5 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L5->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L5->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_L6";
  stub_Z_vs_tpZ_L6 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_L6->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_L6->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "modMaxZ_vs_modMinZ_highPt";
  modMaxZ_vs_modMinZ_highPt = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  modMaxZ_vs_modMinZ_highPt->setAxisTitle("track z [cm]", 1);
  modMaxZ_vs_modMinZ_highPt->setAxisTitle("stub z avg [cm]", 2);

  HistoName = "stub_Z_vs_tpZ_highPt";
  stub_Z_vs_tpZ_highPt = iBooker.book2D(
      HistoName, 
      HistoName,
      psTpZVsStubZ.getParameter<int32_t>("Nbinsx"),
      psTpZVsStubZ.getParameter<double>("xmin"),
      psTpZVsStubZ.getParameter<double>("xmax"),
      psTpZVsStubZ.getParameter<int32_t>("Nbinsy"),
      psTpZVsStubZ.getParameter<double>("ymin"),
      psTpZVsStubZ.getParameter<double>("ymax"));
  stub_Z_vs_tpZ_highPt->setAxisTitle("track z [cm]", 1);
  stub_Z_vs_tpZ_highPt->setAxisTitle("stub z avg [cm]", 2);

  edm::ParameterSet psTrackVsStub = conf_.getParameter<edm::ParameterSet>("TH2TrackVsStub");
  HistoName = "trackBend_vs_stubBend";
  trackBend_vs_stubBend = iBooker.book2D(
      HistoName, 
      HistoName,
      psTrackVsStub.getParameter<int32_t>("Nbinsx"),
      psTrackVsStub.getParameter<double>("xmin"),
      psTrackVsStub.getParameter<double>("xmax"),
      psTrackVsStub.getParameter<int32_t>("Nbinsy"),
      psTrackVsStub.getParameter<double>("ymin"),
      psTrackVsStub.getParameter<double>("ymax"));
  trackBend_vs_stubBend->setAxisTitle("Truth Particle Bend", 1);
  trackBend_vs_stubBend->setAxisTitle("Stub Bend", 2);

  // 2D: barrel trackBend vs. stubBend
  HistoName = "barrel_trackBend_vs_stubBend";
  barrel_trackBend_vs_stubBend = iBooker.book2D(
      HistoName, 
      HistoName,
      psTrackVsStub.getParameter<int32_t>("Nbinsx"),
      psTrackVsStub.getParameter<double>("xmin"),
      psTrackVsStub.getParameter<double>("xmax"),
      psTrackVsStub.getParameter<int32_t>("Nbinsy"),
      psTrackVsStub.getParameter<double>("ymin"),
      psTrackVsStub.getParameter<double>("ymax"));
  barrel_trackBend_vs_stubBend->setAxisTitle("Track Bend", 1);
  barrel_trackBend_vs_stubBend->setAxisTitle("Stub Bend", 2);

  HistoName = "barrel_trackBend_vs_stubBend_L1";
  barrel_trackBend_vs_stubBend_L1 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTrackVsStub.getParameter<int32_t>("Nbinsx"),
      psTrackVsStub.getParameter<double>("xmin"),
      psTrackVsStub.getParameter<double>("xmax"),
      psTrackVsStub.getParameter<int32_t>("Nbinsy"),
      psTrackVsStub.getParameter<double>("ymin"),
      psTrackVsStub.getParameter<double>("ymax"));
  barrel_trackBend_vs_stubBend_L1->setAxisTitle("Truth Particle Bend", 1);
  barrel_trackBend_vs_stubBend_L1->setAxisTitle("Stub Bend", 2);

  HistoName = "barrel_trackBend_vs_stubBend_L2";
  barrel_trackBend_vs_stubBend_L2 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTrackVsStub.getParameter<int32_t>("Nbinsx"),
      psTrackVsStub.getParameter<double>("xmin"),
      psTrackVsStub.getParameter<double>("xmax"),
      psTrackVsStub.getParameter<int32_t>("Nbinsy"),
      psTrackVsStub.getParameter<double>("ymin"),
      psTrackVsStub.getParameter<double>("ymax"));
  barrel_trackBend_vs_stubBend_L2->setAxisTitle("Truth Particle Bend", 1);
  barrel_trackBend_vs_stubBend_L2->setAxisTitle("Stub Bend", 2);

  HistoName = "barrel_trackBend_vs_stubBend_L3";
  barrel_trackBend_vs_stubBend_L3 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTrackVsStub.getParameter<int32_t>("Nbinsx"),
      psTrackVsStub.getParameter<double>("xmin"),
      psTrackVsStub.getParameter<double>("xmax"),
      psTrackVsStub.getParameter<int32_t>("Nbinsy"),
      psTrackVsStub.getParameter<double>("ymin"),
      psTrackVsStub.getParameter<double>("ymax"));
  barrel_trackBend_vs_stubBend_L3->setAxisTitle("Truth Particle Bend", 1);
  barrel_trackBend_vs_stubBend_L3->setAxisTitle("Stub Bend", 2);

  HistoName = "barrel_trackBend_vs_stubBend_L4";
  barrel_trackBend_vs_stubBend_L4 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTrackVsStub.getParameter<int32_t>("Nbinsx"),
      psTrackVsStub.getParameter<double>("xmin"),
      psTrackVsStub.getParameter<double>("xmax"),
      psTrackVsStub.getParameter<int32_t>("Nbinsy"),
      psTrackVsStub.getParameter<double>("ymin"),
      psTrackVsStub.getParameter<double>("ymax"));
  barrel_trackBend_vs_stubBend_L4->setAxisTitle("Truth Particle Bend", 1);
  barrel_trackBend_vs_stubBend_L4->setAxisTitle("Stub Bend", 2);

  HistoName = "barrel_trackBend_vs_stubBend_L5";
  barrel_trackBend_vs_stubBend_L5 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTrackVsStub.getParameter<int32_t>("Nbinsx"),
      psTrackVsStub.getParameter<double>("xmin"),
      psTrackVsStub.getParameter<double>("xmax"),
      psTrackVsStub.getParameter<int32_t>("Nbinsy"),
      psTrackVsStub.getParameter<double>("ymin"),
      psTrackVsStub.getParameter<double>("ymax"));
  barrel_trackBend_vs_stubBend_L5->setAxisTitle("Truth Particle Bend", 1);
  barrel_trackBend_vs_stubBend_L5->setAxisTitle("Stub Bend", 2);

  HistoName = "barrel_trackBend_vs_stubBend_L6";
  barrel_trackBend_vs_stubBend_L6 = iBooker.book2D(
      HistoName, 
      HistoName,
      psTrackVsStub.getParameter<int32_t>("Nbinsx"),
      psTrackVsStub.getParameter<double>("xmin"),
      psTrackVsStub.getParameter<double>("xmax"),
      psTrackVsStub.getParameter<int32_t>("Nbinsy"),
      psTrackVsStub.getParameter<double>("ymin"),
      psTrackVsStub.getParameter<double>("ymax"));
  barrel_trackBend_vs_stubBend_L6->setAxisTitle("Truth Particle Bend", 1);
  barrel_trackBend_vs_stubBend_L6->setAxisTitle("Stub Bend", 2);

  // 2D: endcap trackBend vs. stubBend
  HistoName = "endcap trackBend_vs_stubBend";
  endcap_trackBend_vs_stubBend = iBooker.book2D(
      HistoName, 
      HistoName,
      psTrackVsStub.getParameter<int32_t>("Nbinsx"),
      psTrackVsStub.getParameter<double>("xmin"),
      psTrackVsStub.getParameter<double>("xmax"),
      psTrackVsStub.getParameter<int32_t>("Nbinsy"),
      psTrackVsStub.getParameter<double>("ymin"),
      psTrackVsStub.getParameter<double>("ymax"));
  endcap_trackBend_vs_stubBend->setAxisTitle("Track Bend", 1);
  endcap_trackBend_vs_stubBend->setAxisTitle("Stub Bend", 2);

  // 2D: endcap fw trackBend vs. stubBend
  HistoName = "forward endcap trackBend_vs_stubBend";
  endcap_fw_trackBend_vs_stubBend = iBooker.book2D(
      HistoName, 
      HistoName,
      psTrackVsStub.getParameter<int32_t>("Nbinsx"),
      psTrackVsStub.getParameter<double>("xmin"),
      psTrackVsStub.getParameter<double>("xmax"),
      psTrackVsStub.getParameter<int32_t>("Nbinsy"),
      psTrackVsStub.getParameter<double>("ymin"),
      psTrackVsStub.getParameter<double>("ymax"));
  endcap_fw_trackBend_vs_stubBend->setAxisTitle("Track Bend", 1);
  endcap_fw_trackBend_vs_stubBend->setAxisTitle("Stub Bend", 2);

  // 2D: endcap fw trackBend vs. stubBend
  HistoName = "backward endcap trackBend_vs_stubBend";
  endcap_bw_trackBend_vs_stubBend = iBooker.book2D(
      HistoName, 
      HistoName,
      psTrackVsStub.getParameter<int32_t>("Nbinsx"),
      psTrackVsStub.getParameter<double>("xmin"),
      psTrackVsStub.getParameter<double>("xmax"),
      psTrackVsStub.getParameter<int32_t>("Nbinsy"),
      psTrackVsStub.getParameter<double>("ymin"),
      psTrackVsStub.getParameter<double>("ymax"));
  endcap_bw_trackBend_vs_stubBend->setAxisTitle("Track Bend", 1);
  endcap_bw_trackBend_vs_stubBend->setAxisTitle("Stub Bend", 2);

  edm::ParameterSet psMaxRvsMinR = conf_.getParameter<edm::ParameterSet>("TH2MaxRvsMinR");
  HistoName = "modMaxR_vs_modMinR";
  modMaxR_vs_modMinR = iBooker.book2D(
      HistoName, 
      HistoName,
      psMaxRvsMinR.getParameter<int32_t>("Nbinsx"),
      psMaxRvsMinR.getParameter<double>("xmin"),
      psMaxRvsMinR.getParameter<double>("xmax"),
      psMaxRvsMinR.getParameter<int32_t>("Nbinsy"),
      psMaxRvsMinR.getParameter<double>("ymin"),
      psMaxRvsMinR.getParameter<double>("ymax"));
  modMaxR_vs_modMinR->setAxisTitle("track z [cm]", 1);
  modMaxR_vs_modMinR->setAxisTitle("stub z avg [cm]", 2);



}  // end of method

DEFINE_FWK_MODULE(OuterTrackerMonitorTrackingParticles);