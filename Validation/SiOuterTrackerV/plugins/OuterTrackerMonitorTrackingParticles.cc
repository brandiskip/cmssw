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
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"

#include "OuterTrackerMonitorTrackingParticles.h"

#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"

//
// constructors and destructor
//
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
  getTokenTrackerGeom_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();
  edm::InputTag L1StubInputTag = conf_.getParameter<edm::InputTag>("L1StubInputTag");
  ttStubToken_ = consumes<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > >(L1StubInputTag);
  L1Tk_minNStub = conf_.getParameter<int>("L1Tk_minNStub");         // min number of stubs in the track
  L1Tk_maxChi2dof = conf_.getParameter<double>("L1Tk_maxChi2dof");  // maximum chi2/dof of the track
  TP_minNStub = conf_.getParameter<int>("TP_minNStub");             // min number of stubs in the tracking particle to
  //min number of layers with stubs in the tracking particle to consider matching
  TP_minNLayersStub = conf_.getParameter<int>("TP_minNLayersStub");
  TP_minPt = conf_.getParameter<double>("TP_minPt");      // min pT to consider matching
  TP_maxEta = conf_.getParameter<double>("TP_maxEta");    // max eta to consider matching
  TP_maxVtxZ = conf_.getParameter<double>("TP_maxVtxZ");  // max vertZ (or z0) to consider matching
}

OuterTrackerMonitorTrackingParticles::~OuterTrackerMonitorTrackingParticles() = default;

// member functions

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

  //-------------------------------------------------------------------------------------------------
  // declare L1StubInputTag as an edm::InputTag variable 
  // create a placeholder for an input tag that will be used to specify the collection of data to be processed
  edm::InputTag L1StubInputTag("TTStubsFromPhase2TrackerDigis","StubAccepted");

  // retrieves a configuration parameter named "L1StubInputTag" of type edm::InputTag from the iConfig object and assigns its value to the L1StubInputTag variable
  // this used to read iConfig.getParameter but I had to change it to conf_.getParameter because of the constructor (need to check out exactly why)
  // L1StubInputTag = conf_.getParameter<edm::InputTag>("L1StubInputTag");

  // L1 Stubs
  edm::Handle<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > > TTStubHandle;
  iEvent.getByToken(ttStubToken_, TTStubHandle);

  // more for TTStubs
  edm::ESHandle<TrackerGeometry> tGeomHandle = iSetup.getHandle(getTokenTrackerGeom_);
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();

  // Geometries
  const TrackerTopology *const tTopo = &iSetup.getData(m_topoToken);

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
      float matchStub_eta = -999;
      float myTP_phi = -999;

      //-------------------------------------------------------------------------------------------------      
      // Loop over L1 stubs
      //-------------------------------------------------------------------------------------------------      
      for (auto gd = theTrackerGeom->dets().begin(); gd != theTrackerGeom->dets().end(); gd++) {
        DetId detid = (*gd)->geographicalId();
        if (detid.subdetId() != StripSubdetector::TOB && detid.subdetId() != StripSubdetector::TID)
          continue;
        if (!tTopo->isLower(detid))
          continue;                              // loop on the stacks: choose the lower arbitrarily
        DetId stackDetid = tTopo->stack(detid);  // Stub module detid

        if (TTStubHandle->find(stackDetid) == TTStubHandle->end())
          continue;

        // Get the DetSets of the Clusters
        edmNew::DetSet<TTStub<Ref_Phase2TrackerDigi_> > stubs = (*TTStubHandle)[stackDetid];
      
        // loop over stubs
        for (auto stubIter = stubs.begin(); stubIter != stubs.end(); ++stubIter) {
          edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> >, TTStub<Ref_Phase2TrackerDigi_> > tempStubPtr = edmNew::makeRefTo(TTStubHandle, stubIter);
  
          // matched stub to tracking particle? (access TrackingParticle object)
          // find the associated TrackingParticle corresponding to a given stub
          edm::Ptr<TrackingParticle> match_stub = MCTruthTTStubHandle->findTrackingParticlePtr(tempStubPtr);
            
          // if the following is not null, it means that a valid associated tracking particle was found for the stub
          if (match_stub.isNull() == false) {
            // retrieve the event ID (event number) associated with the tracking particle
            int tmp_eventid = match_stub->eventId().event();

            if (tmp_eventid > 0)
              continue;  // this means stub from pileup track

            matchStub_eta = match_stub->p4().eta();
            myTP_phi = match_stub->p4().phi();

            float stub_eta_res = matchStub_eta - tmp_tp_eta;
            float stub_phi_res = myTP_phi - tmp_tp_phi;

            stub_res_eta->Fill(stub_eta_res);
            stub_res_phi->Fill(stub_phi_res);
            check_matchStub_eta->Fill(matchStub_eta);
            check_tmp_tp_eta->Fill(tmp_tp_eta);
                        
            std::cout << "matchStub_eta value: " << matchStub_eta << std::endl;
            std::cout << "tmp_tp_eta value: " << tmp_tp_eta << std::endl;

            // Fill resolution plots for different abs(eta) bins:
            // (0, 0.4), (0.4, 0.7), (0.7, 1.0), (1.0, 1.2), (1.2, 1.6), (1.6, 2.0), (2.0, 2.4)
            if (std::fabs(tmp_tp_eta) >= 0 && std::fabs(tmp_tp_eta) < 0.4) {
              stub_reseta_eta0to0p4->Fill(stub_eta_res);
              stub_resphi_eta0to0p4->Fill(stub_phi_res);
          } else if (std::fabs(tmp_tp_eta) >= 0.4 && std::fabs(tmp_tp_eta) < 0.7) {
              stub_reseta_eta0p4to0p7->Fill(stub_eta_res);
              stub_resphi_eta0p4to0p7->Fill(stub_phi_res); 
          } else if (std::fabs(tmp_tp_eta) >= 0.7 && std::fabs(tmp_tp_eta) < 1.0) {
              stub_reseta_eta0p7to1->Fill(stub_eta_res);
              stub_resphi_eta0p7to1->Fill(stub_phi_res); 
          } else if (std::fabs(tmp_tp_eta) >= 1.0 && std::fabs(tmp_tp_eta) < 1.2) {
              stub_reseta_eta1to1p2->Fill(stub_eta_res);
              stub_resphi_eta1to1p2->Fill(stub_phi_res);
          } else if (std::fabs(tmp_tp_eta) >= 1.2 && std::fabs(tmp_tp_eta) < 1.6) {
              stub_reseta_eta1p2to1p6->Fill(stub_eta_res);
              stub_resphi_eta1p2to1p6->Fill(stub_phi_res); 
          } else if (std::fabs(tmp_tp_eta) >= 1.6 && std::fabs(tmp_tp_eta) < 2.0) {
              stub_reseta_eta1p6to2->Fill(stub_eta_res);
              stub_resphi_eta1p6to2->Fill(stub_phi_res);
          } else if (std::fabs(tmp_tp_eta) >= 2.0 && std::fabs(tmp_tp_eta) < 2.4) {
              stub_reseta_eta2to2p4->Fill(stub_eta_res);
              stub_resphi_eta2to2p4->Fill(stub_phi_res);
            }      
          }
        } // end loop over stubs
      } // end loop over L1 stubs


//-------------------------------------------------------------------------------------------------

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

      std::cout << "tmp_tp_eta outside value: " << tmp_tp_eta << std::endl;
      // Eta and pT histograms for resolution
      float pt_diff = tmp_matchtrk_pt - tmp_tp_pt;
      float test_pt_diff = tmp_matchtrk_pt - tmp_tp_pt;
      float pt_res = pt_diff / tmp_tp_pt;
      float eta_res = tmp_matchtrk_eta - tmp_tp_eta;
      float phi_res = tmp_matchtrk_phi - tmp_tp_phi;
      float VtxZ_res = tmp_matchtrk_VtxZ - tmp_tp_VtxZ;
      float d0_res = tmp_matchtrk_d0 - tmp_tp_d0;

      // fill total resolution histograms
      res_pt->Fill(pt_diff);
      test_res_pt->Fill(test_pt_diff);
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
}  // end of method

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

  // test full pT
  edm::ParameterSet test_psRes_pt = conf_.getParameter<edm::ParameterSet>("test_TH1Res_pt");
  HistoName = "test_res_pt";
  test_res_pt = iBooker.book1D(HistoName,
                          HistoName,
                          test_psRes_pt.getParameter<int32_t>("Nbinsx"),
                          test_psRes_pt.getParameter<double>("xmin"),
                          test_psRes_pt.getParameter<double>("xmax"));
  test_res_pt->setAxisTitle("p_{T} [GeV]", 1);
  test_res_pt->setAxisTitle("# test_tracking particles", 2);

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

  // Full stub eta
  edm::ParameterSet psStubRes_eta = conf_.getParameter<edm::ParameterSet>("TH1StubRes_eta");
  HistoName = "stub_res_eta";
  stub_res_eta = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_eta.getParameter<int32_t>("Nbinsx"),
                           psStubRes_eta.getParameter<double>("xmin"),
                           psStubRes_eta.getParameter<double>("xmax"));
  stub_res_eta->setAxisTitle("#eta", 1);
  stub_res_eta->setAxisTitle("# tracking particles", 2);

  
  // Full stub phi
  edm::ParameterSet psStubRes_phi = conf_.getParameter<edm::ParameterSet>("TH1StubRes_phi");
  HistoName = "stub_res_phi";
  stub_res_phi = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_phi.getParameter<int32_t>("Nbinsx"),
                           psStubRes_phi.getParameter<double>("xmin"),
                           psStubRes_phi.getParameter<double>("xmax"));
  stub_res_phi->setAxisTitle("#phi", 1);
  stub_res_phi->setAxisTitle("# tracking particles", 2);

  // check to see what myTP_eta histo looks like
  edm::ParameterSet check_psRes_eta = conf_.getParameter<edm::ParameterSet>("check_TH1Res_eta");
  HistoName = "check_matchStub_eta";
  check_matchStub_eta = iBooker.book1D(HistoName,
                           HistoName,
                           check_psRes_eta.getParameter<int32_t>("Nbinsx"),
                           check_psRes_eta.getParameter<double>("xmin"),
                           check_psRes_eta.getParameter<double>("xmax"));
  check_matchStub_eta->setAxisTitle("#eta", 1);
  check_matchStub_eta->setAxisTitle("# tracking particles", 2);

  // check to see what tmp_tp_eta histo looks like
  edm::ParameterSet check_tmp_psRes_eta = conf_.getParameter<edm::ParameterSet>("check_tmp_TH1Res_eta");
  HistoName = "check_tmp_tp_eta";
  check_tmp_tp_eta = iBooker.book1D(HistoName,
                           HistoName,
                           check_psRes_eta.getParameter<int32_t>("Nbinsx"),
                           check_psRes_eta.getParameter<double>("xmin"),
                           check_psRes_eta.getParameter<double>("xmax"));
  check_tmp_tp_eta->setAxisTitle("#eta", 1);
  check_tmp_tp_eta->setAxisTitle("# tracking particles", 2);

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

  // Eta parts (for resolution)
  // Stub Eta 1 (0 to 0.4)
  HistoName = "stub_reseta_eta0to0p4";
  stub_reseta_eta0to0p4 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_eta.getParameter<int32_t>("Nbinsx"),
                           psStubRes_eta.getParameter<double>("xmin"),
                           psStubRes_eta.getParameter<double>("xmax"));
  stub_reseta_eta0to0p4->setAxisTitle("#eta_{trk} - #eta_{tp}", 1);
  stub_reseta_eta0to0p4->setAxisTitle("# tracking particles", 2);
  
  // Stub Eta 2 (0.4 to 0.7)
  HistoName = "stub_reseta_eta0p4to0p7";
  stub_reseta_eta0p4to0p7 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_eta.getParameter<int32_t>("Nbinsx"),
                           psStubRes_eta.getParameter<double>("xmin"),
                           psStubRes_eta.getParameter<double>("xmax"));
  stub_reseta_eta0p4to0p7->setAxisTitle("#eta_{trk} - #eta_{tp}", 1);
  stub_reseta_eta0p4to0p7->setAxisTitle("# tracking particles", 2);

  // Stub Eta 3 (0.7 to 1.0)
  HistoName = "stub_reseta_eta0p7to1";
  stub_reseta_eta0p7to1 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_eta.getParameter<int32_t>("Nbinsx"),
                           psStubRes_eta.getParameter<double>("xmin"),
                           psStubRes_eta.getParameter<double>("xmax"));
  stub_reseta_eta0p7to1->setAxisTitle("#eta_{trk} - #eta_{tp}", 1);
  stub_reseta_eta0p7to1->setAxisTitle("# tracking particles", 2);

  // Stub Eta 4 (1.0 to 1.2)
  HistoName = "stub_reseta_eta1to1p2";
  stub_reseta_eta1to1p2 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_eta.getParameter<int32_t>("Nbinsx"),
                           psStubRes_eta.getParameter<double>("xmin"),
                           psStubRes_eta.getParameter<double>("xmax"));
  stub_reseta_eta1to1p2->setAxisTitle("#eta_{trk} - #eta_{tp}", 1);
  stub_reseta_eta1to1p2->setAxisTitle("# tracking particles", 2);

  // Stub Eta 5 (1.2 to 1.6)
  HistoName = "stub_reseta_eta1p2to1p6";
  stub_reseta_eta1p2to1p6 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_eta.getParameter<int32_t>("Nbinsx"),
                           psStubRes_eta.getParameter<double>("xmin"),
                           psStubRes_eta.getParameter<double>("xmax"));
  stub_reseta_eta1p2to1p6->setAxisTitle("#eta_{trk} - #eta_{tp}", 1);
  stub_reseta_eta1p2to1p6->setAxisTitle("# tracking particles", 2);

  // Stub Eta 6 (1.2 to 2.0)
  HistoName = "stub_reseta_eta1p6to2";
  stub_reseta_eta1p6to2 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_eta.getParameter<int32_t>("Nbinsx"),
                           psStubRes_eta.getParameter<double>("xmin"),
                           psStubRes_eta.getParameter<double>("xmax"));
  stub_reseta_eta1p6to2->setAxisTitle("#eta_{trk} - #eta_{tp}", 1);
  stub_reseta_eta1p6to2->setAxisTitle("# tracking particles", 2);

  // Stub Eta 7 (2.0 to 2.4)
  HistoName = "stub_reseta_eta2to2p4";
  stub_reseta_eta2to2p4 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_eta.getParameter<int32_t>("Nbinsx"),
                           psStubRes_eta.getParameter<double>("xmin"),
                           psStubRes_eta.getParameter<double>("xmax"));
  stub_reseta_eta2to2p4->setAxisTitle("#eta_{trk} - #eta_{tp}", 1);
  stub_reseta_eta2to2p4->setAxisTitle("# tracking particles", 2);

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

  // Phi parts (for resolution)
  // Eta 1 (0 to 0.4)
  HistoName = "stub_resphi_eta0to0p4";
  stub_resphi_eta0to0p4 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_phi.getParameter<int32_t>("Nbinsx"),
                           psStubRes_phi.getParameter<double>("xmin"),
                           psStubRes_phi.getParameter<double>("xmax"));
  stub_resphi_eta0to0p4->setAxisTitle("#phi", 1);
  stub_resphi_eta0to0p4->setAxisTitle("# tracking particles", 2);

  // Eta 2 (0.4 to 0.7)
  HistoName = "stub_resphi_eta0p4to0p7";
  stub_resphi_eta0p4to0p7 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_phi.getParameter<int32_t>("Nbinsx"),
                           psStubRes_phi.getParameter<double>("xmin"),
                           psStubRes_phi.getParameter<double>("xmax"));
  stub_resphi_eta0p4to0p7->setAxisTitle("#phi", 1);
  stub_resphi_eta0p4to0p7->setAxisTitle("# tracking particles", 2);

  // Eta 3 (0.7 to 1.0)
  HistoName = "stub_resphi_eta0p7to1";
  stub_resphi_eta0p7to1 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_phi.getParameter<int32_t>("Nbinsx"),
                           psStubRes_phi.getParameter<double>("xmin"),
                           psStubRes_phi.getParameter<double>("xmax"));
  stub_resphi_eta0p7to1->setAxisTitle("#phi", 1);
  stub_resphi_eta0p7to1->setAxisTitle("# tracking particles", 2);

  // Eta 4 (1.0 to 1.2)
  HistoName = "stub_resphi_eta1to1p2";
  stub_resphi_eta1to1p2 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_phi.getParameter<int32_t>("Nbinsx"),
                           psStubRes_phi.getParameter<double>("xmin"),
                           psStubRes_phi.getParameter<double>("xmax"));
  stub_resphi_eta1to1p2->setAxisTitle("#phi", 1);
  stub_resphi_eta1to1p2->setAxisTitle("# tracking particles", 2);

  // Eta 5 (1.2 to 1.6)
  HistoName = "stub_resphi_eta1p2to1p6";
  stub_resphi_eta1p2to1p6 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_phi.getParameter<int32_t>("Nbinsx"),
                           psStubRes_phi.getParameter<double>("xmin"),
                           psStubRes_phi.getParameter<double>("xmax"));
  stub_resphi_eta1p2to1p6->setAxisTitle("#phi", 1);
  stub_resphi_eta1p2to1p6->setAxisTitle("# tracking particles", 2);

  // Eta 6 (1.6 to 2.0)
  HistoName = "stub_resphi_eta1p6to2";
  stub_resphi_eta1p6to2 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_phi.getParameter<int32_t>("Nbinsx"),
                           psStubRes_phi.getParameter<double>("xmin"),
                           psStubRes_phi.getParameter<double>("xmax"));
  stub_resphi_eta1p6to2->setAxisTitle("#phi", 1);
  stub_resphi_eta1p6to2->setAxisTitle("# tracking particles", 2);

  // Eta 6 (2.0 to 2.4)
  HistoName = "stub_resphi_eta2to2p4";
  stub_resphi_eta2to2p4 = iBooker.book1D(HistoName,
                           HistoName,
                           psStubRes_phi.getParameter<int32_t>("Nbinsx"),
                           psStubRes_phi.getParameter<double>("xmin"),
                           psStubRes_phi.getParameter<double>("xmax"));
  stub_resphi_eta2to2p4->setAxisTitle("#phi", 1);
  stub_resphi_eta2to2p4->setAxisTitle("# tracking particles", 2);

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

}  // end of method

DEFINE_FWK_MODULE(OuterTrackerMonitorTrackingParticles);
