#ifndef SiOuterTrackerV_OuterTrackerMonitorTrackingParticles_h
#define SiOuterTrackerV_OuterTrackerMonitorTrackingParticles_h

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include <memory>
#include <string>
#include <vector>

class OuterTrackerMonitorTrackingParticles : public DQMEDAnalyzer {
public:
  explicit OuterTrackerMonitorTrackingParticles(const edm::ParameterSet &);
  //explicit OuterTrackerMonitorTrackingParticles(const edm::ParameterSet &, const trklet::Settings&);
  ~OuterTrackerMonitorTrackingParticles() override;
  void analyze(const edm::Event &, const edm::EventSetup &) override;
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  float phiOverBendCorrection(bool, float, float, const TrackerTopology*, uint32_t, const GeomDetUnit*, const GeomDetUnit*);
  std::vector<double> getTPDerivedCoords(edm::Ptr<TrackingParticle> my_tp, bool, double modMaxZ, double modMinZ, float stub_r) const;


  // Number of stubs
  MonitorElement *Stub_Barrel = nullptr;       // TTStub per layer
  MonitorElement *hist_stub_z = nullptr;
  MonitorElement *hist_stub_z_barrel = nullptr;
  MonitorElement *hist_stub_z_endcap = nullptr;

  // 1D correction factor
  MonitorElement *hist_deltaZ = nullptr;
  MonitorElement *hist_deltaR = nullptr;
  MonitorElement *hist_tiltAngle = nullptr;
  MonitorElement *hist_tp_phi = nullptr;
  MonitorElement *hist_stub_phi = nullptr;
  MonitorElement *hist_cosTiltAngle = nullptr;
  MonitorElement *hist_sinTiltAngle = nullptr;

  // fake rate
  MonitorElement *TotalStubs = nullptr;
  MonitorElement *FakeStubs = nullptr;

  // 2D correction factor
  MonitorElement *hist_tiltAngle_vs_Z0 = nullptr;
  MonitorElement *hist_deltaR_vs_deltaZ =nullptr;
  MonitorElement *hist_Z0_vs_deltaZ = nullptr;
  MonitorElement *hist_R0_vs_deltaR = nullptr;

  // Tracking particle distributions
  MonitorElement *trackParts_Eta = nullptr;
  MonitorElement *trackParts_Phi = nullptr;
  MonitorElement *trackParts_Pt = nullptr;
  MonitorElement *TP_z0 = nullptr;
  MonitorElement *hist_tp_z = nullptr;
  MonitorElement *hist_tp_z_barrel = nullptr;
  MonitorElement *hist_tp_z_endcap = nullptr;

  // pT and eta for efficiency plots
  MonitorElement *gen_clusters = nullptr;               // denominator
  MonitorElement *gen_clusters_zoom = nullptr;          // denominator
  MonitorElement *tp_pt = nullptr;                      // denominator
  MonitorElement *tp_pt_zoom = nullptr;                 // denominator
  MonitorElement *tp_eta = nullptr;                     // denominator
  MonitorElement *tp_d0 = nullptr;                      // denominator
  MonitorElement *tp_VtxR = nullptr;                    // denominator (also known as vxy)
  MonitorElement *tp_VtxZ = nullptr;                    // denominator
  MonitorElement *gen_clusters_if_stub = nullptr;       // numerator
  MonitorElement *gen_clusters_if_stub_zoom = nullptr;  // numerator
  MonitorElement *match_tp_pt = nullptr;                // numerator
  MonitorElement *match_tp_pt_zoom = nullptr;           // numerator
  MonitorElement *match_tp_eta = nullptr;               // numerator
  MonitorElement *match_tp_d0 = nullptr;                // numerator
  MonitorElement *match_tp_VtxR = nullptr;              // numerator (also known as vxy)
  MonitorElement *match_tp_VtxZ = nullptr;              // numerator

  // 2D plots
  MonitorElement *trackPhi_vs_stubPhi = nullptr;
  MonitorElement *trackPhi_vs_stubPhi_barrel = nullptr;
  MonitorElement *trackPhi_vs_stubPhi_endcap = nullptr;
  MonitorElement *trackBend_vs_stubBend = nullptr;
  MonitorElement *coordsBC = nullptr;
  MonitorElement *stub_maxZ_vs_minZ = nullptr;
  MonitorElement *modMaxZ_vs_modMinZ = nullptr;
  MonitorElement *stub_Z_vs_tpZ = nullptr;
  MonitorElement *stub_maxZ_vs_minZ_highPt = nullptr;
  MonitorElement *modMaxZ_vs_modMinZ_highPt = nullptr;
  MonitorElement *stub_Z_vs_tpZ_highPt = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L1 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L1_tilted = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L1_tilted_pT_2to3 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L1_tilted_pT_3to5 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L1_tilted_pT_5to10 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L2 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L2_tilted = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L2_tilted_pT_2to3 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L2_tilted_pT_3to5 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L2_tilted_pT_5to10 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L3 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L3_tilted = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L3_tilted_pT_2to3 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L3_tilted_pT_3to5 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L3_tilted_pT_5to10 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L4 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L5 = nullptr;
  MonitorElement *stub_Z_vs_tpZ_L6 = nullptr;
  MonitorElement *modMaxR_vs_modMinR = nullptr;
  MonitorElement *barrel_trackBend_vs_stubBend = nullptr;
  MonitorElement *barrel_trackBend_vs_stubBend_L1 = nullptr;
  MonitorElement *barrel_trackBend_vs_stubBend_L2 = nullptr;
  MonitorElement *barrel_trackBend_vs_stubBend_L3 = nullptr;
  MonitorElement *barrel_trackBend_vs_stubBend_L4 = nullptr;
  MonitorElement *barrel_trackBend_vs_stubBend_L5 = nullptr;
  MonitorElement *barrel_trackBend_vs_stubBend_L6 = nullptr;
  MonitorElement *endcap_trackBend_vs_stubBend = nullptr;
  MonitorElement *endcap_fw_trackBend_vs_stubBend = nullptr;
  MonitorElement *endcap_bw_trackBend_vs_stubBend = nullptr;
 
  // 1D stub bend resolution plots
  MonitorElement *bend_res = nullptr;
  MonitorElement *bend_res_barrel = nullptr;
  MonitorElement *bend_res_barrel_L1 = nullptr;
  MonitorElement *bend_res_barrel_L2 = nullptr;
  MonitorElement *bend_res_barrel_L3 = nullptr;
  MonitorElement *bend_res_barrel_L4 = nullptr;
  MonitorElement *bend_res_barrel_L5 = nullptr;
  MonitorElement *bend_res_barrel_L6 = nullptr;
  MonitorElement *bend_res_endcap = nullptr;
  MonitorElement *bend_res_fw_endcap = nullptr;
  MonitorElement *bend_res_fw_endcap_D1 = nullptr;
  MonitorElement *bend_res_fw_endcap_D2 = nullptr;
  MonitorElement *bend_res_fw_endcap_D3 = nullptr;
  MonitorElement *bend_res_fw_endcap_D4 = nullptr;
  MonitorElement *bend_res_fw_endcap_D5 = nullptr;
  MonitorElement *bend_res_bw_endcap = nullptr;
  MonitorElement *bend_res_bw_endcap_D1 = nullptr;
  MonitorElement *bend_res_bw_endcap_D2 = nullptr;
  MonitorElement *bend_res_bw_endcap_D3 = nullptr;
  MonitorElement *bend_res_bw_endcap_D4 = nullptr;
  MonitorElement *bend_res_bw_endcap_D5 = nullptr;
  MonitorElement *stub_res_phi = nullptr;
  MonitorElement *stub_z_res = nullptr;
  MonitorElement *stub_z_res_g0 = nullptr;
  MonitorElement *stub_z_res_l0 = nullptr;
  MonitorElement *z_res_barrel = nullptr;
  MonitorElement *z_res_barrel_L1 = nullptr;
  MonitorElement *z_res_barrel_L2 = nullptr;
  MonitorElement *z_res_barrel_L3 = nullptr;
  MonitorElement *z_res_barrel_L4 = nullptr;
  MonitorElement *z_res_barrel_L5 = nullptr;
  MonitorElement *z_res_barrel_L6 = nullptr;
  MonitorElement *z_res_isPS = nullptr;
  MonitorElement *z_res_isPS_barrel = nullptr;
  MonitorElement *z_res_isPS_endcap = nullptr;
  MonitorElement *z_res_is2S = nullptr;
  MonitorElement *z_res_is2S_barrel = nullptr;
  MonitorElement *z_res_is2S_endcap = nullptr;
  MonitorElement *z_res_endcap = nullptr;
  MonitorElement *stub_phi_res = nullptr;
  MonitorElement *stub_phi_res_barrel = nullptr;
  MonitorElement *stub_phi_res_isPS = nullptr;
  MonitorElement *stub_phi_res_isPS_barrel = nullptr;
  MonitorElement *stub_phi_res_isPS_endcap = nullptr;
  MonitorElement *stub_phi_res_is2S = nullptr;
  MonitorElement *stub_phi_res_is2S_barrel = nullptr;
  MonitorElement *stub_phi_res_is2S_endcap = nullptr;
  MonitorElement *stub_phi_res_endcap = nullptr;
  
  // 1D stub and associated tp plots
  MonitorElement *barrelHistogram_genuine = nullptr;
  MonitorElement *endcapHistogram_genuine = nullptr;
  MonitorElement *endcap_disc_Fw_genuine = nullptr;
  MonitorElement *endcap_disc_Bw_genuine = nullptr;
  MonitorElement *stub_rawBend = nullptr;
  MonitorElement *stub_bendOffset = nullptr;
  MonitorElement *stub_inClusPos = nullptr;
  MonitorElement *stub_bendFE = nullptr;
  MonitorElement *track_bend = nullptr;
  MonitorElement *numOfStubs = nullptr;
  MonitorElement *bend_of_tp = nullptr;
  MonitorElement *TP_pT = nullptr;
  MonitorElement *TP_pT_bendres_g1p5 = nullptr;
  MonitorElement *TP_pT_bendres_0_to_1p5 = nullptr;
  MonitorElement *TP_eta_bendres_g1p5 = nullptr;
  MonitorElement *TP_eta_bendres_0_to_1p5 = nullptr;
  MonitorElement *TP_dxy_bendres_g1p5 = nullptr;
  MonitorElement *TP_dxy_bendres_0_to_1p5 = nullptr;

  // 1D intermediate resolution plots (pT and eta)
  MonitorElement *res_eta = nullptr;    // for all eta and pT
  MonitorElement *res_pt = nullptr;     // for all eta and pT
  MonitorElement *res_ptRel = nullptr;  // for all eta and pT (delta(pT)/pT)
  MonitorElement *respt_eta0to0p7_pt2to3 = nullptr;
  MonitorElement *respt_eta0p7to1_pt2to3 = nullptr;
  MonitorElement *respt_eta1to1p2_pt2to3 = nullptr;
  MonitorElement *respt_eta1p2to1p6_pt2to3 = nullptr;
  MonitorElement *respt_eta1p6to2_pt2to3 = nullptr;
  MonitorElement *respt_eta2to2p4_pt2to3 = nullptr;
  MonitorElement *respt_eta0to0p7_pt3to8 = nullptr;
  MonitorElement *respt_eta0p7to1_pt3to8 = nullptr;
  MonitorElement *respt_eta1to1p2_pt3to8 = nullptr;
  MonitorElement *respt_eta1p2to1p6_pt3to8 = nullptr;
  MonitorElement *respt_eta1p6to2_pt3to8 = nullptr;
  MonitorElement *respt_eta2to2p4_pt3to8 = nullptr;
  MonitorElement *respt_eta0to0p7_pt8toInf = nullptr;
  MonitorElement *respt_eta0p7to1_pt8toInf = nullptr;
  MonitorElement *respt_eta1to1p2_pt8toInf = nullptr;
  MonitorElement *respt_eta1p2to1p6_pt8toInf = nullptr;
  MonitorElement *respt_eta1p6to2_pt8toInf = nullptr;
  MonitorElement *respt_eta2to2p4_pt8toInf = nullptr;
  MonitorElement *reseta_eta0to0p7 = nullptr;
  MonitorElement *reseta_eta0p7to1 = nullptr;
  MonitorElement *reseta_eta1to1p2 = nullptr;
  MonitorElement *reseta_eta1p2to1p6 = nullptr;
  MonitorElement *reseta_eta1p6to2 = nullptr;
  MonitorElement *reseta_eta2to2p4 = nullptr;
  MonitorElement *resphi_eta0to0p7 = nullptr;
  MonitorElement *resphi_eta0p7to1 = nullptr;
  MonitorElement *resphi_eta1to1p2 = nullptr;
  MonitorElement *resphi_eta1p2to1p6 = nullptr;
  MonitorElement *resphi_eta1p6to2 = nullptr;
  MonitorElement *resphi_eta2to2p4 = nullptr;
  MonitorElement *resVtxZ_eta0to0p7 = nullptr;
  MonitorElement *resVtxZ_eta0p7to1 = nullptr;
  MonitorElement *resVtxZ_eta1to1p2 = nullptr;
  MonitorElement *resVtxZ_eta1p2to1p6 = nullptr;
  MonitorElement *resVtxZ_eta1p6to2 = nullptr;
  MonitorElement *resVtxZ_eta2to2p4 = nullptr;

  // For d0
  MonitorElement *resd0_eta0to0p7 = nullptr;
  MonitorElement *resd0_eta0p7to1 = nullptr;
  MonitorElement *resd0_eta1to1p2 = nullptr;
  MonitorElement *resd0_eta1p2to1p6 = nullptr;
  MonitorElement *resd0_eta1p6to2 = nullptr;
  MonitorElement *resd0_eta2to2p4 = nullptr;

private:
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> m_topoToken;
  edm::ParameterSet conf_;
  // Declares a token, ttStubToken, that will be used to retrieve a collection of TTStub objects associated with Phase2TrackerDigi data
  edm::EDGetTokenT<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > > ttStubToken_;
  edm::EDGetTokenT<std::vector<TrackingParticle>> trackingParticleToken_;
  edm::EDGetTokenT<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>>
      ttClusterMCTruthToken_;  // MC truth association map for clusters
  edm::EDGetTokenT<TTStubAssociationMap<Ref_Phase2TrackerDigi_>>
      ttStubMCTruthToken_;  // MC truth association map for stubs
  edm::EDGetTokenT<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>>
      ttTrackMCTruthToken_;  // MC truth association map for tracks
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> getTokenTrackerGeom_;
  const edm::ESInputTag magneticFieldInputTag_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
  //trklet::Settings settings_;
  struct MatchedClusterInfo {
    GlobalPoint coordsB;
    GlobalPoint coordsC;
    unsigned int widthB;
    unsigned int widthC;
  };
  int clustersWithSingleStub = 0;
  int clustersWithMultipleStubs = 0;


  int L1Tk_minNStub;
  double L1Tk_maxChi2dof;
  int TP_minNStub;
  int TP_minNLayersStub;
  double TP_minPt;
  double TP_maxEta;
  double TP_maxVtxZ;
  std::string topFolderName_;
};
#endif
