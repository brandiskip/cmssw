// Package:    Validation/SiTrackerPhase2V
// Class:      Phase2OTHarvestTracks

/**
 * This class is part of the Phase 2 Tracker validation framework and performs
 * the harvesting step for L1 tracks validation. It processes histograms
 * created during the earlier validation steps to calculate resolutions for
 * tracking performance. Specifically it extracts standard deviations from
 * resolution ingredients for Nominal, Extended PromptParticles, and Extended DisplacedParticles
 * L1 tracks.
 * 
 * It also forms the track-level rate plots from the ingredients booked by the
 * analyzer: the fake fraction (trk_fake_vs_pt / trk_all_vs_pt) and duplicate
 * fraction (trk_duplicate_vs_pt / trk_all_vs_pt) with binomial errors, and the
 * per-event track rate (trk_all_vs_pt / nEvents), for both the nominal and
 * extended L1 track collections. Definitions follow
 * L1Trigger/TrackFindingTracklet/test/L1TrackNtuplePlot.C.
 *
 * Usage:
 * To generate histograms from this code, run the test configuration files
 * provided in the DQM/SiTrackerPhase2/test directory. The generated histograms
 * can then be analyzed or visualized.
 */

// Updated by: Brandi Skipworth, 2026

#include <algorithm>
#include <vector>
#include <string>

#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Validation/SiTrackerPhase2V/interface/TrackerPhase2HistUtil.h"

class Phase2OTHarvestTracks : public DQMEDHarvester {
public:
  explicit Phase2OTHarvestTracks(const edm::ParameterSet &);
  ~Phase2OTHarvestTracks() override;
  void dqmEndJob(DQMStore::IBooker &ibooker, DQMStore::IGetter &igetter) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  // ----------member data ---------------------------
  std::string topFolderName_;
};

Phase2OTHarvestTracks::Phase2OTHarvestTracks(const edm::ParameterSet &iConfig)
    : topFolderName_(iConfig.getParameter<std::string>("TopFolderName")) {}

Phase2OTHarvestTracks::~Phase2OTHarvestTracks() {}

// ------------ method called once each job just after ending the event loop
// ------------
void Phase2OTHarvestTracks::dqmEndJob(DQMStore::IBooker &ibooker, DQMStore::IGetter &igetter) {
  using namespace edm;

  float eta_bins[] = {0.0, 0.7, 1.0, 1.2, 1.6, 2.0, 2.4};
  int eta_binnum = 6;
  std::string eta_ranges[6] = {"eta0to0p7", "eta0p7to1", "eta1to1p2", "eta1p2to1p6", "eta1p6to2", "eta2to2p4"};

  // nominal collection
  std::vector<MonitorElement *> respt_pt2to3(6, nullptr), respt_pt3to8(6, nullptr), respt_pt8toInf(6, nullptr);
  std::vector<MonitorElement *> mereseta_vect(6, nullptr), meresphi_vect(6, nullptr), meresz0_vect(6, nullptr);

  std::vector<MonitorElement *> prompt_respt_pt2to3(6, nullptr), prompt_respt_pt3to8(6, nullptr),
      prompt_respt_pt8toInf(6, nullptr);
  std::vector<MonitorElement *> prompt_mereseta_vect(6, nullptr), prompt_meresphi_vect(6, nullptr),
      prompt_meresz0_vect(6, nullptr), prompt_meresd0_vect(6, nullptr);

  std::vector<MonitorElement *> displaced_respt_pt2to3(6, nullptr), displaced_respt_pt3to8(6, nullptr),
      displaced_respt_pt8toInf(6, nullptr);
  std::vector<MonitorElement *> displaced_mereseta_vect(6, nullptr), displaced_meresphi_vect(6, nullptr),
      displaced_meresz0_vect(6, nullptr), displaced_meresd0_vect(6, nullptr);

  for (int i = 0; i < 6; i++) {
    std::string resIng = topFolderName_ + "/Nominal_L1TF/FinalResolution/Ingredients/";
    respt_pt2to3[i] = igetter.get(resIng + "respt_" + eta_ranges[i] + "_pt2to3");
    respt_pt3to8[i] = igetter.get(resIng + "respt_" + eta_ranges[i] + "_pt3to8");
    respt_pt8toInf[i] = igetter.get(resIng + "respt_" + eta_ranges[i] + "_pt8toInf");
    mereseta_vect[i] = igetter.get(resIng + "reseta_" + eta_ranges[i]);
    meresphi_vect[i] = igetter.get(resIng + "resphi_" + eta_ranges[i]);
    meresz0_vect[i] = igetter.get(resIng + "resz0_" + eta_ranges[i]);

    std::string promptIng = topFolderName_ + "/Extended_L1TF/PromptParticles/FinalResolution/Ingredients/";
    prompt_respt_pt2to3[i] = igetter.get(promptIng + "respt_prompt_" + eta_ranges[i] + "_pt2to3");
    prompt_respt_pt3to8[i] = igetter.get(promptIng + "respt_prompt_" + eta_ranges[i] + "_pt3to8");
    prompt_respt_pt8toInf[i] = igetter.get(promptIng + "respt_prompt_" + eta_ranges[i] + "_pt8toInf");
    prompt_mereseta_vect[i] = igetter.get(promptIng + "reseta_prompt_" + eta_ranges[i]);
    prompt_meresphi_vect[i] = igetter.get(promptIng + "resphi_prompt_" + eta_ranges[i]);
    prompt_meresz0_vect[i] = igetter.get(promptIng + "resz0_prompt_" + eta_ranges[i]);
    prompt_meresd0_vect[i] = igetter.get(promptIng + "resd0_prompt_" + eta_ranges[i]);

    std::string dispIng = topFolderName_ + "/Extended_L1TF/DisplacedParticles/FinalResolution/Ingredients/";
    displaced_respt_pt2to3[i] = igetter.get(dispIng + "respt_displaced_" + eta_ranges[i] + "_pt2to3");
    displaced_respt_pt3to8[i] = igetter.get(dispIng + "respt_displaced_" + eta_ranges[i] + "_pt3to8");
    displaced_respt_pt8toInf[i] = igetter.get(dispIng + "respt_displaced_" + eta_ranges[i] + "_pt8toInf");
    displaced_mereseta_vect[i] = igetter.get(dispIng + "reseta_displaced_" + eta_ranges[i]);
    displaced_meresphi_vect[i] = igetter.get(dispIng + "resphi_displaced_" + eta_ranges[i]);
    displaced_meresz0_vect[i] = igetter.get(dispIng + "resz0_displaced_" + eta_ranges[i]);
    displaced_meresd0_vect[i] = igetter.get(dispIng + "resd0_displaced_" + eta_ranges[i]);
  }

  // Final Resolutions - Nominal
  ibooker.setCurrentFolder(topFolderName_ + "/Nominal_L1TF/FinalResolution");
  phase2tkutil::fillResolutionFromVec(
      respt_pt2to3,
      ibooker.book1D("pTResVsEta_2-3",
                     "Nominal L1 tracks p_{T} resolution vs |#eta| (2-3 GeV);|#eta|;#sigma(#Delta p_{T}/p_{T})",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      respt_pt3to8,
      ibooker.book1D("pTResVsEta_3-8",
                     "Nominal L1 tracks p_{T} resolution vs |#eta| (3-8 GeV);|#eta|;#sigma(#Delta p_{T}/p_{T})",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      respt_pt8toInf,
      ibooker.book1D("pTResVsEta_8-inf",
                     "Nominal L1 tracks p_{T} resolution vs |#eta| (> 8 GeV);|#eta|;#sigma(#Delta p_{T}/p_{T})",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      mereseta_vect,
      ibooker.book1D("EtaResolution",
                     "Nominal L1 tracks #eta resolution vs |#eta|;|#eta|;#sigma(#Delta#eta)",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      meresphi_vect,
      ibooker.book1D("PhiResolution",
                     "Nominal L1 tracks #phi resolution vs |#eta|;|#eta|;#sigma(#Delta#phi)",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      meresz0_vect,
      ibooker.book1D("z0Resolution",
                     "Nominal L1 tracks z_{0} resolution vs |#eta|;|#eta|;#sigma(#Deltaz_{0}) [cm]",
                     eta_binnum,
                     eta_bins),
      "");

  // Final Resolutions - PromptParticles
  ibooker.setCurrentFolder(topFolderName_ + "/Extended_L1TF/PromptParticles/FinalResolution");
  phase2tkutil::fillResolutionFromVec(
      prompt_respt_pt2to3,
      ibooker.book1D("pTResVsEta_2-3_prompt",
                     "Prompt particles p_{T} resolution vs |#eta| (2-3 GeV);|#eta|;#sigma(#Delta p_{T}/p_{T})",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      prompt_respt_pt3to8,
      ibooker.book1D("pTResVsEta_3-8_prompt",
                     "Prompt particles p_{T} resolution vs |#eta| (3-8 GeV);|#eta|;#sigma(#Delta p_{T}/p_{T})",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      prompt_respt_pt8toInf,
      ibooker.book1D("pTResVsEta_8-inf_prompt",
                     "Prompt particles p_{T} resolution vs |#eta| (> 8 GeV);|#eta|;#sigma(#Delta p_{T}/p_{T})",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      prompt_mereseta_vect,
      ibooker.book1D("EtaResolution_prompt",
                     "Prompt particles #eta resolution vs |#eta|;|#eta|;#sigma(#Delta#eta)",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      prompt_meresphi_vect,
      ibooker.book1D("PhiResolution_prompt",
                     "Prompt particles #phi resolution vs |#eta|;|#eta|;#sigma(#Delta#phi)",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      prompt_meresz0_vect,
      ibooker.book1D("z0Resolution_prompt",
                     "Prompt particles z_{0} resolution vs |#eta|;|#eta|;#sigma(#Deltaz_{0}) [cm]",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      prompt_meresd0_vect,
      ibooker.book1D("d0Resolution_prompt",
                     "Prompt particles d_{0} resolution vs |#eta|;|#eta|;#sigma(#Deltad_{0}) [cm]",
                     eta_binnum,
                     eta_bins),
      "");

  // Final Resolutions - DisplacedParticles
  ibooker.setCurrentFolder(topFolderName_ + "/Extended_L1TF/DisplacedParticles/FinalResolution");
  phase2tkutil::fillResolutionFromVec(
      displaced_respt_pt2to3,
      ibooker.book1D("pTResVsEta_2-3_displaced",
                     "Displaced particles p_{T} resolution vs |#eta| (2-3 GeV);|#eta|;#sigma(#Delta p_{T}/p_{T})",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      displaced_respt_pt3to8,
      ibooker.book1D("pTResVsEta_3-8_displaced",
                     "Displaced particles p_{T} resolution vs |#eta| (3-8 GeV);|#eta|;#sigma(#Delta p_{T}/p_{T})",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      displaced_respt_pt8toInf,
      ibooker.book1D("pTResVsEta_8-inf_displaced",
                     "Displaced particles p_{T} resolution vs |#eta| (> 8 GeV);|#eta|;#sigma(#Delta p_{T}/p_{T})",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      displaced_mereseta_vect,
      ibooker.book1D("EtaResolution_displaced",
                     "Displaced particles #eta resolution vs |#eta|;|#eta|;#sigma(#Delta#eta)",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      displaced_meresphi_vect,
      ibooker.book1D("PhiResolution_displaced",
                     "Displaced particles #phi resolution vs |#eta|;|#eta|;#sigma(#Delta#phi)",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      displaced_meresz0_vect,
      ibooker.book1D("z0Resolution_displaced",
                     "Displaced particles z_{0} resolution vs |#eta|;|#eta|;#sigma(#Deltaz_{0}) [cm]",
                     eta_binnum,
                     eta_bins),
      "");
  phase2tkutil::fillResolutionFromVec(
      displaced_meresd0_vect,
      ibooker.book1D("d0Resolution_displaced",
                     "Displaced particles d_{0} resolution vs |#eta|;|#eta|;#sigma(#Deltad_{0}) [cm]",
                     eta_binnum,
                     eta_bins),
      "");

  // ---------------------------------------------------------------
  // Track rate, fake fraction, duplicate fraction
  // (ingredients filled in Phase2OTValidateTracks; definitions follow
  //  L1Trigger/TrackFindingTracklet/test/L1TrackNtuplePlot.C)

  // number of processed events, from the once-per-event TP multiplicity histogram
  double nEvents = 0;
  if (MonitorElement *meNevt = igetter.get(topFolderName_ + "/trackParticles/trackParts_Num"))
    nEvents = meNevt->getTH1F()->GetEntries();

  std::string ratesNom = topFolderName_ + "/Nominal_L1TF/TrackRateIngredients/";
  ibooker.setCurrentFolder(topFolderName_ + "/Nominal_L1TF/FinalRates");
  phase2tkutil::bookFraction(ibooker,
                             igetter,
                             ratesNom + "trk_fake_vs_pt",
                             ratesNom + "trk_all_vs_pt",
                             "fakeFraction_vs_pt",
                             "Nominal L1 tracks fake fraction;track p_{T} [GeV];fake fraction");
  phase2tkutil::bookFraction(ibooker,
                             igetter,
                             ratesNom + "trk_duplicate_vs_pt",
                             ratesNom + "trk_all_vs_pt",
                             "duplicateFraction_vs_pt",
                             "Nominal L1 tracks duplicate fraction;track p_{T} [GeV];duplicate fraction");
  phase2tkutil::bookPerEventRate(ibooker,
                                 igetter,
                                 ratesNom + "trk_all_vs_pt",
                                 "trackRate_vs_pt",
                                 "Nominal L1 tracks per event;track p_{T} [GeV];tracks / event",
                                 nEvents);

  std::string ratesExt = topFolderName_ + "/Extended_L1TF/TrackRateIngredients/";
  ibooker.setCurrentFolder(topFolderName_ + "/Extended_L1TF/FinalRates");
  phase2tkutil::bookFraction(ibooker,
                             igetter,
                             ratesExt + "trk_fake_vs_pt_extended",
                             ratesExt + "trk_all_vs_pt_extended",
                             "fakeFraction_vs_pt",
                             "Extended L1 tracks fake fraction;track p_{T} [GeV];fake fraction");
  phase2tkutil::bookFraction(ibooker,
                             igetter,
                             ratesExt + "trk_duplicate_vs_pt_extended",
                             ratesExt + "trk_all_vs_pt_extended",
                             "duplicateFraction_vs_pt",
                             "Extended L1 tracks duplicate fraction;track p_{T} [GeV];duplicate fraction");
  phase2tkutil::bookPerEventRate(ibooker,
                                 igetter,
                                 ratesExt + "trk_all_vs_pt_extended",
                                 "trackRate_vs_pt",
                                 "Extended L1 tracks per event;track p_{T} [GeV];tracks / event",
                                 nEvents);
}  // end dqmEndJob

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
void Phase2OTHarvestTracks::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("TopFolderName", "TrackerPhase2OTL1TrackV");
  descriptions.add("Phase2OTHarvestTracks", desc);
}
DEFINE_FWK_MODULE(Phase2OTHarvestTracks);
