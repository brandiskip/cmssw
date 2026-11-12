# ==============================================================================
# Phase-2 Outer Tracker DQM Efficiency Client Configuration
#
# This configuration defines DQMGenericClient modules to calculate final
# tracking and stub efficiencies for the Phase-2 Outer Tracker. It computes
# ratios for Nominal L1 tracks and Extended L1 tracks (prompt and displaced tracking particles),
# along with Barrel and Endcap stubs, reading from EfficiencyIngredients and
# outputting to FinalEfficiency folders.
#
# Author: Brandi Skipworth
# ==============================================================================

import FWCore.ParameterSet.Config as cms

# ============================================================================
# TRACK EFFICIENCY CLIENT
# ============================================================================
_topFolder = "TrackerPhase2OTL1TrackV"

phase2OTEffClient = cms.EDProducer(
    "DQMGenericClient",
    subDirs=cms.untracked.vstring(
        f"{_topFolder}/Nominal_L1TF",
        f"{_topFolder}/Extended_L1TF/PromptParticles",
        f"{_topFolder}/Extended_L1TF/DisplacedParticles",
    ),
    # Output: FinalEfficiency/Name
    # Input:  FinalEfficiency/Ingredients/Name
    efficiency=cms.vstring(
        # Nominal
        "FinalEfficiency/EtaEfficiency        'Nominal L1 tracks #eta efficiency;tracking particle #eta;Efficiency'      FinalEfficiency/Ingredients/match_tp_eta       FinalEfficiency/Ingredients/tp_eta",
        "FinalEfficiency/PtEfficiency         'Nominal L1 tracks p_{T} efficiency;p_{T} [GeV];Efficiency'                FinalEfficiency/Ingredients/match_tp_pt        FinalEfficiency/Ingredients/tp_pt",
        "FinalEfficiency/PtEfficiencyZoom     'Nominal L1 tracks p_{T} efficiency (zoom);p_{T} [GeV];Efficiency'         FinalEfficiency/Ingredients/match_tp_pt_zoom   FinalEfficiency/Ingredients/tp_pt_zoom",
        "FinalEfficiency/d0Efficiency         'Nominal L1 tracks d_{0} efficiency;d_{0} [cm];Efficiency'                 FinalEfficiency/Ingredients/match_tp_d0        FinalEfficiency/Ingredients/tp_d0",
        "FinalEfficiency/LxyEfficiency        'Nominal L1 tracks L_{xy} efficiency;L_{xy} [cm];Efficiency'               FinalEfficiency/Ingredients/match_tp_Lxy       FinalEfficiency/Ingredients/tp_Lxy",
        "FinalEfficiency/z0Efficiency         'Nominal L1 tracks z_{0} efficiency;z_{0} [cm];Efficiency'                 FinalEfficiency/Ingredients/match_tp_z0        FinalEfficiency/Ingredients/tp_z0",

        # Prompt particles
        "FinalEfficiency/EtaEfficiency        'Prompt particles #eta efficiency;tracking particle #eta;Efficiency'       FinalEfficiency/Ingredients/match_prompt_tp_eta       FinalEfficiency/Ingredients/tp_eta_for_prompt",
        "FinalEfficiency/PtEfficiency         'Prompt particles p_{T} efficiency;p_{T} [GeV];Efficiency'                 FinalEfficiency/Ingredients/match_prompt_tp_pt        FinalEfficiency/Ingredients/tp_pt_for_prompt",
        "FinalEfficiency/PtEfficiencyZoom     'Prompt particles p_{T} efficiency (zoom);p_{T} [GeV];Efficiency'          FinalEfficiency/Ingredients/match_prompt_tp_pt_zoom   FinalEfficiency/Ingredients/tp_pt_zoom_for_prompt",
        "FinalEfficiency/d0Efficiency         'Prompt particles d_{0} efficiency;d_{0} [cm];Efficiency'                  FinalEfficiency/Ingredients/match_prompt_tp_d0        FinalEfficiency/Ingredients/tp_d0_for_prompt",
        "FinalEfficiency/LxyEfficiency        'Prompt particles L_{xy} efficiency;L_{xy} [cm];Efficiency'                FinalEfficiency/Ingredients/match_prompt_tp_Lxy       FinalEfficiency/Ingredients/tp_Lxy_for_prompt",
        "FinalEfficiency/z0Efficiency         'Prompt particles z_{0} efficiency;z_{0} [cm];Efficiency'                  FinalEfficiency/Ingredients/match_prompt_tp_z0        FinalEfficiency/Ingredients/tp_z0_for_prompt",

        # Displaced particles
        "FinalEfficiency/EtaEfficiency        'Displaced particles #eta efficiency;tracking particle #eta;Efficiency'    FinalEfficiency/Ingredients/match_displaced_tp_eta       FinalEfficiency/Ingredients/tp_eta_for_dis",
        "FinalEfficiency/PtEfficiency         'Displaced particles p_{T} efficiency;p_{T} [GeV];Efficiency'              FinalEfficiency/Ingredients/match_displaced_tp_pt        FinalEfficiency/Ingredients/tp_pt_for_dis",
        "FinalEfficiency/PtEfficiencyZoom     'Displaced particles p_{T} efficiency (zoom);p_{T} [GeV];Efficiency'       FinalEfficiency/Ingredients/match_displaced_tp_pt_zoom   FinalEfficiency/Ingredients/tp_pt_zoom_for_dis",
        "FinalEfficiency/d0Efficiency         'Displaced particles d_{0} efficiency;d_{0} [cm];Efficiency'               FinalEfficiency/Ingredients/match_displaced_tp_d0        FinalEfficiency/Ingredients/tp_d0_for_dis",
        "FinalEfficiency/LxyEfficiency        'Displaced particles L_{xy} efficiency;L_{xy} [cm];Efficiency'             FinalEfficiency/Ingredients/match_displaced_tp_Lxy       FinalEfficiency/Ingredients/tp_Lxy_for_dis",
        "FinalEfficiency/z0Efficiency         'Displaced particles z_{0} efficiency;z_{0} [cm];Efficiency'               FinalEfficiency/Ingredients/match_displaced_tp_z0        FinalEfficiency/Ingredients/tp_z0_for_dis",
   ),
    resolution=cms.vstring(),
    efficiencyProfile=cms.untracked.vstring(),
    verbose=cms.untracked.uint32(0)
)

# ============================================================================
# STUB EFFICIENCY CLIENT
# ============================================================================
_topFolderStubs = "TrackerPhase2OTStubV"

phase2OTStubEffClient = cms.EDProducer(
    "DQMGenericClient",
    subDirs=cms.untracked.vstring(
        f"{_topFolderStubs}",
    ),
    # Path logic: FinalEfficiency/Name vs EfficiencyIngredients/Name
    efficiency=cms.vstring(
        "FinalEfficiency/StubEfficiencyBarrel       'Stub Efficiency Barrel;tracking particle p_{T} [GeV];Efficiency'        EfficiencyIngredients/gen_clusters_if_stub_barrel         EfficiencyIngredients/gen_clusters_barrel",
        "FinalEfficiency/StubEfficiencyZoomBarrel   'Stub Efficiency Zoom Barrel;tracking particle p_{T} [GeV];Efficiency'   EfficiencyIngredients/gen_clusters_if_stub_zoom_barrel    EfficiencyIngredients/gen_clusters_zoom_barrel",
        "FinalEfficiency/StubEfficiencyEndcaps      'Stub Efficiency Endcaps;tracking particle p_{T} [GeV];Efficiency'       EfficiencyIngredients/gen_clusters_if_stub_endcaps        EfficiencyIngredients/gen_clusters_endcaps",
        "FinalEfficiency/StubEfficiencyZoomEndcaps  'Stub Efficiency Zoom Endcaps;tracking particle p_{T} [GeV];Efficiency'  EfficiencyIngredients/gen_clusters_if_stub_zoom_endcaps   EfficiencyIngredients/gen_clusters_zoom_endcaps",
    ),
    resolution=cms.vstring(),
    efficiencyProfile=cms.untracked.vstring(),
    verbose=cms.untracked.uint32(0)
)

# Sequence
phase2OTEffClientSeq = cms.Sequence(
    phase2OTEffClient +
    phase2OTStubEffClient
)