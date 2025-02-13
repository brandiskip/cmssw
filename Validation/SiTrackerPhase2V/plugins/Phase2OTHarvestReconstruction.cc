// Package:    Validation/SiTrackerPhase2V
// Class:      Phase2OTHarvestReconstruction

/**
 * This class is part of the Phase 2 Tracker validation framework and performs
 * the harvesting step for tracking particle validation. It processes histograms
 * created during the earlier validation steps to calculate efficiencies and
 * resolutions for stub reconstruction and tracking performance.
 *
 * Usage:
 * To generate histograms from this code, run the test configuration files provided
 * in the DQM/SiTrackerPhase2/test directory. The generated histograms can then be
 * analyzed or visualized.
 */

// Updated by: Brandi Skipworth, 2025

#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <algorithm>

class Phase2OTHarvestReconstruction : public DQMEDHarvester {
public:
  explicit Phase2OTHarvestReconstruction(const edm::ParameterSet &);
  ~Phase2OTHarvestReconstruction() override;
  void dqmEndJob(DQMStore::IBooker &ibooker, DQMStore::IGetter &igetter) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  // ----------member data ---------------------------
  std::string topFolderName_;
};

Phase2OTHarvestReconstruction::Phase2OTHarvestReconstruction(const edm::ParameterSet &iConfig)
    : topFolderName_(iConfig.getParameter<std::string>("TopFolderName")) {}

Phase2OTHarvestReconstruction::~Phase2OTHarvestReconstruction() {}

// ------------ method called once each job just after ending the event loop
// ------------
void Phase2OTHarvestReconstruction::dqmEndJob(DQMStore::IBooker &ibooker, DQMStore::IGetter &igetter) {
  using namespace edm;

  float eta_bins[] = {0.0, 0.7, 1.0, 1.2, 1.6, 2.0, 2.4};
  int eta_binnum = 6;

  // Find all monitor elements for histograms
  MonitorElement *meN_clus_barrel = igetter.get(topFolderName_ + "/EfficiencyIngredients/gen_clusters_if_stub_barrel");
  MonitorElement *meD_clus_barrel = igetter.get(topFolderName_ + "/EfficiencyIngredients/gen_clusters_barrel");
  MonitorElement *meN_clus_zoom_barrel =
      igetter.get(topFolderName_ + "/EfficiencyIngredients/gen_clusters_if_stub_zoom_barrel");
  MonitorElement *meD_clus_zoom_barrel =
      igetter.get(topFolderName_ + "/EfficiencyIngredients/gen_clusters_zoom_barrel");
  MonitorElement *meN_clus_endcaps =
      igetter.get(topFolderName_ + "/EfficiencyIngredients/gen_clusters_if_stub_endcaps");
  MonitorElement *meD_clus_endcaps = igetter.get(topFolderName_ + "/EfficiencyIngredients/gen_clusters_endcaps");
  MonitorElement *meN_clus_zoom_endcaps =
      igetter.get(topFolderName_ + "/EfficiencyIngredients/gen_clusters_if_stub_zoom_endcaps");
  MonitorElement *meD_clus_zoom_endcaps =
      igetter.get(topFolderName_ + "/EfficiencyIngredients/gen_clusters_zoom_endcaps");
  MonitorElement *meN_eta = igetter.get(topFolderName_ + "/EfficiencyIngredients/match_tp_eta");
  MonitorElement *meD_eta = igetter.get(topFolderName_ + "/EfficiencyIngredients/tp_eta");
  MonitorElement *meN_pt = igetter.get(topFolderName_ + "/EfficiencyIngredients/match_tp_pt");
  MonitorElement *meD_pt = igetter.get(topFolderName_ + "/EfficiencyIngredients/tp_pt");
  MonitorElement *meN_pt_zoom = igetter.get(topFolderName_ + "/EfficiencyIngredients/match_tp_pt_zoom");
  MonitorElement *meD_pt_zoom = igetter.get(topFolderName_ + "/EfficiencyIngredients/tp_pt_zoom");
  MonitorElement *meN_d0 = igetter.get(topFolderName_ + "/EfficiencyIngredients/match_tp_d0");
  MonitorElement *meD_d0 = igetter.get(topFolderName_ + "/EfficiencyIngredients/tp_d0");
  MonitorElement *meN_VtxR = igetter.get(topFolderName_ + "/EfficiencyIngredients/match_tp_VtxR");
  MonitorElement *meD_VtxR = igetter.get(topFolderName_ + "/EfficiencyIngredients/tp_VtxR");
  MonitorElement *meN_VtxZ = igetter.get(topFolderName_ + "/EfficiencyIngredients/match_tp_VtxZ");
  MonitorElement *meD_VtxZ = igetter.get(topFolderName_ + "/EfficiencyIngredients/tp_VtxZ");

  std::string eta_ranges[6] = {"eta0to0p7", "eta0p7to1", "eta1to1p2", "eta1p2to1p6", "eta1p6to2", "eta2to2p4"};

  std::vector<MonitorElement*> respt_pt2to3 = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  std::vector<MonitorElement*> respt_pt3to8 = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  std::vector<MonitorElement*> respt_pt8toInf = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  std::vector<MonitorElement*> mereseta_vect = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  std::vector<MonitorElement*> meresphi_vect = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  std::vector<MonitorElement*> meresVtxZ_vect = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  std::vector<MonitorElement*> meresd0_vect = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

  for (int i=0; i<5; i++){
    respt_pt2to3[i] = igetter.get(topFolderName_ + "/ResolutionIngredients/respt_" + eta_ranges[i] + "_pt2to3");
    respt_pt3to8[i] = igetter.get(topFolderName_ + "/ResolutionIngredients/respt_" + eta_ranges[i] + "_pt3to8");
    respt_pt8toInf[i] = igetter.get(topFolderName_ + "/ResolutionIngredients/respt_" + eta_ranges[i] + "_pt8toInf");
    mereseta_vect[i] = igetter.get(topFolderName_ + "/ResolutionIngredients/reseta_" + eta_ranges[i]);
    meresphi_vect[i] = igetter.get(topFolderName_ + "/ResolutionIngredients/resphi_" + eta_ranges[i]);
    meresVtxZ_vect[i] = igetter.get(topFolderName_ + "/ResolutionIngredients/resVtxZ_" + eta_ranges[i]);
    meresd0_vect[i] = igetter.get(topFolderName_ + "/ResolutionIngredients/resd0_" + eta_ranges[i]);
  }

  if (meN_clus_barrel && meD_clus_barrel) {
    // Get the numerator and denominator histograms
    TH1F *numerator = meN_clus_barrel->getTH1F();
    TH1F *denominator = meD_clus_barrel->getTH1F();
    numerator->Sumw2();
    denominator->Sumw2();

    // Set the current directory
    igetter.setCurrentFolder(topFolderName_ + "/FinalEfficiency");

    // Book the new histogram to contain the results
    MonitorElement *me_effic_clus_barrel = ibooker.book1D("StubEfficiencyBarrel",
                                                          "Stub Efficiency Barrel",
                                                          numerator->GetNbinsX(),
                                                          numerator->GetXaxis()->GetXmin(),
                                                          numerator->GetXaxis()->GetXmax());

    // Calculate the efficiency
    me_effic_clus_barrel->getTH1F()->Divide(numerator, denominator, 1., 1., "B");
    me_effic_clus_barrel->setAxisTitle("tracking particle pT [GeV]");
    me_effic_clus_barrel->getTH1F()->GetYaxis()->SetTitle("Efficiency");
    me_effic_clus_barrel->getTH1F()->SetMaximum(1.1);
    me_effic_clus_barrel->getTH1F()->SetMinimum(0.0);
    me_effic_clus_barrel->getTH1F()->SetStats(false);
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for stub efficiency barrel cannot be found!\n";
  }

  if (meN_clus_zoom_barrel && meD_clus_zoom_barrel) {
    // Get the numerator and denominator histograms
    TH1F *numerator_zoom = meN_clus_zoom_barrel->getTH1F();
    TH1F *denominator_zoom = meD_clus_zoom_barrel->getTH1F();
    numerator_zoom->Sumw2();
    denominator_zoom->Sumw2();

    // Set the current directory
    igetter.setCurrentFolder(topFolderName_ + "/FinalEfficiency");

    // Book the new histogram to contain the results
    MonitorElement *me_effic_clus_zoom_barrel = ibooker.book1D("StubEfficiencyZoomBarrel",
                                                               "Stub Efficiency Zoom Barrel",
                                                               numerator_zoom->GetNbinsX(),
                                                               numerator_zoom->GetXaxis()->GetXmin(),
                                                               numerator_zoom->GetXaxis()->GetXmax());

    // Calculate the efficiency
    me_effic_clus_zoom_barrel->getTH1F()->Divide(numerator_zoom, denominator_zoom, 1., 1., "B");
    me_effic_clus_zoom_barrel->setAxisTitle("tracking particle pT [GeV]");
    me_effic_clus_zoom_barrel->getTH1F()->GetYaxis()->SetTitle("Efficiency");
    me_effic_clus_zoom_barrel->getTH1F()->SetMaximum(1.1);
    me_effic_clus_zoom_barrel->getTH1F()->SetMinimum(0.0);
    me_effic_clus_zoom_barrel->getTH1F()->SetStats(false);
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for stub zoom barrel efficiency cannot be found!\n";
  }

  if (meN_clus_endcaps && meD_clus_endcaps) {
    // Get the numerator and denominator histograms
    TH1F *numerator = meN_clus_endcaps->getTH1F();
    TH1F *denominator = meD_clus_endcaps->getTH1F();
    numerator->Sumw2();
    denominator->Sumw2();

    // Set the current directory
    igetter.setCurrentFolder(topFolderName_ + "/FinalEfficiency");

    // Book the new histogram to contain the results
    MonitorElement *me_effic_clus_endcaps = ibooker.book1D("StubEfficiencyEndcaps",
                                                           "Stub Efficiency Endcaps",
                                                           numerator->GetNbinsX(),
                                                           numerator->GetXaxis()->GetXmin(),
                                                           numerator->GetXaxis()->GetXmax());

    // Calculate the efficiency
    me_effic_clus_endcaps->getTH1F()->Divide(numerator, denominator, 1., 1., "B");
    me_effic_clus_endcaps->setAxisTitle("tracking particle pT [GeV]");
    me_effic_clus_endcaps->getTH1F()->GetYaxis()->SetTitle("Efficiency");
    me_effic_clus_endcaps->getTH1F()->SetMaximum(1.1);
    me_effic_clus_endcaps->getTH1F()->SetMinimum(0.0);
    me_effic_clus_endcaps->getTH1F()->SetStats(false);
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for stub efficiency endcaps cannot be found!\n";
  }

  if (meN_clus_zoom_endcaps && meD_clus_zoom_endcaps) {
    // Get the numerator and denominator histograms
    TH1F *numerator_zoom = meN_clus_zoom_endcaps->getTH1F();
    TH1F *denominator_zoom = meD_clus_zoom_endcaps->getTH1F();
    numerator_zoom->Sumw2();
    denominator_zoom->Sumw2();

    // Set the current directory
    igetter.setCurrentFolder(topFolderName_ + "/FinalEfficiency");

    // Book the new histogram to contain the results
    MonitorElement *me_effic_clus_zoom_endcaps = ibooker.book1D("StubEfficiencyZoomEndcaps",
                                                                "Stub Efficiency Zoom Endcaps",
                                                                numerator_zoom->GetNbinsX(),
                                                                numerator_zoom->GetXaxis()->GetXmin(),
                                                                numerator_zoom->GetXaxis()->GetXmax());

    // Calculate the efficiency
    me_effic_clus_zoom_endcaps->getTH1F()->Divide(numerator_zoom, denominator_zoom, 1., 1., "B");
    me_effic_clus_zoom_endcaps->setAxisTitle("tracking particle pT [GeV]");
    me_effic_clus_zoom_endcaps->getTH1F()->GetYaxis()->SetTitle("Efficiency");
    me_effic_clus_zoom_endcaps->getTH1F()->SetMaximum(1.1);
    me_effic_clus_zoom_endcaps->getTH1F()->SetMinimum(0.0);
    me_effic_clus_zoom_endcaps->getTH1F()->SetStats(false);
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for stub zoom endcaps efficiency cannot be found!\n";
  }

  if (meN_eta && meD_eta) {
    // Get the numerator and denominator histograms
    TH1F *numerator = meN_eta->getTH1F();
    TH1F *denominator = meD_eta->getTH1F();
    numerator->Sumw2();
    denominator->Sumw2();

    // Set the current directory
    igetter.setCurrentFolder(topFolderName_ + "/FinalEfficiency");

    // Book the new histogram to contain the results
    MonitorElement *me_effic_eta = ibooker.book1D("EtaEfficiency",
                                                  "#eta efficiency",
                                                  numerator->GetNbinsX(),
                                                  numerator->GetXaxis()->GetXmin(),
                                                  numerator->GetXaxis()->GetXmax());

    // Calculate the efficiency
    me_effic_eta->getTH1F()->Divide(numerator, denominator, 1., 1., "B");
    me_effic_eta->setAxisTitle("tracking particle #eta");
    me_effic_eta->getTH1F()->GetYaxis()->SetTitle("Efficiency");
    me_effic_eta->getTH1F()->SetMaximum(1.0);
    me_effic_eta->getTH1F()->SetMinimum(0.0);
    me_effic_eta->getTH1F()->SetStats(false);
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for eta efficiency cannot be found!\n";
  }

  if (meN_pt && meD_pt) {
    // Get the numerator and denominator histograms
    TH1F *numerator2 = meN_pt->getTH1F();
    numerator2->Sumw2();
    TH1F *denominator2 = meD_pt->getTH1F();
    denominator2->Sumw2();

    // Set the current directory
    igetter.setCurrentFolder(topFolderName_ + "/FinalEfficiency");

    // Book the new histogram to contain the results
    MonitorElement *me_effic_pt = ibooker.book1D("PtEfficiency",
                                                 "p_{T} efficiency",
                                                 numerator2->GetNbinsX(),
                                                 numerator2->GetXaxis()->GetXmin(),
                                                 numerator2->GetXaxis()->GetXmax());

    // Calculate the efficiency
    me_effic_pt->getTH1F()->Divide(numerator2, denominator2, 1., 1., "B");
    me_effic_pt->setAxisTitle("Tracking particle p_{T} [GeV]");
    me_effic_pt->getTH1F()->GetYaxis()->SetTitle("Efficiency");
    me_effic_pt->getTH1F()->SetMaximum(1.0);
    me_effic_pt->getTH1F()->SetMinimum(0.0);
    me_effic_pt->getTH1F()->SetStats(false);
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for pT efficiency cannot be found!\n";
  }

  if (meN_pt_zoom && meD_pt_zoom) {
    // Get the numerator and denominator histograms
    TH1F *numerator2_zoom = meN_pt_zoom->getTH1F();
    numerator2_zoom->Sumw2();
    TH1F *denominator2_zoom = meD_pt_zoom->getTH1F();
    denominator2_zoom->Sumw2();

    // Set the current directory
    igetter.setCurrentFolder(topFolderName_ + "/FinalEfficiency");

    // Book the new histogram to contain the results
    MonitorElement *me_effic_pt_zoom = ibooker.book1D("PtEfficiency_zoom",
                                                      "p_{T} efficiency",
                                                      numerator2_zoom->GetNbinsX(),
                                                      numerator2_zoom->GetXaxis()->GetXmin(),
                                                      numerator2_zoom->GetXaxis()->GetXmax());

    // Calculate the efficiency
    me_effic_pt_zoom->getTH1F()->Divide(numerator2_zoom, denominator2_zoom, 1., 1., "B");
    me_effic_pt_zoom->setAxisTitle("Tracking particle p_{T} [GeV]");
    me_effic_pt_zoom->getTH1F()->GetYaxis()->SetTitle("Efficiency");
    me_effic_pt_zoom->getTH1F()->SetMaximum(1.0);
    me_effic_pt_zoom->getTH1F()->SetMinimum(0.0);
    me_effic_pt_zoom->getTH1F()->SetStats(false);
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for zoom pT efficiency cannot be found!\n";
  }

  if (meN_d0 && meD_d0) {
    // Get the numerator and denominator histograms
    TH1F *numerator5 = meN_d0->getTH1F();
    numerator5->Sumw2();
    TH1F *denominator5 = meD_d0->getTH1F();
    denominator5->Sumw2();

    // Set the current directory
    igetter.setCurrentFolder(topFolderName_ + "/FinalEfficiency");

    // Book the new histogram to contain the results
    MonitorElement *me_effic_d0 = ibooker.book1D("d0Efficiency",
                                                 "d_{0} efficiency",
                                                 numerator5->GetNbinsX(),
                                                 numerator5->GetXaxis()->GetXmin(),
                                                 numerator5->GetXaxis()->GetXmax());

    // Calculate the efficiency
    me_effic_d0->getTH1F()->Divide(numerator5, denominator5, 1., 1., "B");
    me_effic_d0->setAxisTitle("Tracking particle d_{0} [cm]");
    me_effic_d0->getTH1F()->GetYaxis()->SetTitle("Efficiency");
    me_effic_d0->getTH1F()->SetMaximum(1.0);
    me_effic_d0->getTH1F()->SetMinimum(0.0);
    me_effic_d0->getTH1F()->SetStats(false);
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for d0 efficiency cannot be found!\n";
  }

  if (meN_VtxR && meD_VtxR) {
    // Get the numerator and denominator histograms
    TH1F *numerator6 = meN_VtxR->getTH1F();
    numerator6->Sumw2();
    TH1F *denominator6 = meD_VtxR->getTH1F();
    denominator6->Sumw2();

    // Set the current directory
    igetter.setCurrentFolder(topFolderName_ + "/FinalEfficiency");

    // Book the new histogram to contain the results
    MonitorElement *me_effic_VtxR = ibooker.book1D("VtxREfficiency",
                                                   "Vtx R efficiency",
                                                   numerator6->GetNbinsX(),
                                                   numerator6->GetXaxis()->GetXmin(),
                                                   numerator6->GetXaxis()->GetXmax());

    // Calculate the efficiency
    me_effic_VtxR->getTH1F()->Divide(numerator6, denominator6, 1., 1., "B");
    me_effic_VtxR->setAxisTitle("Tracking particle VtxR [cm]");
    me_effic_VtxR->getTH1F()->GetYaxis()->SetTitle("Efficiency");
    me_effic_VtxR->getTH1F()->SetMaximum(1.0);
    me_effic_VtxR->getTH1F()->SetMinimum(0.0);
    me_effic_VtxR->getTH1F()->SetStats(false);
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for VtxR efficiency cannot be found!\n";
  }

  if (meN_VtxZ && meD_VtxZ) {
    // Get the numerator and denominator histograms
    TH1F *numerator7 = meN_VtxZ->getTH1F();
    numerator7->Sumw2();
    TH1F *denominator7 = meD_VtxZ->getTH1F();
    denominator7->Sumw2();

    // Set the current directory
    igetter.setCurrentFolder(topFolderName_ + "/FinalEfficiency");

    // Book the new histogram to contain the results
    MonitorElement *me_effic_VtxZ = ibooker.book1D("VtxZEfficiency",
                                                   "Vtx Z efficiency",
                                                   numerator7->GetNbinsX(),
                                                   numerator7->GetXaxis()->GetXmin(),
                                                   numerator7->GetXaxis()->GetXmax());

    // Calculate the efficiency
    me_effic_VtxZ->getTH1F()->Divide(numerator7, denominator7, 1., 1., "B");
    me_effic_VtxZ->setAxisTitle("Tracking particle VtxZ [cm]");
    me_effic_VtxZ->getTH1F()->GetYaxis()->SetTitle("Efficiency");
    me_effic_VtxZ->getTH1F()->SetMaximum(1.0);
    me_effic_VtxZ->getTH1F()->SetMinimum(0.0);
    me_effic_VtxZ->getTH1F()->SetStats(false);
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for VtxZ efficiency cannot be found!\n";
  }

  if (std::find(respt_pt2to3.begin(), respt_pt2to3.end(), nullptr)==respt_pt2to3.end()) {
    // Set the current directoy
    igetter.setCurrentFolder(topFolderName_ + "/FinalResolution");

    // Grab the histograms
    std::vector<TH1F *> vResPt1 = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    for (int i=0; i<6; i++){
      vResPt1[i] = respt_pt2to3[i]->getTH1F();
    }

    // Book the new histogram to contain the results
    MonitorElement *me_res_pt1 =
        ibooker.book1D("pTResVsEta_2-3", "p_{T} resolution vs |#eta|, for p_{T}: 2-3 GeV", eta_binnum, eta_bins);
    TH1F *resPt1 = me_res_pt1->getTH1F();
    resPt1->GetXaxis()->SetTitle("tracking particle |#eta|");
    resPt1->GetYaxis()->SetTitle("#sigma(#Deltap_{T}/p_{T})");
    resPt1->SetMinimum(0.0);
    resPt1->SetStats(false);

    for (int i = 0; i < 6; i++) {
      resPt1->SetBinContent(i + 1, vResPt1[i]->GetStdDev());
      resPt1->SetBinError(i + 1, vResPt1[i]->GetStdDevError());
    }
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for pT resolution (2-3) cannot be found!\n";
  }

  if (std::find(respt_pt3to8.begin(), respt_pt3to8.end(), nullptr)==respt_pt3to8.end()) {
    // Set the current directoy
    igetter.setCurrentFolder(topFolderName_ + "/FinalResolution");

    // Grab the histograms
    std::vector<TH1F *> vResPt2 = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    for (int i=0; i<6; i++){
      vResPt2[i] = respt_pt3to8[i]->getTH1F();
    }

    // Book the new histogram to contain the results
    MonitorElement *me_res_pt2 =
        ibooker.book1D("pTResVsEta_3-8", "p_{T} resolution vs |#eta|, for p_{T}: 3-8 GeV", eta_binnum, eta_bins);
    TH1F *resPt2 = me_res_pt2->getTH1F();
    resPt2->GetXaxis()->SetTitle("tracking particle |#eta|");
    resPt2->GetYaxis()->SetTitle("#sigma(#Deltap_{T}/p_{T})");
    resPt2->SetMinimum(0.0);
    resPt2->SetStats(false);

    for (int i = 0; i < 6; i++) {
      resPt2->SetBinContent(i + 1, vResPt2[i]->GetStdDev());
      resPt2->SetBinError(i + 1, vResPt2[i]->GetStdDevError());
    }
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for pT resolution (3-8) cannot be found!\n";
  }

  if (std::find(respt_pt8toInf.begin(), respt_pt8toInf.end(), nullptr)==respt_pt8toInf.end()) {
    // Set the current directoy
    igetter.setCurrentFolder(topFolderName_ + "/FinalResolution");

    // Grab the histograms
    std::vector<TH1F *> vResPt3 = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    for (int i=0; i<6; i++){
      vResPt3[i] = respt_pt8toInf[i]->getTH1F();
    }

    // Book the new histogram to contain the results
    MonitorElement *me_res_pt3 =
        ibooker.book1D("pTResVsEta_8-inf", "p_{T} resolution vs |#eta|, for p_{T}: >8 GeV", eta_binnum, eta_bins);
    TH1F *resPt3 = me_res_pt3->getTH1F();
    resPt3->GetXaxis()->SetTitle("tracking particle |#eta|");
    resPt3->GetYaxis()->SetTitle("#sigma(#Deltap_{T}/p_{T})");
    resPt3->SetMinimum(0.0);
    resPt3->SetStats(false);

    for (int i = 0; i < 6; i++) {
      resPt3->SetBinContent(i + 1, vResPt3[i]->GetStdDev());
      resPt3->SetBinError(i + 1, vResPt3[i]->GetStdDevError());
    }
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for pT resolution (8-inf) cannot be found!\n";
  }

  if (std::find(mereseta_vect.begin(), mereseta_vect.end(), nullptr)==mereseta_vect.end()) {
    // Set the current directoy
    igetter.setCurrentFolder(topFolderName_ + "/FinalResolution");

    // Grab the histograms
    std::vector<TH1F *> vResEta = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    for (int i=0; i<6; i++){
      vResEta[i] = mereseta_vect[i]->getTH1F();
    }

    // Book the new histogram to contain the results
    MonitorElement *me_res_eta = ibooker.book1D("EtaResolution", "#eta resolution vs |#eta|", eta_binnum, eta_bins);
    TH1F *resEta = me_res_eta->getTH1F();
    resEta->GetXaxis()->SetTitle("tracking particle |#eta|");
    resEta->GetYaxis()->SetTitle("#sigma(#Delta#eta)");
    resEta->SetMinimum(0.0);
    resEta->SetStats(false);

    for (int i = 0; i < 6; i++) {
      resEta->SetBinContent(i + 1, vResEta[i]->GetStdDev());
      resEta->SetBinError(i + 1, vResEta[i]->GetStdDevError());
    }
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for eta resolution cannot be found!\n";
  }

  if (std::find(meresphi_vect.begin(), meresphi_vect.end(), nullptr)==meresphi_vect.end()) {
    // Set the current directoy
    igetter.setCurrentFolder(topFolderName_ + "/FinalResolution");

    // Grab the histograms
    std::vector<TH1F *> vResPhi = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    for (int i=0; i<6; i++){
      vResPhi[i] = meresphi_vect[i]->getTH1F();
    }

    // Book the new histogram to contain the results
    MonitorElement *me_res_phi = ibooker.book1D("PhiResolution", "#phi resolution vs |#eta|", eta_binnum, eta_bins);
    TH1F *resPhi = me_res_phi->getTH1F();
    resPhi->GetXaxis()->SetTitle("tracking particle |#eta|");
    resPhi->GetYaxis()->SetTitle("#sigma(#Delta#phi)");
    resPhi->SetMinimum(0.0);
    resPhi->SetStats(false);

    for (int i = 0; i < 6; i++) {
      resPhi->SetBinContent(i + 1, vResPhi[i]->GetStdDev());
      resPhi->SetBinError(i + 1, vResPhi[i]->GetStdDevError());
    }
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for phi resolution cannot be found!\n";
  }

  if (std::find(meresVtxZ_vect.begin(), meresVtxZ_vect.end(), nullptr)==meresVtxZ_vect.end()) {
    // Set the current directoy
    igetter.setCurrentFolder(topFolderName_ + "/FinalResolution");

    // Grab the histograms
    std::vector<TH1F *> vResVtxZ = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    for (int i=0; i<6; i++){
      vResVtxZ[i] = meresVtxZ_vect[i]->getTH1F();
    }

    // Book the new histogram to contain the results
    MonitorElement *me_res_VtxZ = ibooker.book1D("VtxZResolution", "VtxZ resolution vs |#eta|", eta_binnum, eta_bins);
    TH1F *resVtxZ = me_res_VtxZ->getTH1F();
    resVtxZ->GetXaxis()->SetTitle("tracking particle |#eta|");
    resVtxZ->GetYaxis()->SetTitle("#sigma(#DeltaVtxZ) [cm]");
    resVtxZ->SetMinimum(0.0);
    resVtxZ->SetStats(false);

    for (int i = 0; i < 6; i++) {
      resVtxZ->SetBinContent(i + 1, vResVtxZ[i]->GetStdDev());
      resVtxZ->SetBinError(i + 1, vResVtxZ[i]->GetStdDevError());
    }
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for VtxZ resolution cannot be found!\n";
  }

  if (std::find(meresd0_vect.begin(), meresd0_vect.end(), nullptr)==meresd0_vect.end()) {
    // Set the current directoy
    igetter.setCurrentFolder(topFolderName_ + "/FinalResolution");

    // Grab the histograms
    std::vector<TH1F *> vResD0 = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
    for (int i=0; i<6; i++){
      vResD0[i] = meresd0_vect[i]->getTH1F();
    }

    // Book the new histogram to contain the results
    MonitorElement *me_res_d0 = ibooker.book1D("d0Resolution", "d_{0} resolution vs |#eta|", eta_binnum, eta_bins);
    TH1F *resd0 = me_res_d0->getTH1F();
    resd0->GetXaxis()->SetTitle("tracking particle |#eta|");
    resd0->GetYaxis()->SetTitle("#sigma(#Deltad_{0}) [cm]");
    resd0->SetMinimum(0.0);
    resd0->SetStats(false);

    for (int i = 0; i < 6; i++) {
      resd0->SetBinContent(i + 1, vResD0[i]->GetStdDev());
      resd0->SetBinError(i + 1, vResD0[i]->GetStdDevError());
    }
  }  // if ME found
  else {
    edm::LogWarning("DataNotFound") << "Monitor elements for d0 resolution cannot be found!\n";
  }
}  // end dqmEndJob

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
void Phase2OTHarvestReconstruction::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("TopFolderName", "TrackerPhase2OTL1TrackV");
  descriptions.add("Phase2OTHarvestReconstruction", desc);
}
DEFINE_FWK_MODULE(Phase2OTHarvestReconstruction);
