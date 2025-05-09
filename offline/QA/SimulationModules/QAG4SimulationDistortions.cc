
#include "QAG4SimulationDistortions.h"

#include <fun4all/SubsysReco.h>

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/TrackSeed.h>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TTree.h>

#include <Acts/Definitions/Algebra.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>  // for pair
#include <vector>

namespace
{

  // square
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }

  // radius
  template <class T>
  T get_r(const T& x, const T& y)
  {
    return std::sqrt(square(x) + square(y));
  }

  template <class T>
  inline constexpr T deltaPhi(const T& phi)
  {
    if (phi > M_PI)
    {
      return phi - 2. * M_PI;
    }
    else if (phi <= -M_PI)
    {
      return phi + 2. * M_PI;
    }
    else
    {
      return phi;
    }
  }

  /// return number of clusters of a given type that belong to a tracks
  template <int type>
  int count_clusters(const std::vector<TrkrDefs::cluskey>& keys)
  {
    return std::count_if(keys.begin(), keys.end(),
                         [](const TrkrDefs::cluskey& key)
                         { return TrkrDefs::getTrkrId(key) == type; });
  }
}  // namespace

//____________________________________________________________________________..
QAG4SimulationDistortions::QAG4SimulationDistortions(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
QAG4SimulationDistortions::~QAG4SimulationDistortions() = default;

//____________________________________________________________________________..
int QAG4SimulationDistortions::Init(PHCompositeNode* /*unused*/)
{

  // reset counters
  m_total_tracks = 0;
  m_accepted_tracks = 0;

  m_total_states = 0;
  m_accepted_states = 0;

  // histogram manager
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1* h(nullptr);

  h = new TH2F(TString(get_histo_prefix()) + "betadz", ";tan#beta; #Deltaz [cm]", 100, -0.5, 0.5, 100, -0.5, 0.5);

  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "alphardphi", ";tan#alpha; r#Delta#phi [cm]", 100, -0.5, 0.5, 100, -0.5, 0.5);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "rphiResid", ";r [cm]; #Deltar#phi [cm]", 60, 20, 80, 500, -2, 2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "zResid", ";z [cm]; #Deltaz [cm]", 200, -100, 100, 1000, -2, 2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "etaResid", ";#eta;#Delta#eta", 20, -1, 1, 500, -0.2, 0.2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "etaResidLayer", ";r [cm]; #Delta#eta", 60, 20, 80, 500, -0.2, 0.2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "zResidLayer", ";r [cm]; #Deltaz [cm]", 60, 20, 80, 1000, -2, 2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "deltarphi_layer", ";layer; r.#Delta#phi_{track-cluster} (cm)", 57, 0, 57, 500, -2, 2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "deltaz_layer", ";layer; #Deltaz_{track-cluster} (cm)", 57, 0, 57, 100, -2, 2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "statez_pulls", "layer; #Deltaz_{track-cluster}/#sigma_{z}^{state}", 57, 0, 57, 100, -5, 5);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "staterphi_pulls", "layer; #Deltar#phi_{track-cluster}/#sigma_{rphi}^{state}", 57, 0, 57, 100, -5, 5);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "clusz_pulls", "layer; #Deltaz_{track-cluster}/#sigma_{z}^{clus}", 57, 0, 57, 100, -5, 5);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "clusrphi_pulls", "layer; #Deltar#phi_{track-cluster}/#sigma_{r#phi}^{clus}", 57, 0, 57, 100, -5, 5);
  hm->registerHisto(h);

  TTree* t(nullptr);


  t = new TTree(TString(get_histo_prefix()) + "residTree", "tpc residual info");
  t->Branch("tanAlpha", &m_tanAlpha, "tanAlpha/F");
  t->Branch("tanBeta", &m_tanBeta, "tanBeta/F");
  t->Branch("drphi", &m_drphi, "drphi/F");
  t->Branch("dz", &m_dz, "dz/F");
  t->Branch("clusR", &m_clusR, "clusR/F");
  t->Branch("clusPhi", &m_clusPhi, "clusPhi/F");
  t->Branch("clusZ", &m_clusZ, "clusZ/F");
  t->Branch("statePhi", &m_statePhi, "statePhi/F");
  t->Branch("stateZ", &m_stateZ, "stateZ/F");
  t->Branch("stateR", &m_stateR, "stateR/F");
  t->Branch("stateRPhiErr", &m_stateRPhiErr, "stateRPhiErr/F");
  t->Branch("stateZErr", &m_stateZErr, "stateZErr/F");
  t->Branch("clusRPhiErr", &m_clusRPhiErr, "clusRPhiErr/F");
  t->Branch("clusZErr", &m_clusZErr, "clusZErr/F");
  t->Branch("cluskey", &m_cluskey, "cluskey/l");
  t->Branch("event", &m_event, "event/I");

  hm->registerHisto(t);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int QAG4SimulationDistortions::InitRun(PHCompositeNode* topNode)
{
  // track map
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackmapname);

  // cluster map
  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  // load geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  // load distortion corrections
  m_globalPositionWrapper.loadNodes(topNode);


  if (!m_trackMap || !m_clusterContainer || !m_tGeometry)
  {
    std::cout << PHWHERE << "Necessary distortion container not on node tree. Bailing."
              << std::endl;

    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int QAG4SimulationDistortions::process_event(PHCompositeNode* /*unused*/)
{
  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  auto h_beta = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "betadz"));
  assert(h_beta);

  auto h_alpha = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "alphardphi"));
  assert(h_alpha);

  auto h_rphiResid = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "rphiResid"));
  assert(h_rphiResid);

  auto h_zResid = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "zResid"));
  assert(h_zResid);

  auto h_etaResid = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "etaResid"));
  assert(h_etaResid);

  auto h_etaResidLayer = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "etaResidLayer"));
  assert(h_etaResidLayer);

  auto h_zResidLayer = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "zResidLayer"));
  assert(h_zResidLayer);

  auto h_deltarphi_layer = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "deltarphi_layer"));
  assert(h_deltarphi_layer);

  auto h_deltaz_layer = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "deltaz_layer"));
  assert(h_deltaz_layer);

  auto h_statez_pulls = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "statez_pulls"));
  assert(h_statez_pulls);

  auto h_staterphi_pulls = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "staterphi_pulls"));
  assert(h_staterphi_pulls);

  auto h_clusz_pulls = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "clusz_pulls"));
  assert(h_clusz_pulls);

  auto h_clusrphi_pulls = dynamic_cast<TH1*>(hm->getHisto(get_histo_prefix() + "clusrphi_pulls"));
  assert(h_clusrphi_pulls);

  auto t_tree = dynamic_cast<TTree*>(hm->getHisto(get_histo_prefix() + "residTree"));
  assert(t_tree);

  std::cout << "QAG4SimulationDistortions::process_event - tracks: " << m_trackMap->size() << std::endl;
  for (const auto& [key, track] : *m_trackMap)
  {

    // total track counter
    ++m_total_tracks;

    // get track crossing and check
    const auto crossing = track->get_crossing();
    if(crossing == SHRT_MAX)
    {
      std::cout << "QAG4SimulationDistortions::process_event - invalid crossing. Track skipped." << std::endl;
      continue;
    }

    // check track quality
    if (!checkTrack(track))
    {
      continue;
    }

    // get seeeds
    auto tpcSeed = track->get_tpc_seed();
    auto siliconSeed = track->get_silicon_seed();

    /// Should have never been added to the map...
    if (!tpcSeed || !siliconSeed)
    {
      continue;
    }

    // accepted track counter
    ++m_accepted_tracks;

    for (auto iter = track->begin_states(); iter != track->end_states(); ++iter)
    {

      ++m_total_states;

      auto& state = iter->second;
      const auto ckey = state->get_cluskey();
      const auto trkrId = TrkrDefs::getTrkrId(ckey);

      if( trkrId != TrkrDefs::tpcId )
      { continue; }

      ++m_accepted_states;

      auto cluster = m_clusterContainer->findCluster(ckey);

      const auto clusGlobPosition = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cluster, crossing);

      const float clusR = get_r(clusGlobPosition(0), clusGlobPosition(1));
      const float clusPhi = std::atan2(clusGlobPosition(1), clusGlobPosition(0));
      const float clusZ = clusGlobPosition(2);

      // cluster errors
      const float clusRPhiErr = cluster->getRPhiError();
      const float clusZErr = cluster->getZError();

      const Acts::Vector3 stateGlobPosition = Acts::Vector3(state->get_x(),
                                                            state->get_y(),
                                                            state->get_z());
      const Acts::Vector3 stateGlobMom = Acts::Vector3(state->get_px(),
                                                       state->get_py(),
                                                       state->get_pz());

      const float stateRPhiErr = state->get_rphi_error();
      const float stateZErr = state->get_z_error();

      const float stateR = get_r(stateGlobPosition(0), stateGlobPosition(1));

      const auto dr = clusR - stateR;
      const auto trackDrDt = (stateGlobPosition(0) * stateGlobMom(0) + stateGlobPosition(1) * stateGlobMom(1)) / stateR;
      const auto trackDxDr = stateGlobMom(0) / trackDrDt;
      const auto trackDyDr = stateGlobMom(1) / trackDrDt;
      const auto trackDzDr = stateGlobMom(2) / trackDrDt;

      const auto trackX = stateGlobPosition(0) + dr * trackDxDr;
      const auto trackY = stateGlobPosition(1) + dr * trackDyDr;
      const auto trackZ = stateGlobPosition(2) + dr * trackDzDr;
      const float statePhi = std::atan2(trackY, trackX);
      const float stateZ = trackZ;

      // Calculate residuals
      const float drphi = clusR * deltaPhi(clusPhi - statePhi);
      const float dz = clusZ - stateZ;

      const auto trackPPhi = -stateGlobMom(0) * std::sin(statePhi) + stateGlobMom(1) * std::cos(statePhi);
      const auto trackPR = stateGlobMom(0) * std::cos(statePhi) + stateGlobMom(1) * std::sin(statePhi);
      const auto trackPZ = stateGlobMom(2);

      const auto trackAlpha = -trackPPhi / trackPR;
      const auto trackBeta = -trackPZ / trackPR;
      const auto trackEta = std::atanh(stateGlobMom(2) / stateGlobMom.norm());
      const auto clusEta = std::atanh(clusZ / clusGlobPosition.norm());

      h_alpha->Fill(trackAlpha, drphi);
      h_beta->Fill(trackBeta, dz);
      h_rphiResid->Fill(clusR, drphi);
      h_zResid->Fill(stateZ, dz);
      h_etaResid->Fill(trackEta, clusEta - trackEta);
      h_zResidLayer->Fill(clusR, dz);
      h_etaResidLayer->Fill(clusR, clusEta - trackEta);

      const auto layer = TrkrDefs::getLayer(ckey);
      h_deltarphi_layer->Fill(layer, drphi);
      h_deltaz_layer->Fill(layer, dz);

      h_statez_pulls->Fill(layer, dz / stateZErr);
      h_staterphi_pulls->Fill(layer, drphi / stateRPhiErr);
      h_clusz_pulls->Fill(layer, dz / clusZErr);
      h_clusrphi_pulls->Fill(layer, drphi / clusRPhiErr);

      m_tanAlpha = trackAlpha;
      m_tanBeta = trackBeta;
      m_drphi = drphi;
      m_dz = dz;
      m_clusR = clusR;
      m_clusPhi = clusPhi;
      m_clusZ = clusZ;
      m_statePhi = statePhi;
      m_stateZ = stateZ;
      m_stateR = stateR;
      m_stateRPhiErr = stateRPhiErr;
      m_stateZErr = stateZErr;
      m_clusRPhiErr = clusRPhiErr;
      m_clusZErr = clusZErr;
      m_cluskey = ckey;
      t_tree->Fill();
    }
  }

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
int QAG4SimulationDistortions::End(PHCompositeNode* /*topNode*/)
{
  // print counters
  std::cout
      << "QAG4SimulationDistortions::End -"
      << " track statistics total: " << m_total_tracks
      << " accepted: " << m_accepted_tracks
      << " fraction: " << 100. * m_accepted_tracks / m_total_tracks << "%"
      << std::endl;

  std::cout
      << "QAG4SimulationDistortions::End -"
      << " state statistics total: " << m_total_states
      << " accepted: " << m_accepted_states << " fraction: "
      << 100. * m_accepted_states / m_total_states << "%"
      << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________________
bool QAG4SimulationDistortions::checkTrack(SvtxTrack* track)
{

  if (track->get_pt() < 0.5)
  {
    return false;
  }

  // ignore tracks with too few mvtx, intt and micromegas hits
  const auto cluster_keys(get_cluster_keys(track));
  if (count_clusters<TrkrDefs::mvtxId>(cluster_keys) < 2)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::inttId>(cluster_keys) < 2)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::micromegasId>(cluster_keys) < 2)
  {
    return false;
  }

  return true;
}

std::vector<TrkrDefs::cluskey> QAG4SimulationDistortions::get_cluster_keys(SvtxTrack* track)
{
  std::vector<TrkrDefs::cluskey> out;
  for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
  {
    if (seed)
    {
      std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
    }
  }
  return out;
}
