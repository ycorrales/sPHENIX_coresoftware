#include "KFParticle_truthAndDetTools.h"

#include "KFParticle_Tools.h"  // for KFParticle_Tools

#include <g4eval/SvtxEvalStack.h>   // for SvtxEvalStack
#include <g4eval/SvtxTrackEval.h>   // for SvtxTrackEval
#include <g4eval/SvtxTruthEval.h>   // for SvtxTruthEval
#include <g4eval/SvtxVertexEval.h>  // for SvtxVertexEval

#include <globalvertex/SvtxVertex.h>         // for SvtxVertex
#include <globalvertex/SvtxVertexMap.h>      // for SvtxVertexMap, SvtxVer...
#include <trackbase/InttDefs.h>              // for getLadderPhiId, getLad...
#include <trackbase/MvtxDefs.h>              // for getChipId, getStaveId
#include <trackbase/TpcDefs.h>               // for getSectorId, getSide
#include <trackbase/TrkrCluster.h>           // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>  // for TrkrClusterContainer
#include <trackbase/TrkrDefs.h>              // for getLayer, getTrkrId
#include <trackbase_historic/SvtxPHG4ParticleMap.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::...
#include <trackbase_historic/SvtxTrackMap.h>  // for SvtxTrackMap, SvtxTrac...

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <g4main/PHG4Particle.h>            // for PHG4Particle
#include <g4main/PHG4TruthInfoContainer.h>  // for PHG4TruthInfoContainer
#include <g4main/PHG4VtxPoint.h>            // for PHG4VtxPoint

#include <phhepmc/PHHepMCGenEvent.h>     // for PHHepMCGenEvent
#include <phhepmc/PHHepMCGenEventMap.h>  // for PHHepMCGenEventMap

#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>        // for getClass

#include <KFParticle.h>  // for KFParticle
#include <TString.h>     // for TString, operator+
#include <TTree.h>       // for TTree

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>   // for GenEvent::particle_con...
#include <HepMC/GenVertex.h>  // for GenVertex::particle_it...
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h>    // for GenParticle
#include <HepMC/IteratorRange.h>  // for parents
#include <HepMC/SimpleVector.h>   // for FourVector

#include <algorithm>  // for max, find
#include <cmath>      // for pow, sqrt
#include <cstdlib>    // for NULL, abs
#include <iostream>   // for operator<<, endl, basi...
#include <iterator>   // for end, begin
#include <map>        // for _Rb_tree_iterator, map
#include <memory>     // for allocator_traits<>::va...
#include <utility>    // for pair

class PHNode;

std::map<std::string, int> Use =
    {
        {"MVTX", 1},
        {"INTT", 1},
        {"TPC", 1},
        {"TPOT", 1},
        {"EMCAL", 0},
        {"OHCAL", 0},
        {"IHCAL", 0}};

KFParticle_truthAndDetTools::KFParticle_truthAndDetTools()
  : m_svtx_evalstack(nullptr)
{
}  // Constructor

KFParticle_truthAndDetTools::~KFParticle_truthAndDetTools() = default;  // Destructor

SvtxTrack *KFParticle_truthAndDetTools::getTrack(unsigned int track_id, SvtxTrackMap *trackmap)
{
  SvtxTrack *matched_track = nullptr;

  for (auto &iter : *trackmap)
  {
    if (iter.first == track_id)
    {
      matched_track = iter.second;
    }
  }

  return matched_track;
}

GlobalVertex *KFParticle_truthAndDetTools::getVertex(unsigned int vertex_id, GlobalVertexMap *vertexmap)
{
  GlobalVertex *matched_vertex = vertexmap->get(vertex_id);

  return matched_vertex;
}

PHG4Particle *KFParticle_truthAndDetTools::getTruthTrack(SvtxTrack *thisTrack, PHCompositeNode *topNode)
{
  /*
   * There are two methods for getting the truth rack from the reco track
   * 1. (recommended) Use the reco -> truth tables (requires SvtxPHG4ParticleMap). Introduced Summer of 2022
   * 2. Get truth track via nClusters. Older method and will work with older DSTs
   */

  PHG4Particle *particle = nullptr;

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("SvtxPHG4ParticleMap"));
  if (findNode)
  {
    findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("G4TruthInfo"));
    if (findNode)
    {
      m_truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    }
    else
    {
      std::cout << "KFParticle truth matching: G4TruthInfo does not exist" << std::endl;
    }

    SvtxPHG4ParticleMap *dst_reco_truth_map = findNode::getClass<SvtxPHG4ParticleMap>(topNode, "SvtxPHG4ParticleMap");
    std::map<float, std::set<int>> truth_set = dst_reco_truth_map->get(thisTrack->get_id());
    if (truth_set.size() > 0)
    {
      std::pair<float, std::set<int>> best_weight = *truth_set.rbegin();
      int best_truth_id = *best_weight.second.rbegin();
      particle = m_truthinfo->GetParticle(best_truth_id);
    }
  }
  else
  {
    std::cout << __FILE__ << ": SvtxPHG4ParticleMap not found, reverting to max_truth_particle_by_nclusters()" << std::endl;

    if (!m_svtx_evalstack)
    {
      m_svtx_evalstack = new SvtxEvalStack(topNode);
      // clustereval = m_svtx_evalstack->get_cluster_eval();
      // hiteval = m_svtx_evalstack->get_hit_eval();
      trackeval = m_svtx_evalstack->get_track_eval();
      trutheval = m_svtx_evalstack->get_truth_eval();
      vertexeval = m_svtx_evalstack->get_vertex_eval();
    }

    m_svtx_evalstack->next_event(topNode);

    particle = trackeval->max_truth_particle_by_nclusters(thisTrack);
  }
  return particle;
}

void KFParticle_truthAndDetTools::initializeTruthBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number, bool m_constrain_to_vertex_truthMatch)
{
  m_tree->Branch(TString(daughter_number) + "_true_ID", &m_true_daughter_id[daughter_id], TString(daughter_number) + "_true_ID/I");
  if (m_constrain_to_vertex_truthMatch)
  {
    m_tree->Branch(TString(daughter_number) + "_true_IP", &m_true_daughter_ip[daughter_id], TString(daughter_number) + "_true_IP/F");
    m_tree->Branch(TString(daughter_number) + "_true_IP_xy", &m_true_daughter_ip_xy[daughter_id], TString(daughter_number) + "_true_IP_xy/F");
  }
  m_tree->Branch(TString(daughter_number) + "_true_px", &m_true_daughter_px[daughter_id], TString(daughter_number) + "_true_px/F");
  m_tree->Branch(TString(daughter_number) + "_true_py", &m_true_daughter_py[daughter_id], TString(daughter_number) + "_true_py/F");
  m_tree->Branch(TString(daughter_number) + "_true_pz", &m_true_daughter_pz[daughter_id], TString(daughter_number) + "_true_pz/F");
  m_tree->Branch(TString(daughter_number) + "_true_p", &m_true_daughter_p[daughter_id], TString(daughter_number) + "_true_p/F");
  m_tree->Branch(TString(daughter_number) + "_true_pT", &m_true_daughter_pt[daughter_id], TString(daughter_number) + "_true_pT/F");
  m_tree->Branch(TString(daughter_number) + "_true_EV_x", &m_true_daughter_vertex_x[daughter_id], TString(daughter_number) + "_true_EV_x/F");
  m_tree->Branch(TString(daughter_number) + "_true_EV_y", &m_true_daughter_vertex_y[daughter_id], TString(daughter_number) + "_true_EV_y/F");
  m_tree->Branch(TString(daughter_number) + "_true_EV_z", &m_true_daughter_vertex_z[daughter_id], TString(daughter_number) + "_true_EV_z/F");
  if (m_constrain_to_vertex_truthMatch)
  {
    m_tree->Branch(TString(daughter_number) + "_true_PV_x", &m_true_daughter_pv_x[daughter_id], TString(daughter_number) + "_true_PV_x/F");
    m_tree->Branch(TString(daughter_number) + "_true_PV_y", &m_true_daughter_pv_y[daughter_id], TString(daughter_number) + "_true_PV_y/F");
    m_tree->Branch(TString(daughter_number) + "_true_PV_z", &m_true_daughter_pv_z[daughter_id], TString(daughter_number) + "_true_PV_z/F");
  }
  m_tree->Branch(TString(daughter_number) + "_true_track_history_PDG_ID", &m_true_daughter_track_history_PDG_ID[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_PDG_mass", &m_true_daughter_track_history_PDG_mass[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_px", &m_true_daughter_track_history_px[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_py", &m_true_daughter_track_history_py[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_pz", &m_true_daughter_track_history_pz[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_pE", &m_true_daughter_track_history_pE[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_pT", &m_true_daughter_track_history_pT[daughter_id]);
}

void KFParticle_truthAndDetTools::fillTruthBranch(PHCompositeNode *topNode, TTree * /*m_tree*/, const KFParticle &daughter, int daughter_id, const KFParticle &kfvertex, bool m_constrain_to_vertex_truthMatch)
{
  float true_px, true_py, true_pz, true_p, true_pt;

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_trk_map_node_name_nTuple));

  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }

  if (m_use_mbd_vertex_truth)
  {
    findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("MbdVertexMap"));
    if (findNode)
    {
      dst_mbdvertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
    }
    else
    {
      std::cout << "KFParticle truth matching: " << m_vtx_map_node_name_nTuple << " does not exist" << std::endl;
    }
  }
  else
  {
    findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_vtx_map_node_name_nTuple));
    if (findNode)
    {
      dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vtx_map_node_name_nTuple);
    }
    else
    {
      std::cout << "KFParticle truth matching: " << m_vtx_map_node_name_nTuple << " does not exist" << std::endl;
    }
  }

  auto globalvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!globalvertexmap)
  {
    std::cout << "KFParticle truth matching: GlobalVertexMap does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);
  g4particle = getTruthTrack(track, topNode);

  bool isParticleValid = g4particle == nullptr ? false : true;

  true_px = isParticleValid ? (Float_t) g4particle->get_px() : 0.;
  true_py = isParticleValid ? (Float_t) g4particle->get_py() : 0.;
  true_pz = isParticleValid ? (Float_t) g4particle->get_pz() : 0.;
  true_p = sqrt(pow(true_px, 2) + pow(true_py, 2) + pow(true_pz, 2));
  true_pt = sqrt(pow(true_px, 2) + pow(true_py, 2));

  m_true_daughter_px[daughter_id] = true_px;
  m_true_daughter_py[daughter_id] = true_py;
  m_true_daughter_pz[daughter_id] = true_pz;
  m_true_daughter_p[daughter_id] = true_p;
  m_true_daughter_pt[daughter_id] = true_pt;
  m_true_daughter_id[daughter_id] = isParticleValid ? g4particle->get_pid() : 0;

  if (!m_svtx_evalstack)
  {
    m_svtx_evalstack = new SvtxEvalStack(topNode);
    // clustereval = m_svtx_evalstack->get_cluster_eval();
    // hiteval = m_svtx_evalstack->get_hit_eval();
    trackeval = m_svtx_evalstack->get_track_eval();
    trutheval = m_svtx_evalstack->get_truth_eval();
    vertexeval = m_svtx_evalstack->get_vertex_eval();
  }

  if (isParticleValid)
  {
    g4vertex_point = trutheval->get_vertex(g4particle);
  }

  m_true_daughter_vertex_x[daughter_id] = isParticleValid ? g4vertex_point->get_x() : 0.;
  m_true_daughter_vertex_y[daughter_id] = isParticleValid ? g4vertex_point->get_y() : 0.;
  m_true_daughter_vertex_z[daughter_id] = isParticleValid ? g4vertex_point->get_z() : 0.;

  if (m_constrain_to_vertex_truthMatch)
  {
    // Calculate true DCA
    GlobalVertex *recoVertex = getVertex(kfvertex.Id(), globalvertexmap);
    GlobalVertex::VTXTYPE whichVtx = m_use_mbd_vertex_truth ? GlobalVertex::MBD : GlobalVertex::SVTX;
    auto svtxviter = recoVertex->find_vertexes(whichVtx);

    // check that it contains a track vertex
    if (svtxviter == recoVertex->end_vertexes())
    {
      std::string vtxType = m_use_mbd_vertex_truth ? "MBD" : "silicon";
      std::cout << "Have a global vertex with no " << vtxType << " vertex... shouldn't happen in KFParticle_truthAndDetTools::fillTruthBranch..." << std::endl;
    }

    auto svtxvertexvector = svtxviter->second;
    MbdVertex *mbdvertex = nullptr;
    SvtxVertex *svtxvertex = nullptr;

    for (auto &vertex_iter : svtxvertexvector)
    {
      if (m_use_mbd_vertex_truth)
      {
        mbdvertex = dst_mbdvertexmap->find(vertex_iter->get_id())->second;
      }
      else
      {
        svtxvertex = dst_vertexmap->find(vertex_iter->get_id())->second;
      }
    }

    PHG4VtxPoint *truePoint = nullptr;
    if (m_use_mbd_vertex_truth)
    {
      std::set<PHG4VtxPoint*> truePointSet = vertexeval->all_truth_points(mbdvertex);
      truePoint = *truePointSet.begin();
    }
    else
    {
      truePoint = vertexeval->max_truth_point_by_ntracks(svtxvertex);
    }

    if (truePoint == nullptr && isParticleValid)
    {
      //PHG4Particle *g4mother = m_truthinfo->GetParticle(g4particle->get_parent_id());
      PHG4Particle *g4mother = m_truthinfo->GetPrimaryParticle(g4particle->get_parent_id());
      if (!g4mother)
      {
        std::cout << "KFParticle truth matching: True mother not found!\n";
        std::cout << "Your truth track DCA will be measured wrt a reconstructed vertex!" << std::endl;
        truePoint = nullptr;
      }
      else
      {
        truePoint = m_truthinfo->GetVtx(g4mother->get_vtx_id());  // Note, this may not be the PV for a decay with tertiaries
      }
    }

    KFParticle trueKFParticleVertex;

    float f_vertexParameters[6] = {0};

    if (truePoint == nullptr)
    {
      std::cout << "KFParticle truth matching: This event has no PHG4VtxPoint information!\n";
      std::cout << "Your truth track DCA will be measured wrt a reconstructed vertex!" << std::endl;

      f_vertexParameters[0] = recoVertex->get_x();
      f_vertexParameters[1] = recoVertex->get_y();
      f_vertexParameters[2] = recoVertex->get_z();
    }
    else
    {
      f_vertexParameters[0] = truePoint->get_x();
      f_vertexParameters[1] = truePoint->get_y();
      f_vertexParameters[2] = truePoint->get_z();
    }

    float f_vertexCovariance[21] = {0};

    trueKFParticleVertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);

    KFParticle trueKFParticle;

    float f_trackParameters[6] = {m_true_daughter_vertex_x[daughter_id],
                                  m_true_daughter_vertex_y[daughter_id],
                                  m_true_daughter_vertex_z[daughter_id],
                                  true_px,
                                  true_py,
                                  true_pz};

    float f_trackCovariance[21] = {0};

    trueKFParticle.Create(f_trackParameters, f_trackCovariance, 1, -1);

    m_true_daughter_ip[daughter_id] = trueKFParticle.GetDistanceFromVertex(trueKFParticleVertex);
    m_true_daughter_ip_xy[daughter_id] = trueKFParticle.GetDistanceFromVertexXY(trueKFParticleVertex);

    m_true_daughter_pv_x[daughter_id] = truePoint == nullptr ? -99. : truePoint->get_x();
    m_true_daughter_pv_y[daughter_id] = truePoint == nullptr ? -99. : truePoint->get_y();
    m_true_daughter_pv_z[daughter_id] = truePoint == nullptr ? -99. : truePoint->get_z();
  }
}

void KFParticle_truthAndDetTools::fillGeant4Branch(PHG4Particle *particle, int daughter_id)
{
  Float_t pT = sqrt(pow(particle->get_px(), 2) + pow(particle->get_py(), 2));

  m_true_daughter_track_history_PDG_ID[daughter_id].push_back(particle->get_pid());
  m_true_daughter_track_history_PDG_mass[daughter_id].push_back(0);
  m_true_daughter_track_history_px[daughter_id].push_back((Float_t) particle->get_px());
  m_true_daughter_track_history_py[daughter_id].push_back((Float_t) particle->get_py());
  m_true_daughter_track_history_pz[daughter_id].push_back((Float_t) particle->get_pz());
  m_true_daughter_track_history_pE[daughter_id].push_back((Float_t) particle->get_e());
  m_true_daughter_track_history_pT[daughter_id].push_back((Float_t) pT);
}

void KFParticle_truthAndDetTools::fillHepMCBranch(HepMC::GenParticle *particle, int daughter_id)
{
  const HepMC::FourVector &myFourVector = particle->momentum();

  m_true_daughter_track_history_PDG_ID[daughter_id].push_back(particle->pdg_id());
  m_true_daughter_track_history_PDG_mass[daughter_id].push_back((Float_t) particle->generatedMass());
  m_true_daughter_track_history_px[daughter_id].push_back((Float_t) myFourVector.px());
  m_true_daughter_track_history_py[daughter_id].push_back((Float_t) myFourVector.py());
  m_true_daughter_track_history_pz[daughter_id].push_back((Float_t) myFourVector.pz());
  m_true_daughter_track_history_pE[daughter_id].push_back((Float_t) myFourVector.e());
  m_true_daughter_track_history_pT[daughter_id].push_back((Float_t) myFourVector.perp());
}

int KFParticle_truthAndDetTools::getHepMCInfo(PHCompositeNode *topNode, TTree * /*m_tree*/, const KFParticle &daughter, int daughter_id)
{
  // Make dummy particle for null pointers and missing nodes
  HepMC::GenParticle *dummyParticle = new HepMC::GenParticle();
  HepMC::FourVector dummyFourVector(0, 0, 0, 0);
  dummyParticle->set_momentum(dummyFourVector);
  dummyParticle->set_pdg_id(0);
  dummyParticle->set_generated_mass(0.);

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_trk_map_node_name_nTuple));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);
  g4particle = getTruthTrack(track, topNode);

  bool isParticleValid = g4particle == nullptr ? false : true;

  if (!isParticleValid)
  {
    std::cout << "KFParticle truth matching: this track is a ghost" << std::endl;
    fillHepMCBranch(dummyParticle, daughter_id);
    return 0;
  }

  m_geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!m_geneventmap)
  {
    std::cout << "KFParticle truth matching: Missing node PHHepMCGenEventMap" << std::endl;
    std::cout << "You will have no mother information" << std::endl;
    fillHepMCBranch(dummyParticle, daughter_id);
    return 0;
  }

  m_genevt = m_geneventmap->get(1);
  if (!m_genevt)
  {
    std::cout << "KFParticle truth matching: Missing node PHHepMCGenEvent" << std::endl;
    std::cout << "You will have no mother information" << std::endl;
    fillHepMCBranch(dummyParticle, daughter_id);
    return 0;
  }

  // Start by looking for our particle in the Geant record
  // Any decay that Geant4 handles will not be in the HepMC record
  // This can happen if you limit the decay volume in the generator
  if (g4particle->get_parent_id() != 0)
  {
    findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("G4TruthInfo"));
    if (findNode)
    {
      m_truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    }
    else
    {
      std::cout << "KFParticle truth matching: G4TruthInfo does not exist" << std::endl;
    }
    while (g4particle->get_parent_id() != 0)
    {
      g4particle = m_truthinfo->GetParticle(g4particle->get_parent_id());
      fillGeant4Branch(g4particle, daughter_id);
    }
  }

  HepMC::GenEvent *theEvent = m_genevt->getEvent();
  HepMC::GenParticle *prevParticle = nullptr;

  // int forbiddenPDGIDs[] = {21, 22};  //Stop tracing history when we reach quarks, gluons and photons
  int forbiddenPDGIDs[] = {0};  // 20230921 - Request made to have gluon information to see about gluon-splitting

  for (HepMC::GenEvent::particle_const_iterator p = theEvent->particles_begin(); p != theEvent->particles_end(); ++p)
  {
    if (((*p)->barcode() == g4particle->get_barcode()))
    {
      prevParticle = *p;
      while (!prevParticle->is_beam())
      {
        bool breakOut = false;
        for (HepMC::GenVertex::particle_iterator mother = prevParticle->production_vertex()->particles_begin(HepMC::parents);
             mother != prevParticle->production_vertex()->particles_end(HepMC::parents); ++mother)
        {
          if (std::find(std::begin(forbiddenPDGIDs), std::end(forbiddenPDGIDs),
                        abs((*mother)->pdg_id())) != std::end(forbiddenPDGIDs))
          {
            breakOut = true;
            break;
          }

          fillHepMCBranch((*mother), daughter_id);
          prevParticle = *mother;
        }
        if (breakOut)
        {
          break;
        }
      }
    }
  }

  return 0;
}  // End of function

void KFParticle_truthAndDetTools::initializeCaloBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number)
{
  m_tree->Branch(TString(daughter_number) + "_EMCAL_DeltaPhi", &detector_emcal_deltaphi[daughter_id], TString(daughter_number) + "_EMCAL_DeltaPhi/F");
  m_tree->Branch(TString(daughter_number) + "_EMCAL_DeltaEta", &detector_emcal_deltaeta[daughter_id], TString(daughter_number) + "_EMCAL_DeltaEta/F");
  m_tree->Branch(TString(daughter_number) + "_EMCAL_energy_3x3", &detector_emcal_energy_3x3[daughter_id], TString(daughter_number) + "_EMCAL_energy_3x3/F");
  m_tree->Branch(TString(daughter_number) + "_EMCAL_energy_5x5", &detector_emcal_energy_5x5[daughter_id], TString(daughter_number) + "_EMCAL_energy_5x5/F");
  m_tree->Branch(TString(daughter_number) + "_EMCAL_energy_cluster", &detector_emcal_cluster_energy[daughter_id], TString(daughter_number) + "_EMCAL_energy_cluster/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_DeltaPhi", &detector_ihcal_deltaphi[daughter_id], TString(daughter_number) + "_IHCAL_DeltaPhi/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_DeltaEta", &detector_ihcal_deltaeta[daughter_id], TString(daughter_number) + "_IHCAL_DeltaEta/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_energy_3x3", &detector_ihcal_energy_3x3[daughter_id], TString(daughter_number) + "_IHCAL_energy_3x3/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_energy_5x5", &detector_ihcal_energy_5x5[daughter_id], TString(daughter_number) + "_IHCAL_energy_5x5/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_energy_cluster", &detector_ihcal_cluster_energy[daughter_id], TString(daughter_number) + "_IHCAL_energy_cluster/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_DeltaPhi", &detector_ohcal_deltaphi[daughter_id], TString(daughter_number) + "_OHCAL_DeltaEta/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_DeltaEta", &detector_ohcal_deltaeta[daughter_id], TString(daughter_number) + "_OHCAL_DeltaEta/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_energy_3x3", &detector_ohcal_energy_3x3[daughter_id], TString(daughter_number) + "_OHCAL_energy_3x3/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_energy_5x5", &detector_ohcal_energy_5x5[daughter_id], TString(daughter_number) + "_OHCAL_energy_5x5/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_energy_cluster", &detector_ohcal_cluster_energy[daughter_id], TString(daughter_number) + "_OHCAL_energy_cluster/F");
}

void KFParticle_truthAndDetTools::fillCaloBranch(PHCompositeNode *topNode,
                                                 TTree * /*m_tree*/, const KFParticle &daughter, int daughter_id)
{
  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_trk_map_node_name_nTuple));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);

  detector_emcal_deltaphi[daughter_id] = track->get_cal_dphi(SvtxTrack::CAL_LAYER(1));
  detector_emcal_deltaeta[daughter_id] = track->get_cal_deta(SvtxTrack::CAL_LAYER(1));
  detector_emcal_energy_3x3[daughter_id] = track->get_cal_energy_3x3(SvtxTrack::CAL_LAYER(1));
  detector_emcal_energy_5x5[daughter_id] = track->get_cal_energy_5x5(SvtxTrack::CAL_LAYER(1));
  detector_emcal_cluster_energy[daughter_id] = track->get_cal_cluster_e(SvtxTrack::CAL_LAYER(1));

  detector_ihcal_deltaphi[daughter_id] = track->get_cal_dphi(SvtxTrack::CAL_LAYER(2));
  detector_ihcal_deltaeta[daughter_id] = track->get_cal_deta(SvtxTrack::CAL_LAYER(2));
  detector_ihcal_energy_3x3[daughter_id] = track->get_cal_energy_3x3(SvtxTrack::CAL_LAYER(2));
  detector_ihcal_energy_5x5[daughter_id] = track->get_cal_energy_5x5(SvtxTrack::CAL_LAYER(2));
  detector_ihcal_cluster_energy[daughter_id] = track->get_cal_cluster_e(SvtxTrack::CAL_LAYER(2));

  detector_ohcal_deltaphi[daughter_id] = track->get_cal_dphi(SvtxTrack::CAL_LAYER(3));
  detector_ohcal_deltaeta[daughter_id] = track->get_cal_deta(SvtxTrack::CAL_LAYER(3));
  detector_ohcal_energy_3x3[daughter_id] = track->get_cal_energy_3x3(SvtxTrack::CAL_LAYER(3));
  detector_ohcal_energy_5x5[daughter_id] = track->get_cal_energy_5x5(SvtxTrack::CAL_LAYER(3));
  detector_ohcal_cluster_energy[daughter_id] = track->get_cal_cluster_e(SvtxTrack::CAL_LAYER(3));
}

void KFParticle_truthAndDetTools::initializeDetectorBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number)
{
  m_tree->Branch(TString(daughter_number) + "_residual_x", &residual_x[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_residual_y", &residual_y[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_residual_z", &residual_z[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_layer", &detector_layer[daughter_id]);

  for (auto const &subdetector : Use)
  {
    if (subdetector.second)
    {
      initializeSubDetectorBranches(m_tree, subdetector.first, daughter_id, daughter_number);
    }
  }
}

void KFParticle_truthAndDetTools::initializeSubDetectorBranches(TTree *m_tree, const std::string &detectorName, int daughter_id, const std::string &daughter_number)
{
  if (detectorName == "MVTX")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_staveID", &mvtx_staveID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_chipID", &mvtx_chipID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_nHits", &detector_nHits_MVTX[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_nStates", &detector_nStates_MVTX[daughter_id]);
  }
  if (detectorName == "INTT")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_ladderZID", &intt_ladderZID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_ladderPhiID", &intt_ladderPhiID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_nHits", &detector_nHits_INTT[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_nStates", &detector_nStates_INTT[daughter_id]);
  }
  if (detectorName == "TPC")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_sectorID", &tpc_sectorID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_side", &tpc_side[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_nHits", &detector_nHits_TPC[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_nStates", &detector_nStates_TPC[daughter_id]);
  }
  if (detectorName == "TPOT")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_nHits", &detector_nHits_TPOT[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_nStates", &detector_nStates_TPOT[daughter_id]);
  }
}

void KFParticle_truthAndDetTools::fillDetectorBranch(PHCompositeNode *topNode,
                                                     TTree * /*m_tree*/, const KFParticle &daughter, int daughter_id)
{
  PHNodeIterator nodeIter(topNode);

  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_trk_map_node_name_nTuple));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }

  std::string geoName = "ActsGeometry";
  findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(geoName));
  if (findNode)
  {
    geometry = findNode::getClass<ActsGeometry>(topNode, geoName);
  }
  else
  {
    std::cout << "KFParticle detector info: " << geoName << " does not exist" << std::endl;
  }

  findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("TRKR_CLUSTER"));
  if (findNode)
  {
    dst_clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }
  else
  {
    std::cout << "KFParticle detector info: TRKR_CLUSTER does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);
  detector_nHits_MVTX[daughter_id] = 0;
  detector_nHits_INTT[daughter_id] = 0;
  detector_nHits_TPC[daughter_id] = 0;
  detector_nHits_TPOT[daughter_id] = 0;

  TrackSeed *silseed = track->get_silicon_seed();
  TrackSeed *tpcseed = track->get_tpc_seed();

  if (silseed)
  {
    for (auto cluster_iter = silseed->begin_cluster_keys(); cluster_iter != silseed->end_cluster_keys(); ++cluster_iter)
    {
      const auto &cluster_key = *cluster_iter;
      const auto trackerID = TrkrDefs::getTrkrId(cluster_key);
  
      detector_layer[daughter_id].push_back(TrkrDefs::getLayer(cluster_key));
  
      unsigned int staveId, chipId, ladderZId, ladderPhiId, sectorId, side;
      staveId = chipId = ladderZId = ladderPhiId = sectorId = side = std::numeric_limits<unsigned int>::quiet_NaN();
  
      if (Use["MVTX"] && trackerID == TrkrDefs::mvtxId)
      {
        staveId = MvtxDefs::getStaveId(cluster_key);
        chipId = MvtxDefs::getChipId(cluster_key);
        ++detector_nHits_MVTX[daughter_id];
      }
      else if (Use["INTT"] && trackerID == TrkrDefs::inttId)
      {
        ladderZId = InttDefs::getLadderZId(cluster_key);
        ladderPhiId = InttDefs::getLadderPhiId(cluster_key);
        ++detector_nHits_INTT[daughter_id];
      }
  
      mvtx_staveID[daughter_id].push_back(staveId);
      mvtx_chipID[daughter_id].push_back(chipId);
      intt_ladderZID[daughter_id].push_back(ladderZId);
      intt_ladderPhiID[daughter_id].push_back(ladderPhiId);
      tpc_sectorID[daughter_id].push_back(sectorId);
      tpc_side[daughter_id].push_back(side);
    }
  }

  if (tpcseed)
  {
    for (auto cluster_iter = tpcseed->begin_cluster_keys(); cluster_iter != tpcseed->end_cluster_keys(); ++cluster_iter)
    {
      const auto &cluster_key = *cluster_iter;
      const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

      detector_layer[daughter_id].push_back(TrkrDefs::getLayer(cluster_key));
  
      unsigned int staveId, chipId, ladderZId, ladderPhiId, sectorId, side;
      staveId = chipId = ladderZId = ladderPhiId = sectorId = side = std::numeric_limits<unsigned int>::quiet_NaN();
  
      if (Use["TPC"] && trackerID == TrkrDefs::tpcId)
      {
        sectorId = TpcDefs::getSectorId(cluster_key);
        side = TpcDefs::getSide(cluster_key);
        ++detector_nHits_TPC[daughter_id];
      }
      else if (Use["TPOT"] && trackerID == TrkrDefs::micromegasId)
      {
        ++detector_nHits_TPOT[daughter_id];
      }
  
      mvtx_staveID[daughter_id].push_back(staveId);
      mvtx_chipID[daughter_id].push_back(chipId);
      intt_ladderZID[daughter_id].push_back(ladderZId);
      intt_ladderPhiID[daughter_id].push_back(ladderPhiId);
      tpc_sectorID[daughter_id].push_back(sectorId);
      tpc_side[daughter_id].push_back(side);
    }
  }

  for (auto state_iter = track->begin_states();
       state_iter != track->end_states();
       ++state_iter)
  {
    SvtxTrackState* tstate = state_iter->second;
    if (tstate->get_pathlength() != 0) //The first track state is an extrapolation so has no cluster
    {
      auto stateckey = tstate->get_cluskey();
      TrkrCluster *cluster = dst_clustermap->findCluster(stateckey);
      auto global = geometry->getGlobalPosition(stateckey, cluster);
  
      residual_x[daughter_id].push_back(global.x() - tstate->get_x());
      residual_y[daughter_id].push_back(global.y() - tstate->get_y());
      residual_z[daughter_id].push_back(global.z() - tstate->get_z());

      uint8_t id = TrkrDefs::getTrkrId(stateckey);
    
      switch (id)
      {
        case TrkrDefs::mvtxId:
          ++detector_nStates_MVTX[daughter_id];
          break;
        case TrkrDefs::inttId:
          ++detector_nStates_INTT[daughter_id];
          break;
        case TrkrDefs::tpcId:
          ++detector_nStates_TPC[daughter_id];
          break;
        case TrkrDefs::micromegasId:
          ++detector_nStates_TPOT[daughter_id];
          break;
        default:
         std::cout << "Cluster key doesnt match a tracking system, this shouldn't happen" << std::endl;
         break; 
      }
    }
  }
}

int KFParticle_truthAndDetTools::getPVID(PHCompositeNode *topNode, const KFParticle& kfpvertex)
{
  PHNodeIterator nodeIter(topNode);

  if (m_use_mbd_vertex_truth)
  {
    PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("MbdVertexMap"));
    if (findNode)
    {
      dst_mbdvertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
      MbdVertex* m_dst_vertex = dst_mbdvertexmap->get(kfpvertex.Id());
      return m_dst_vertex->get_beam_crossing();
    }
    else
    {
      std::cout << "KFParticle vertex matching: " << m_vtx_map_node_name_nTuple << " does not exist" << std::endl;
    }
  }
  else
  {
    PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_vtx_map_node_name_nTuple));
    if (findNode)
    {
      dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vtx_map_node_name_nTuple);
      SvtxVertex* m_dst_vertex = dst_vertexmap->get(kfpvertex.Id());
      return m_dst_vertex->get_beam_crossing();
    }
    else
    {
      std::cout << "KFParticle vertex matching: " << m_vtx_map_node_name_nTuple << " does not exist" << std::endl;
    }
  }

  return -100;
}

void KFParticle_truthAndDetTools::allPVInfo(PHCompositeNode *topNode,
                                            TTree * /*m_tree*/,
                                            const KFParticle &motherParticle,
                                            std::vector<KFParticle> daughters,
                                            std::vector<KFParticle> intermediates)
{
  KFParticle_Tools kfpTupleTools;
  std::vector<KFParticle> primaryVertices = kfpTupleTools.makeAllPrimaryVertices(topNode, m_vtx_map_node_name_nTuple);

  for (auto &primaryVertice : primaryVertices)
  {
    allPV_x.push_back(primaryVertice.GetX());
    allPV_y.push_back(primaryVertice.GetY());
    allPV_z.push_back(primaryVertice.GetZ());

    allPV_mother_IP.push_back(motherParticle.GetDistanceFromVertex(primaryVertice));
    allPV_mother_IPchi2.push_back(motherParticle.GetDeviationFromVertex(primaryVertice));

    for (unsigned int j = 0; j < daughters.size(); ++j)
    {
      allPV_daughter_IP[j].push_back(daughters[j].GetDistanceFromVertex(primaryVertice));
      allPV_daughter_IPchi2[j].push_back(daughters[j].GetDeviationFromVertex(primaryVertice));
    }

    for (unsigned int j = 0; j < intermediates.size(); ++j)
    {
      allPV_intermediates_IP[j].push_back(intermediates[j].GetDistanceFromVertex(primaryVertice));
      allPV_intermediates_IPchi2[j].push_back(intermediates[j].GetDeviationFromVertex(primaryVertice));
    }
  }
}

void KFParticle_truthAndDetTools::clearVectors()
{
  for (int i = 0; i < m_num_tracks_nTuple; ++i)
  {
    // Truth vectors
    m_true_daughter_track_history_PDG_ID[i].clear();
    m_true_daughter_track_history_PDG_mass[i].clear();
    m_true_daughter_track_history_px[i].clear();
    m_true_daughter_track_history_py[i].clear();
    m_true_daughter_track_history_pz[i].clear();
    m_true_daughter_track_history_pE[i].clear();
    m_true_daughter_track_history_pT[i].clear();

    // Detector vectors
    residual_x[i].clear();
    residual_y[i].clear();
    residual_z[i].clear();
    detector_layer[i].clear();
    mvtx_staveID[i].clear();
    mvtx_chipID[i].clear();
    intt_ladderZID[i].clear();
    intt_ladderPhiID[i].clear();
    tpc_sectorID[i].clear();
    tpc_side[i].clear();

    detector_nStates_MVTX[i] = 0;
    detector_nStates_INTT[i] = 0;
    detector_nStates_TPC[i] = 0;
    detector_nStates_TPOT[i] = 0;

    // PV vectors
    allPV_daughter_IP[i].clear();
    allPV_daughter_IPchi2[i].clear();
  }

  allPV_x.clear();
  allPV_y.clear();
  allPV_z.clear();
  allPV_z.clear();

  allPV_mother_IP.clear();
  allPV_mother_IPchi2.clear();

  for (int i = 0; i < m_num_intermediate_states_nTuple; ++i)
  {
    allPV_intermediates_IP[i].clear();
    allPV_intermediates_IPchi2[i].clear();
  }
}
